// EcopanGEM in-browser FBA / pFBA engine (ES module).
// Solves genome-scale COBRA-JSON models entirely client-side with glpk.js (WASM).
// Validated against COBRApy 0.27 (identical growth & pFBA sum|v|).
import GLPK from '../vendor/glpk.esm.js';

let _glpk = null;
export async function getGLPK() {
  if (!_glpk) _glpk = await GLPK();
  return _glpk;
}

// Build the LP: max c·v  s.t.  S·v = 0,  lb <= v <= ub.
// mediaBounds (optional): {EX_id: lower_bound}. When given, every exchange not
// listed is closed (lb=0), reproducing the "defined medium" convention.
// knockouts (optional): Set/array of reaction ids forced to lb=ub=0.
function buildLP(glpk, model, mediaBounds, knockouts) {
  const koSet = knockouts ? (knockouts instanceof Set ? knockouts : new Set(knockouts)) : null;
  const metRows = {};
  for (const met of model.metabolites)
    metRows[met.id] = { name: met.id, vars: [], bnds: { type: glpk.GLP_FX, ub: 0, lb: 0 } };

  const objVars = [];
  const bounds = [];
  for (const rxn of model.reactions) {
    let lb = rxn.lower_bound, ub = rxn.upper_bound;
    if (mediaBounds && rxn.id.startsWith('EX_'))
      lb = Object.prototype.hasOwnProperty.call(mediaBounds, rxn.id) ? mediaBounds[rxn.id] : 0;
    if (koSet && koSet.has(rxn.id)) { lb = 0; ub = 0; }
    let type;
    if (lb === ub) type = glpk.GLP_FX;
    else if (lb <= -1e30 && ub >= 1e30) type = glpk.GLP_FR;
    else if (lb <= -1e30) type = glpk.GLP_UP;
    else if (ub >= 1e30) type = glpk.GLP_LO;
    else type = glpk.GLP_DB;
    bounds.push({ name: rxn.id, type, ub, lb });
    for (const [met, coef] of Object.entries(rxn.metabolites))
      if (metRows[met]) metRows[met].vars.push({ name: rxn.id, coef });
    const oc = rxn.objective_coefficient || 0;
    if (oc !== 0) objVars.push({ name: rxn.id, coef: oc });
  }
  return {
    name: 'fba',
    objective: { direction: glpk.GLP_MAX, name: 'growth', vars: objVars },
    subjectTo: Object.values(metRows),
    bounds,
  };
}

export async function runFBA(model, mediaBounds, opts) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, model, mediaBounds, opts && opts.knockouts);
  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result;
  return {
    status: r.status, optimal: r.status === glpk.GLP_OPT,
    growth: r.z, objectiveId: (lp.objective.vars[0] || {}).name, fluxes: r.vars || {},
  };
}

/* Parsimonious FBA: fix biomass at the FBA optimum, then minimize total flux sum|v|.
   |v_i| is linearized with aux vars a_i >= v_i, a_i >= -v_i (min sum a_i).

   This runs at ZERO growth too, and that matters. When the objective is zero every
   feasible flux vector is optimal, so plain FBA returns an arbitrary one: whichever
   vertex the simplex happens to stop at, loops and all. (Measured on a lethal
   knockout: glpk.js returns loops pinned at the +/-1000 bounds, while COBRApy's
   solver returns a tame vector. Both are "optimal". Neither is an answer.) Minimising
   total flux picks one specific, reproducible state out of that space: the cell does
   as little as its bounds allow, which is what a dead cell does. Bailing out here and
   handing back the raw FBA vector, as this used to, meant every lethal knockout was
   reported with an arbitrary flux distribution. */
export async function runPFBA(model, mediaBounds, fbaResult, opts) {
  const glpk = await getGLPK();
  const fba = fbaResult || (await runFBA(model, mediaBounds, opts));
  if (!fba.optimal) return { ...fba, pfba: false };

  const lp = buildLP(glpk, model, mediaBounds, opts && opts.knockouts);
  const objId = fba.objectiveId;
  const gTarget = fba.growth > 1e-9 ? fba.growth * (1 - 1e-6) : 0;
  for (const b of lp.bounds)
    if (b.name === objId) { b.lb = gTarget; b.ub = Math.max(b.ub, fba.growth); b.type = glpk.GLP_DB; }

  const absVars = [], extraCons = [];
  for (const rxn of model.reactions) {
    const a = 'a_' + rxn.id;
    absVars.push({ name: a, coef: 1 });
    lp.bounds.push({ name: a, type: glpk.GLP_LO, lb: 0, ub: 1e30 });
    extraCons.push({ name: 'ap_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: -1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
    extraCons.push({ name: 'an_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: 1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
  }
  lp.subjectTo = lp.subjectTo.concat(extraCons);
  lp.objective = { direction: glpk.GLP_MIN, name: 'total_flux', vars: absVars };

  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result;
  const all = r.vars || {};
  const fluxes = {};
  for (const rxn of model.reactions) fluxes[rxn.id] = all[rxn.id] || 0;
  return {
    status: r.status, optimal: r.status === glpk.GLP_OPT,
    growth: fluxes[objId], objectiveId: objId, totalFlux: r.z, fluxes, pfba: true,
  };
}

// Flux Variability Analysis: min & max flux for each reaction in reactionIds,
// subject to biomass >= fraction * optimum. opts: {fraction, knockouts, onProgress}.
export async function runFVA(model, mediaBounds, reactionIds, opts = {}) {
  const glpk = await getGLPK();
  const fraction = opts.fraction != null ? opts.fraction : 1.0;
  const lp = buildLP(glpk, model, mediaBounds, opts.knockouts);
  const objId = (lp.objective.vars[0] || {}).name;
  const fba = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
  if (fba.status !== glpk.GLP_OPT || !(fba.z > 1e-9)) return { optimal: false, z: fba.z || 0, ranges: {} };
  // fix biomass >= fraction * optimum
  for (const b of lp.bounds) if (b.name === objId) { b.lb = fraction * fba.z; b.ub = Math.max(b.ub, fba.z); b.type = glpk.GLP_DB; }
  const ranges = {};
  let done = 0;
  for (const rid of reactionIds) {
    lp.objective = { direction: glpk.GLP_MIN, name: 'fva', vars: [{ name: rid, coef: 1 }] };
    const mn = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z;
    lp.objective.direction = glpk.GLP_MAX;
    const mx = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z;
    ranges[rid] = { min: mn, max: mx };
    if (opts.onProgress) opts.onProgress(++done, reactionIds.length);
  }
  return { optimal: true, z: fba.z, fraction, ranges };
}

// Dynamic (batch) FBA by static optimization. Substrate uptake follows
// Michaelis-Menten; biomass and tracked metabolite concentrations integrate by
// forward Euler. opts: {substrateEx, substrate0, biomass0, vmax, km, dt, tmax, trackEx, knockouts, onProgress}.
export async function runDFBA(model, mediaBounds, opts = {}) {
  const o = Object.assign({ substrateEx: null, substrate0: 10, biomass0: 0.01,
    vmax: 10, km: 0.5, dt: 0.1, tmax: 15, trackEx: [] }, opts);
  // never hardcode a carbon source: Lactobacillaceae models have no EX_glc__D_e
  if (!o.substrateEx) o.substrateEx = defaultSubstrate(model);
  const media = { ...mediaBounds };
  const conc = { [o.substrateEx]: o.substrate0 };
  for (const e of o.trackEx) if (!(e in conc)) conc[e] = 0;
  const series = { t: [], biomass: [], conc: {}, substrateEx: o.substrateEx };
  for (const e of Object.keys(conc)) series.conc[e] = [];
  let X = o.biomass0, t = 0;
  const nsteps = Math.ceil(o.tmax / o.dt);
  for (let i = 0; i <= nsteps; i++) {
    const S = Math.max(0, conc[o.substrateEx]);
    let vUp = o.vmax * S / (o.km + S);
    if (X > 0 && o.dt > 0) vUp = Math.min(vUp, S / (X * o.dt)); // never consume more than present
    media[o.substrateEx] = -vUp;
    const fba = await runFBA(model, media, { knockouts: o.knockouts });
    const mu = fba.optimal ? fba.growth : 0;
    series.t.push(+t.toFixed(4)); series.biomass.push(X);
    for (const e of Object.keys(conc)) series.conc[e].push(Math.max(0, conc[e]));
    for (const e of Object.keys(conc)) {
      const flux = fba.fluxes[e] || 0;
      conc[e] += flux * X * o.dt; if (conc[e] < 0) conc[e] = 0;
    }
    X += mu * X * o.dt; t += o.dt;
    if (opts.onProgress) opts.onProgress(i, nsteps);
    if (mu <= 1e-9 && S <= 1e-9) break;
  }
  return series;
}

// Production envelope: trade-off between biomass and a target product. Fix the
// product flux at a grid of values from 0 to its maximum and record the max &
// min biomass at each. opts: {points, knockouts}.
export async function productionEnvelope(model, media, productId, opts = {}) {
  const glpk = await getGLPK();
  const points = opts.points || 20;
  const lpP = buildLP(glpk, model, media, opts.knockouts);
  lpP.objective = { direction: glpk.GLP_MAX, name: 'p', vars: [{ name: productId, coef: 1 }] };
  const prodMax = (await glpk.solve(lpP, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z || 0;
  const out = [];
  for (let i = 0; i < points; i++) {
    const v = prodMax * i / (points - 1);
    const lp = buildLP(glpk, model, media, opts.knockouts);
    for (const b of lp.bounds) if (b.name === productId) { b.lb = v; b.ub = v; b.type = glpk.GLP_FX; }
    const mx = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    lp.objective = { ...lp.objective, direction: glpk.GLP_MIN };
    const mn = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    out.push({ product: v, growthMax: mx.status === glpk.GLP_OPT ? mx.z : 0, growthMin: mn.status === glpk.GLP_OPT ? mn.z : 0 });
    if (opts.onProgress) opts.onProgress(i + 1, points);
  }
  return { productId, prodMax, points: out };
}

// Phenotype phase plane: max biomass over a grid of two uptake capacities.
// opts: {n, xMax, yMax, knockouts, onProgress}. Returns {xs, ys, Z} (Z[j][i]).
export async function phasePlane(model, media, xId, yId, opts = {}) {
  const n = opts.n || 20;
  const xMax = opts.xMax != null ? opts.xMax : Math.abs(media[xId] != null ? media[xId] : 20);
  const yMax = opts.yMax != null ? opts.yMax : Math.abs(media[yId] != null ? media[yId] : 20);
  const xs = [], ys = [], Z = [];
  for (let i = 0; i < n; i++) xs.push(xMax * i / (n - 1));
  for (let j = 0; j < n; j++) ys.push(yMax * j / (n - 1));
  let done = 0;
  for (let j = 0; j < n; j++) {
    const row = [];
    for (let i = 0; i < n; i++) {
      const m2 = { ...media, [xId]: -xs[i], [yId]: -ys[j] };
      const fba = await runFBA(model, m2, { knockouts: opts.knockouts });
      row.push(fba.optimal ? fba.growth : 0);
      if (opts.onProgress) opts.onProgress(++done, n * n);
    }
    Z.push(row);
  }
  return { xId, yId, xs, ys, Z };
}

// Single-reaction deletion (essentiality) over a list of reactions.
// Returns [{id, name, subsystem, growth, ratio}]. opts: {knockouts, onProgress}.
export async function essentialityScan(model, mediaBounds, reactionIds, opts = {}) {
  const wt = await runFBA(model, mediaBounds, { knockouts: opts.knockouts });
  const wtG = wt.optimal ? wt.growth : 0;
  const nameById = {}, subById = {};
  model.reactions.forEach(r => { nameById[r.id] = r.name || ''; subById[r.id] = r.subsystem || ''; });
  const base = opts.knockouts ? [...opts.knockouts] : [];
  const out = [];
  let done = 0;
  for (const rid of reactionIds) {
    const fba = await runFBA(model, mediaBounds, { knockouts: base.concat(rid) });
    const g = fba.optimal ? fba.growth : 0;
    out.push({ id: rid, name: nameById[rid], subsystem: subById[rid], growth: g, ratio: wtG > 1e-9 ? g / wtG : 0 });
    if (opts.onProgress) opts.onProgress(++done, reactionIds.length);
  }
  return { wtGrowth: wtG, results: out };
}

/* ═══════════════════════════════════════════════════════════════════════════
   COBRApy parity. glpk.js exposes a linear solver only, no integers and no
   quadratic term, so everything below is formulated as an LP. What that rules
   out is stated honestly in the UI rather than approximated: ROOM, exact
   loopless (add_loopless), gapfilling and OptKnock all need MILP, and
   quadratic MOMA needs QP.
   ═══════════════════════════════════════════════════════════════════════════ */

/* ── Gene-protein-reaction rules ──────────────────────────────────────────────
   A gene knockout is not a reaction knockout. Genes act through the GPR: an
   isozyme (a or b) survives losing one gene, a complex (a and b) does not.
   Recursive-descent over the boolean rule, never eval(). */
export function evalGPR(rule, isOn) {
  if (!rule || !rule.trim()) return true;          // no GPR: spontaneous / orphan, always on
  const toks = rule.replace(/\(/g, ' ( ').replace(/\)/g, ' ) ').split(/\s+/).filter(Boolean);
  let i = 0;
  const atom = () => {
    if (toks[i] === '(') { i++; const v = or_(); if (toks[i] === ')') i++; return v; }
    const g = toks[i++];
    return isOn(g);
  };
  const and_ = () => {
    let v = atom();
    while (toks[i] && toks[i].toLowerCase() === 'and') { i++; const r = atom(); v = v && r; }
    return v;
  };
  const or_ = () => {
    let v = and_();
    while (toks[i] && toks[i].toLowerCase() === 'or') { i++; const r = and_(); v = v || r; }
    return v;
  };
  try { return or_(); } catch (e) { return true; }
}

/** Reactions switched off when these genes are deleted. */
export function reactionsOffForGenes(model, deadGenes) {
  const dead = new Set(deadGenes);
  const off = [];
  for (const r of model.reactions) {
    const rule = r.gene_reaction_rule || '';
    if (!rule.trim()) continue;                    // not under gene control
    if (!evalGPR(rule, g => !dead.has(g))) off.push(r.id);
  }
  return off;
}

export function listGenes(model) {
  return (model.genes || []).map(g => ({ id: g.id, name: g.name || '' }));
}

/** single_gene_deletion. Returns [{gene, name, nOff, growth, ratio, rxns}]. */
export async function singleGeneDeletion(model, media, geneIds, opts = {}) {
  const base = opts.knockouts ? [...opts.knockouts] : [];
  const wt = await runFBA(model, media, { knockouts: base });
  const wtG = wt.optimal ? wt.growth : 0;
  const nameById = {};
  (model.genes || []).forEach(g => { nameById[g.id] = g.name || ''; });
  const out = [];
  let done = 0;
  for (const gid of geneIds) {
    const off = reactionsOffForGenes(model, [gid]);
    const f = off.length ? await runFBA(model, media, { knockouts: base.concat(off) }) : wt;
    const g = f.optimal ? f.growth : 0;
    out.push({ gene: gid, name: nameById[gid] || '', nOff: off.length, rxns: off,
               growth: g, ratio: wtG > 1e-9 ? g / wtG : 0 });
    if (opts.onProgress) opts.onProgress(++done, geneIds.length);
  }
  return { wtGrowth: wtG, results: out };
}

/** double_reaction_deletion / double_gene_deletion. Pairwise over `ids`.
    Synthetic lethality: neither single kills, the pair does. */
export async function doubleDeletion(model, media, ids, kind = 'reaction', opts = {}) {
  const base = opts.knockouts ? [...opts.knockouts] : [];
  const wt = await runFBA(model, media, { knockouts: base });
  const wtG = wt.optimal ? wt.growth : 0;
  const toRxns = (id) => kind === 'gene' ? reactionsOffForGenes(model, [id]) : [id];

  const singles = {};
  for (const id of ids) {
    const f = await runFBA(model, media, { knockouts: base.concat(toRxns(id)) });
    singles[id] = wtG > 1e-9 ? (f.optimal ? f.growth : 0) / wtG : 0;
  }
  const pairs = [], n = ids.length;
  const total = (n * (n - 1)) / 2;
  let done = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const kos = base.concat(toRxns(ids[i]), toRxns(ids[j]));
      const f = await runFBA(model, media, { knockouts: kos });
      const r = wtG > 1e-9 ? (f.optimal ? f.growth : 0) / wtG : 0;
      // synthetic lethal: both singles tolerable, the pair is not
      const sl = r < 0.01 && singles[ids[i]] > 0.5 && singles[ids[j]] > 0.5;
      pairs.push({ a: ids[i], b: ids[j], ratio: r, ra: singles[ids[i]], rb: singles[ids[j]], synthetic: sl });
      if (opts.onProgress) opts.onProgress(++done, total);
    }
  }
  pairs.sort((x, y) => x.ratio - y.ratio);
  return { wtGrowth: wtG, kind, singles, pairs, ids };
}

/* ── Knockdown titration ──────────────────────────────────────────────────────
   A knockout is a single point on a curve. It answers "can the cell live without
   this reaction" but not "how much of it did the cell actually need", and those
   are different questions: a reaction can carry ten times the flux required of it
   and show no phenotype until the last tenth is taken away.

   Here the targeted reactions keep a fraction f of the throughput they carried in
   the wild type, and f is swept from 0 (the knockout) to 1 (the wild type). A
   cliff means there is no slack. A long plateau that falls only at the end means
   the reaction was carrying far more flux than the cell needed, and a partial
   inhibitor would do nothing.

   Growth is monotone non-decreasing in f, because raising f only relaxes bounds:
   that is a free correctness check on the result.
   opts: {wtFluxes (required), steps, onProgress} */
export async function knockdownCurve(model, media, rxnIds, opts = {}) {
  const glpk = await getGLPK();
  const wt = opts.wtFluxes || {};
  const ids = rxnIds.filter(id => Math.abs(wt[id] || 0) > 1e-6);
  if (!ids.length) return { optimal: false, reason: 'no-wt-flux', ids: [], points: [] };

  const lp = buildLP(glpk, model, media, null);   // caps are set here, not by knockout
  const byName = new Map(lp.bounds.map(b => [b.name, b]));
  const orig = new Map(ids.map(id => {
    const b = byName.get(id);
    return [id, b ? { lb: b.lb, ub: b.ub } : null];
  }));

  const steps = Math.max(2, opts.steps || 21);
  const points = [];
  for (let i = 0; i < steps; i++) {
    const f = i / (steps - 1);
    for (const id of ids) {
      const b = byName.get(id), o = orig.get(id);
      if (!b || !o) continue;
      const c = f * Math.abs(wt[id]);                 // allowed magnitude at this step
      b.lb = Math.max(o.lb, -c);                      // never widen past the model's own bounds
      b.ub = Math.min(o.ub, c);
      b.type = (b.lb === b.ub) ? glpk.GLP_FX : glpk.GLP_DB;
    }
    const r = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    points.push({ f, growth: (r.status === glpk.GLP_OPT && r.z > 0) ? r.z : 0 });
    if (opts.onProgress) opts.onProgress(i + 1, steps);
  }
  return { optimal: true, ids, points };
}

/* ── Linear MOMA ──────────────────────────────────────────────────────────────
   After a knockout a cell does not re-optimise growth, it stays as close to its
   old flux state as it can. MOMA minimises the distance to the wild-type flux.
   COBRApy's default MOMA is quadratic; glpk.js has no QP, so this is the linear
   form (cobra.flux_analysis.moma(..., linear=True)), which is the one usually
   used for large models anyway. */
export async function runLMOMA(model, media, refFluxes, opts = {}) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, model, media, opts.knockouts);
  const objId = (lp.objective.vars[0] || {}).name;
  const absVars = [], extra = [];
  for (const rxn of model.reactions) {
    const w = refFluxes[rxn.id] || 0;
    const a = 'd_' + rxn.id;
    absVars.push({ name: a, coef: 1 });
    lp.bounds.push({ name: a, type: glpk.GLP_LO, lb: 0, ub: 1e30 });
    // a >= v - w   and   a >= w - v
    extra.push({ name: 'dp_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: -1 }],
                 bnds: { type: glpk.GLP_LO, lb: -w, ub: 0 } });
    extra.push({ name: 'dn_' + rxn.id, vars: [{ name: a, coef: 1 }, { name: rxn.id, coef: 1 }],
                 bnds: { type: glpk.GLP_LO, lb: w, ub: 0 } });
  }
  lp.subjectTo = lp.subjectTo.concat(extra);
  lp.objective = { direction: glpk.GLP_MIN, name: 'moma', vars: absVars };
  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result, all = r.vars || {};
  const fluxes = {};
  for (const rxn of model.reactions) fluxes[rxn.id] = all[rxn.id] || 0;
  return { status: r.status, optimal: r.status === glpk.GLP_OPT, objectiveId: objId,
           growth: fluxes[objId] || 0, distance: r.z, fluxes, moma: true };
}

/* ── CycleFreeFlux (loopless_solution) ────────────────────────────────────────
   Thermodynamically infeasible internal loops carry flux that no cell could.
   Exact loop removal needs MILP; CycleFreeFlux (Desouki 2015) is the LP form:
   hold the objective and every exchange at the given solution, keep each
   internal reaction's sign, and minimise total internal flux. Same growth, same
   exchanges, no loops. */
export async function looplessSolution(model, media, fluxes, opts = {}) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, model, media, opts.knockouts);
  const objId = (lp.objective.vars[0] || {}).name;
  const minVars = [];
  for (const b of lp.bounds) {
    const v0 = fluxes[b.name] || 0;
    if (b.name.startsWith('EX_') || b.name === objId) {   // pin what the cell exchanges and how fast it grows
      b.lb = v0; b.ub = v0; b.type = glpk.GLP_FX;
      continue;
    }
    if (v0 >= 0) { b.lb = 0; b.ub = Math.max(v0, 0); }     // same direction, may shrink to zero
    else { b.lb = Math.min(v0, 0); b.ub = 0; }
    b.type = (b.lb === b.ub) ? glpk.GLP_FX : glpk.GLP_DB;
    minVars.push({ name: b.name, coef: v0 >= 0 ? 1 : -1 }); // |v| is linear once the sign is fixed
  }
  lp.objective = { direction: glpk.GLP_MIN, name: 'cyclefree', vars: minVars };
  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result, all = r.vars || {};
  const out = {}, removed = [];
  for (const rxn of model.reactions) {
    out[rxn.id] = all[rxn.id] || 0;
    const before = fluxes[rxn.id] || 0;
    if (Math.abs(before) > 1e-6 && Math.abs(out[rxn.id]) < 1e-9) removed.push({ id: rxn.id, was: before });
  }
  removed.sort((a, b) => Math.abs(b.was) - Math.abs(a.was));
  return { optimal: r.status === glpk.GLP_OPT, fluxes: out, removed,
           totalBefore: sumAbs(fluxes, model), totalAfter: sumAbs(out, model) };
}
const sumAbs = (f, model) => model.reactions.reduce((s, r) => s + Math.abs(f[r.id] || 0), 0);

/* ── find_blocked_reactions + model QC ───────────────────────────────────────
   Blocked reactions can carry no flux at all under this medium: dead code in
   the reconstruction. Mass and charge balance are checked without solving. */
export async function findBlockedReactions(model, media, opts = {}) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, model, media, opts.knockouts);
  const ids = (opts.reactionIds || model.reactions.map(r => r.id));
  const blocked = [];
  let done = 0;
  for (const rid of ids) {
    lp.objective = { direction: glpk.GLP_MAX, name: 'b', vars: [{ name: rid, coef: 1 }] };
    const mx = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z || 0;
    lp.objective.direction = glpk.GLP_MIN;
    const mn = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z || 0;
    if (Math.abs(mx) < 1e-9 && Math.abs(mn) < 1e-9) blocked.push(rid);
    if (opts.onProgress) opts.onProgress(++done, ids.length);
  }
  return { blocked, tested: ids.length };
}

/** Structural QC. No LP: pure stoichiometry, the way MEMOTE does it. */
export function modelQC(model) {
  const metById = {};
  for (const m of model.metabolites) metById[m.id] = m;
  const parseF = (f) => {
    const c = {};
    if (!f) return null;
    const re = /([A-Z][a-z]?)(\d*\.?\d*)/g; let m2, seen = false;
    while ((m2 = re.exec(f))) { seen = true; c[m2[1]] = (c[m2[1]] || 0) + (m2[2] ? parseFloat(m2[2]) : 1); }
    return seen ? c : null;
  };
  const massImbalanced = [], chargeImbalanced = [];
  for (const r of model.reactions) {
    if (r.id.startsWith('EX_') || /BIOMASS|biomass/i.test(r.id)) continue;  // exchanges & biomass are open by design
    const tot = {}; let charge = 0, ok = true;
    for (const [mid, coef] of Object.entries(r.metabolites || {})) {
      const met = metById[mid]; if (!met) { ok = false; break; }
      const f = parseF(met.formula); if (!f) { ok = false; break; }
      for (const [el, n] of Object.entries(f)) tot[el] = (tot[el] || 0) + n * coef;
      charge += (met.charge || 0) * coef;
    }
    if (!ok) continue;
    const off = Object.entries(tot).filter(([, v]) => Math.abs(v) > 1e-6);
    if (off.length) massImbalanced.push({ id: r.id, name: r.name || '', delta: Object.fromEntries(off) });
    if (Math.abs(charge) > 1e-6) chargeImbalanced.push({ id: r.id, name: r.name || '', charge });
  }
  // metabolites only ever produced, or only ever consumed
  const prod = {}, cons = {};
  for (const r of model.reactions) {
    for (const [mid, coef] of Object.entries(r.metabolites || {})) {
      const rev = r.lower_bound < 0, fwd = r.upper_bound > 0;
      if ((coef > 0 && fwd) || (coef < 0 && rev)) prod[mid] = (prod[mid] || 0) + 1;
      if ((coef < 0 && fwd) || (coef > 0 && rev)) cons[mid] = (cons[mid] || 0) + 1;
    }
  }
  const deadEnds = model.metabolites.filter(m => !prod[m.id] || !cons[m.id])
    .map(m => ({ id: m.id, name: m.name || '', onlyProduced: !cons[m.id], onlyConsumed: !prod[m.id] }));
  const noFormula = model.metabolites.filter(m => !m.formula).map(m => m.id);
  const usedGenes = new Set();
  model.reactions.forEach(r => (r.gene_reaction_rule || '').split(/[^\w.\-]+/).forEach(g => g && usedGenes.add(g)));
  const orphanGenes = (model.genes || []).filter(g => !usedGenes.has(g.id)).map(g => g.id);
  const noGPR = model.reactions.filter(r => !r.id.startsWith('EX_') && !(r.gene_reaction_rule || '').trim())
    .map(r => r.id);
  return { massImbalanced, chargeImbalanced, deadEnds, noFormula, orphanGenes, noGPR,
           nReactions: model.reactions.length, nMetabolites: model.metabolites.length, nGenes: (model.genes || []).length };
}

/* ── FSEOF ────────────────────────────────────────────────────────────────────
   Flux Scanning with Enforced Objective Flux (Choi 2010). Force the product
   secretion up in steps while still maximising growth, and watch which internal
   reactions have to carry more flux. Those are the amplification targets: the
   reactions a strain engineer should over-express. */
export async function fseof(model, media, productEx, opts = {}) {
  const steps = opts.steps || 10;
  const glpk = await getGLPK();
  const lpP = buildLP(glpk, model, media, opts.knockouts);
  lpP.objective = { direction: glpk.GLP_MAX, name: 'p', vars: [{ name: productEx, coef: 1 }] };
  const prodMax = (await glpk.solve(lpP, { msglev: glpk.GLP_MSG_OFF, presol: true })).result.z || 0;
  if (prodMax <= 1e-9) return { productEx, prodMax: 0, levels: [], targets: [] };

  const wt = await runPFBA(model, media, null, { knockouts: opts.knockouts });
  const levels = [], series = {};
  for (const r of model.reactions) series[r.id] = [];
  let done = 0;
  for (let i = 1; i <= steps; i++) {
    const enforce = prodMax * (i / (steps + 1));
    const lp = buildLP(glpk, model, media, opts.knockouts);
    for (const b of lp.bounds) if (b.name === productEx) { b.lb = enforce; b.ub = 1e30; b.type = glpk.GLP_LO; }
    const fba = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    if (fba.status !== glpk.GLP_OPT) { if (opts.onProgress) opts.onProgress(++done, steps); continue; }
    const f = fba.vars || {};
    levels.push({ enforced: enforce, growth: fba.z });
    for (const r of model.reactions) series[r.id].push(f[r.id] || 0);
    if (opts.onProgress) opts.onProgress(++done, steps);
  }
  // a target increases monotonically with the enforced product flux
  const targets = [];
  for (const r of model.reactions) {
    if (r.id === productEx || r.id.startsWith('EX_')) continue;
    const v = series[r.id];
    if (v.length < 3) continue;
    const first = v[0], last = v[v.length - 1];
    if (Math.abs(last) < 1e-6) continue;
    let up = 0;
    for (let i = 1; i < v.length; i++) if (Math.abs(v[i]) >= Math.abs(v[i - 1]) - 1e-9) up++;
    const mono = up / (v.length - 1);
    const fold = Math.abs(first) > 1e-6 ? Math.abs(last) / Math.abs(first) : Infinity;
    if (mono > 0.85 && Math.abs(last) > Math.abs(first) + 1e-6) {
      targets.push({ id: r.id, name: r.name || '', subsystem: r.subsystem || '',
                     wt: wt.fluxes ? (wt.fluxes[r.id] || 0) : 0, first, last, fold, series: v });
    }
  }
  targets.sort((a, b) => (Math.abs(b.last) - Math.abs(b.first)) - (Math.abs(a.last) - Math.abs(a.first)));
  return { productEx, prodMax, levels, targets };
}

/* ── Flux sampling (ACHR) ────────────────────────────────────────────────────
   Artificial Centering Hit-and-Run. Warm-up points come from optimising random
   objective directions; then each step moves from a point toward the running
   centre along a random chord of the polytope. Same family as COBRApy's OptGP.
   Approximate: with a few hundred warm-up LPs it explores the space, it is not
   a proof of uniformity. */
export async function sampleFluxes(model, media, opts = {}) {
  const glpk = await getGLPK();
  const nWarm = opts.warmup || 120;
  const nSamp = opts.samples || 400;
  const track = opts.trackIds || model.reactions.filter(r => r.id.startsWith('EX_')).map(r => r.id).slice(0, 30);
  const base = buildLP(glpk, model, media, opts.knockouts);
  const objId = (base.objective.vars[0] || {}).name;

  // hold growth at a fraction of optimum so we sample the space a living cell occupies
  const opt = (await glpk.solve(base, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
  if (opt.status !== glpk.GLP_OPT || !(opt.z > 1e-9)) return { ok: false, reason: 'no growth on this medium' };
  const frac = opts.fraction != null ? opts.fraction : 0.9;
  for (const b of base.bounds) if (b.name === objId) { b.lb = frac * opt.z; b.ub = Math.max(b.ub, opt.z); b.type = glpk.GLP_DB; }

  const rxns = model.reactions.map(r => r.id);
  const warm = [];
  let done = 0;
  for (let i = 0; i < nWarm; i++) {
    const lp = JSON.parse(JSON.stringify(base));
    lp.objective = { direction: Math.random() < 0.5 ? glpk.GLP_MAX : glpk.GLP_MIN, name: 'w',
                     vars: rxns.filter(() => Math.random() < 0.15).map(id => ({ name: id, coef: Math.random() * 2 - 1 })) };
    if (!lp.objective.vars.length) lp.objective.vars = [{ name: rxns[i % rxns.length], coef: 1 }];
    const r = (await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true })).result;
    if (r.status === glpk.GLP_OPT) warm.push(r.vars || {});
    if (opts.onProgress) opts.onProgress(++done, nWarm + nSamp);
  }
  if (warm.length < 5) return { ok: false, reason: 'warm-up failed' };

  // hit-and-run between warm-up points, recentring as we go
  const centre = {};
  for (const id of rxns) centre[id] = warm.reduce((s, p) => s + (p[id] || 0), 0) / warm.length;
  const lo = {}, hi = {};
  for (const b of base.bounds) { lo[b.name] = b.lb; hi[b.name] = b.ub; }
  /* Only reactions that can actually move matter. A defined medium fixes every
     closed exchange at lb = ub = 0; with a whisker of numerical noise in the
     direction, such a variable collapses the allowed step to [0, 0] and the walk
     stalls. Free variables only, and a real tolerance on the direction. */
  const free = rxns.filter(id => (hi[id] - lo[id]) > 1e-7);
  const samples = [];
  const clamp = (pt) => { for (const id of free) pt[id] = Math.min(hi[id], Math.max(lo[id], pt[id] || 0)); return pt; };
  let cur = clamp({ ...warm[0] });
  let stalls = 0;
  for (let s = 0; s < nSamp; s++) {
    let t = null, dir = null, norm = 1;
    for (let attempt = 0; attempt < 6 && t === null; attempt++) {
      const p = warm[Math.floor(Math.random() * warm.length)];
      dir = {}; norm = 0;
      for (const id of free) { const d = (p[id] || 0) - centre[id]; dir[id] = d; norm += d * d; }
      norm = Math.sqrt(norm);
      if (norm < 1e-9) continue;
      let tmax = Infinity, tmin = -Infinity;
      for (const id of free) {
        const d = dir[id] / norm;
        if (Math.abs(d) < 1e-9) continue;
        const up = (hi[id] - (cur[id] || 0)) / d;
        const dn = (lo[id] - (cur[id] || 0)) / d;
        const mn = Math.min(up, dn), mx = Math.max(up, dn);
        if (mx < tmax) tmax = mx;
        if (mn > tmin) tmin = mn;
      }
      if (!isFinite(tmax) || !isFinite(tmin) || (tmax - tmin) < 1e-9) continue;
      // stay off the exact face, it is where the numerics go bad
      t = tmin + (0.02 + 0.96 * Math.random()) * (tmax - tmin);
    }
    if (t === null) { cur = clamp({ ...warm[Math.floor(Math.random() * warm.length)] }); stalls++; if (opts.onProgress) opts.onProgress(++done, nWarm + nSamp); continue; }

    const next = { ...cur };
    for (const id of free) next[id] = (cur[id] || 0) + t * dir[id] / norm;
    cur = clamp(next);
    const row = {}; for (const id of track) row[id] = cur[id] || 0;
    samples.push(row);
    for (const id of free) centre[id] = (centre[id] * (warm.length + s) + cur[id]) / (warm.length + s + 1);
    if (opts.onProgress) opts.onProgress(++done, nWarm + nSamp);
  }
  return { ok: true, track, samples, warmup: warm.length, fraction: frac, free: free.length, stalls };
}

/** metabolite.summary(): who makes it and who eats it, at this flux state. */
export function metaboliteSummary(model, fluxes, metId, tol = 1e-8) {
  const producing = [], consuming = [];
  for (const r of model.reactions) {
    const coef = (r.metabolites || {})[metId];
    if (!coef) continue;
    const v = fluxes[r.id] || 0;
    if (Math.abs(v) < tol) continue;
    const rate = coef * v;                          // >0 net produced, <0 net consumed
    const rec = { id: r.id, name: r.name || '', flux: v, rate: Math.abs(rate),
                  subsystem: r.subsystem || '' };
    if (rate > 0) producing.push(rec); else consuming.push(rec);
  }
  producing.sort((a, b) => b.rate - a.rate);
  consuming.sort((a, b) => b.rate - a.rate);
  const total = producing.reduce((s, x) => s + x.rate, 0);
  producing.forEach(x => x.share = total ? x.rate / total : 0);
  const totalC = consuming.reduce((s, x) => s + x.rate, 0);
  consuming.forEach(x => x.share = totalC ? x.rate / totalC : 0);
  return { metId, producing, consuming, turnover: total };
}

/* ── Exchange-id resolution across BiGG naming generations ────────────────────
   EcopanGEM follows current BiGG and writes stereo tags with a DOUBLE underscore
   (EX_glc__D_e). LactoPanGEM predates that change and writes a single one
   (EX_glc_D_e). buildLP closes every exchange a medium does not name, so an
   unresolved compound is not a cosmetic miss: it silently deletes a nutrient.
   Measured on a Lactobacillaceae GEM with BHI medium: exact-match binds 32/56
   compounds and the model does not grow at all (0.000 h-1); resolving binds
   53/56 and it grows at 0.210 h-1. Both verified against COBRApy 0.27.
   Collapsing both spellings to one key makes the lookup work in either
   direction, so a medium written in either generation binds to either model. */
const canonEx = (id) => String(id).replace(/__([LDRS])_e$/, '_$1_e').toLowerCase();

export function exchangeIndex(model) {
  const idx = new Map();
  for (const r of model.reactions) {
    if (!r.id.startsWith('EX_')) continue;
    idx.set(r.id, r.id);
    const k = canonEx(r.id);
    if (!idx.has(k)) idx.set(k, r.id);
  }
  return idx;
}

export function resolveExchange(idx, wanted) {
  if (!wanted) return null;
  return idx.get(wanted) || idx.get(canonEx(wanted)) || null;
}

// Inorganic ions and water. Sparse media (a defined carbon source, a food record)
// often omit these, and without them nothing grows. Offered as an explicit,
// visible option rather than applied silently.
export const MINERALS = ['EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_so4_e', 'EX_nh4_e', 'EX_k_e',
  'EX_na1_e', 'EX_mg2_e', 'EX_ca2_e', 'EX_fe2_e', 'EX_fe3_e', 'EX_cl_e', 'EX_cu2_e',
  'EX_mn2_e', 'EX_zn2_e', 'EX_cobalt2_e', 'EX_mobd_e', 'EX_ni2_e', 'EX_sel_e',
  'EX_slnt_e', 'EX_tungs_e', 'EX_cbl1_e'];

/* Bind a medium to ONE model and report exactly what landed.
   components: [{exchange, lower_bound, name?}]
   opts: {openMinerals:bool, o2:null|number}
   returns {bounds, mapped[], missing[], added[], coverage} */
export function bindMedium(model, components, opts = {}) {
  const idx = exchangeIndex(model);
  const bounds = {}, mapped = [], missing = [], added = [];
  for (const c of components || []) {
    const rid = resolveExchange(idx, c.exchange);
    if (rid) {
      bounds[rid] = Number(c.lower_bound);
      mapped.push({ ...c, resolved: rid, renamed: rid !== c.exchange });
    } else {
      missing.push(c);
    }
  }
  if (opts.openMinerals) {
    for (const w of MINERALS) {
      const rid = resolveExchange(idx, w);
      if (rid && !(rid in bounds)) { bounds[rid] = -1000; added.push(rid); }
    }
  }
  if (opts.o2 != null) {
    const rid = resolveExchange(idx, 'EX_o2_e');
    if (rid) bounds[rid] = Number(opts.o2);
  }
  const n = (components || []).length;
  return { bounds, mapped, missing, added, coverage: n ? mapped.length / n : 1 };
}

// A carbon source that actually exists in THIS model. dFBA used to hardcode
// EX_glc__D_e, which no Lactobacillaceae model has.
export function defaultSubstrate(model) {
  const idx = exchangeIndex(model);
  for (const w of ['EX_glc__D_e', 'EX_glc_e', 'EX_fru_e', 'EX_sucr_e', 'EX_lcts_e']) {
    const rid = resolveExchange(idx, w);
    if (rid) return rid;
  }
  const ex = model.reactions.filter(r => r.id.startsWith('EX_'));
  return ex.length ? ex[0].id : null;
}

// All exchange reactions in a model: {id, name, met, defaultLb, defaultUb}. For the media editor.
export function listExchanges(model) {
  return model.reactions
    .filter(r => r.id.startsWith('EX_'))
    .map(r => ({ id: r.id, name: r.name || r.id, met: Object.keys(r.metabolites || {})[0] || '',
                 defaultLb: r.lower_bound, defaultUb: r.upper_bound }))
    .sort((a, b) => a.id.localeCompare(b.id));
}

// All reactions (id + name) for the knockout search.
export function listReactions(model) {
  return model.reactions.map(r => ({ id: r.id, name: r.name || '', subsystem: r.subsystem || '' }));
}

// Split exchange fluxes into uptake (negative) and secretion (positive), sorted by magnitude.
export function exchangeReport(model, fluxes, tol = 1e-6) {
  const nameById = {};
  for (const r of model.reactions) nameById[r.id] = r.name || r.id;
  const uptake = [], secretion = [];
  for (const [id, v] of Object.entries(fluxes)) {
    if (!id.startsWith('EX_')) continue;
    if (v < -tol) uptake.push({ id, name: nameById[id], flux: v });
    else if (v > tol) secretion.push({ id, name: nameById[id], flux: v });
  }
  uptake.sort((a, b) => a.flux - b.flux);
  secretion.sort((a, b) => b.flux - a.flux);
  return { uptake, secretion };
}
