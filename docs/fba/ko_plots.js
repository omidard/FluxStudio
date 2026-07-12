/* ── Knockout vs wild type: the figures ───────────────────────────────────────
   Five plots, each answering a different question, because "what did the
   knockout do" is really five questions:

     1  Titration      how much of the reaction did the cell actually need?
     2  Pathway impact where in metabolism did the damage land?
     3  Re-routing     what did the whole flux distribution do?
     4  Largest moves  which reactions specifically, by name?
     5  Boundary shift what does the cell now eat and excrete?

   ⚠ THE TRAP THIS MODULE EXISTS TO AVOID. A flux-difference plot is a lie if the
   fluxes came from plain FBA. FBA has alternate optima: many different flux
   vectors give exactly the same growth rate, and the solver returns an arbitrary
   one. Diff two arbitrary choices and you get a picture full of "re-routing" that
   is solver noise, not biology.

   So every change shown here is checked against the alternate optima. For each
   reaction we compute its full feasible flux range at the wild-type optimum and
   at the knockout optimum (that is FVA, run on both). If the two ranges are
   DISJOINT, no choice of solution could have given the same flux in both, so the
   knockout forces the change: it is real. If they OVERLAP, the difference we
   happened to plot could have been the solver picking differently, and the plot
   says so rather than pretending.

   The test is conservative on purpose. FVA explores the full FBA optimal space,
   which is wider than the space pFBA or MOMA actually picks from, so a change can
   be genuine under parsimony and still land in "not distinguishable". We would
   rather under-claim than draw a confident arrow through nothing.
   ────────────────────────────────────────────────────────────────────────── */
import { runFVA, knockdownCurve } from './fba_engine.js';

const TOL = 1e-6;
const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c =>
  ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/* Data surfaces stay white in both themes (--canvas), so these are fixed. */
const C = {
  same: '#94A3B8', up: '#15803D', down: '#B45309',
  off: '#DC2626', on: '#2563EB', reversed: '#7C3AED',
  ink: '#334155', grid: '#E3E8EF', axis: '#64748B', wt: '#94A3B8', ko: '#1E40AF',
  warn: '#B45309',
};
const CLASS_LABEL = {
  up: 'carries more flux', down: 'carries less flux', off: 'switched off',
  on: 'switched on', reversed: 'runs backwards now', same: 'unchanged',
};

const CFG = {
  responsive: true, displaylogo: false,
  modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
  toImageButtonOptions: { format: 'svg', scale: 2 },   // vector out, editable in Illustrator
};
function layout(xt, yt, extra) {
  return Object.assign({
    margin: { l: 64, r: 18, t: 12, b: 52 },
    xaxis: { title: { text: xt || '', font: { size: 12 } }, automargin: true,
             gridcolor: C.grid, zerolinecolor: C.axis, tickfont: { size: 11 } },
    yaxis: { title: { text: yt || '', font: { size: 12 } }, automargin: true,
             gridcolor: C.grid, zerolinecolor: C.axis, tickfont: { size: 11 } },
    paper_bgcolor: '#FFFFFF', plot_bgcolor: '#FFFFFF',
    font: { family: "'Fira Sans', system-ui, sans-serif", size: 12, color: C.ink },
    hoverlabel: { font: { family: "'Fira Sans', system-ui, sans-serif", size: 12 } },
    showlegend: false,
  }, extra || {});
}

/* Log-modulus: signed, and defined at zero, which plain log is not. Fluxes span
   four orders of magnitude and go negative, so a linear axis shows one big point
   and a heap at the origin. */
const L = (x) => Math.sign(x) * Math.log10(1 + Math.abs(x));
function lmAxis(maxAbs) {
  const decades = [0, 1, 10, 100, 1000, 10000].filter(v => v <= Math.max(1, maxAbs) * 3);
  const vals = [...new Set([...decades.map(v => -v), ...decades])].sort((a, b) => a - b);
  return { tickvals: vals.map(L), ticktext: vals.map(String) };
}

function classify(wt, ko) {
  const a = Math.abs(wt) > TOL, b = Math.abs(ko) > TOL;
  if (!a && !b) return 'zero';
  if (a && !b) return 'off';
  if (!a && b) return 'on';
  if (Math.sign(wt) !== Math.sign(ko)) return 'reversed';
  if (Math.abs(ko) > Math.abs(wt) + TOL) return 'up';
  if (Math.abs(ko) < Math.abs(wt) - TOL) return 'down';
  return 'same';
}

/* The enzyme name, never the database id. Both pangenomes carry reaction names; the
   GPR carries locus tags, which are not gene symbols and mean nothing to a reader.
   The two reconstructions name their exchanges differently ("D-Glucose exchange" in
   one, "exchange reaction for D-Glucose" in the other), so both shapes are trimmed
   down to the compound. */
function label(r, max = 34) {
  let n = (r.name || '').trim()
    .replace(/^exchange\s+reaction\s+for\s+/i, '')
    .replace(/\s+exchange$/i, '');
  if (!n) n = r.id.replace(/^EX_/, '').replace(/_e$/, '');
  return n.length > max ? n.slice(0, max - 1) + '…' : n;
}

/* ── Subsystems ───────────────────────────────────────────────────────────────
   LactoPanGEM predates subsystem annotation and carries none. But a subsystem is
   a property of the BiGG reaction, not of the organism, so the assignments the
   EcopanGEM reconstructions carry transfer to any model using the same ids. That
   covers about half of a Lacto model; the rest is reported as unassigned rather
   than guessed at. */
let SUBSYS = null;
async function subsystems() {
  if (!SUBSYS) {
    try {
      const d = await (await fetch('fba/subsystems.json')).json();
      SUBSYS = { names: d.subsystems, map: d.map };
    } catch { SUBSYS = { names: [], map: {} }; }
  }
  return SUBSYS;
}
/* The E. coli reconstructions carry a literal "Unassigned" subsystem. It means the
   same thing as carrying none, so it collapses to null rather than becoming a bar of
   its own next to the reactions we could not place: one unassigned bucket, not two. */
const NO_SUBSYS = new Set(['', 'unassigned', 'biomass', 'other', 'none']);
/* The subsystem field is not clean in either pangenome, in two different ways.

   Some LactoPanGEM reconstructions store it as a stringified numpy array rather than a
   string: "[array(['Purine metabolism'], dtype='<U17')]", with an absent one written
   "[array([], dtype='<U1')]". Left alone, that repr goes straight onto the axis of the
   figure. And the value inside is often a pipe-joined list of KEGG maps, because KEGG
   assigns one enzyme to many pathways: "Glycolysis / Gluconeogenesis|Carbon fixation in
   photosynthetic organisms|Methane metabolism". Take the first, which is the primary. */
function cleanSubsystem(s) {
  if (s == null) return null;
  let t = String(s).trim();
  const arr = t.match(/array\(\s*\[([^\]]*)\]/);        // matched inside array([...]) so a
  if (arr) {                                            // loose quote match cannot return the dtype
    const q = arr[1].match(/'([^']*)'|"([^"]*)"/);
    t = q ? (q[1] || q[2] || '') : '';
  }
  t = t.split('|')[0].trim();
  return (t && !NO_SUBSYS.has(t.toLowerCase())) ? t : null;
}

/* ONE VOCABULARY PER CHART. BiGG subsystems and KEGG pathways name the same biology in
   different words, so filling the gaps in a KEGG-annotated model from the BiGG table
   would put "Pentose phosphate pathway" and "Pentose Phosphate Pathway" on the chart as
   two separate bars. So choose whichever source actually covers this model, and use it
   alone: the model's own annotation if it has one worth the name, otherwise the BiGG
   assignments carried by the reaction identifiers. */
function pickSource(model) {
  let own = 0;
  for (const r of model.reactions) if (cleanSubsystem(r.subsystem)) own++;
  return own >= 0.4 * model.reactions.length ? 'model' : 'bigg';
}
function subOf(r, S, src) {
  if (r.id.startsWith('EX_')) return 'Extracellular exchange';
  if (src === 'model') return cleanSubsystem(r.subsystem);
  const i = S.map[r.id];
  return i == null ? null : cleanSubsystem(S.names[i]);
}

/* ────────────────────────────────────────────────────────────────────────────
   ctx: {model, media, fW, fK, gW, gK, kos, method, note(msg)}
   ──────────────────────────────────────────────────────────────────────────── */
export async function renderKOPlots(ctx) {
  const P = window.Plotly;
  if (!P) return;
  const { model, fW, fK, gW, gK, kos } = ctx;
  const koSet = new Set(kos);
  const byId = new Map(model.reactions.map(r => [r.id, r]));

  // every reaction that moved, excluding the knocked-out ones themselves: those
  // did not re-route, they were removed. Counting them as re-routing inflates it.
  const moved = [];
  for (const r of model.reactions) {
    if (koSet.has(r.id)) continue;
    const a = fW[r.id] || 0, b = fK[r.id] || 0;
    const cls = classify(a, b);
    if (cls === 'zero' || cls === 'same') continue;
    moved.push({ id: r.id, r, wt: a, ko: b, d: b - a, cls });
  }
  moved.sort((x, y) => Math.abs(y.d) - Math.abs(x.d));

  const S = await subsystems();
  const dead = gK <= 1e-9;

  plotTitration(P, ctx);
  plotPathways(P, ctx, moved, S);
  plotScatter(P, ctx, moved, null);
  plotTornado(P, ctx, moved, null);
  plotExchanges(P, ctx);

  /* Now the honesty pass. Only the changes we are going to name get checked, and
     only if the knockout left the cell alive: if it is dead, every flux is zero
     and there is nothing to disambiguate. */
  const verdictBox = document.getElementById('ko-verdict');
  if (dead) {
    verdictBox.innerHTML = `<span class="kv-note">The knockout is lethal, so every
      flux is zero and there is no re-routing to verify.</span>`;
    return;
  }
  if (!moved.length) {
    verdictBox.innerHTML = `<span class="kv-note">No reaction changed flux. The cell
      absorbed the knockout without re-routing anything.</span>`;
    return;
  }

  const check = moved.slice(0, 40).map(m => m.id);
  verdictBox.innerHTML = `<span class="kv-run">Checking the ${check.length} largest
    changes against alternate optima (FVA on both states)…</span>`;
  let verdicts = null;
  try {
    const [vw, vk] = [
      await runFVA(model, ctx.media, check, { knockouts: [] }),
      await runFVA(model, ctx.media, check, { knockouts: kos }),
    ];
    if (vw.optimal && vk.optimal) {
      verdicts = new Map();
      for (const id of check) {
        const w = vw.ranges[id], k = vk.ranges[id];
        if (!w || !k) continue;
        const disjoint = (w.max < k.min - TOL) || (k.max < w.min - TOL);
        verdicts.set(id, { forced: disjoint, wt: w, ko: k });
      }
    }
  } catch (e) { /* leave verdicts null and say so */ }

  if (!verdicts) {
    verdictBox.innerHTML = `<span class="kv-note">Could not verify against alternate
      optima. Treat the differences below as one solution among several.</span>`;
    return;
  }

  const forced = [...verdicts.values()].filter(v => v.forced).length;
  const amb = verdicts.size - forced;
  verdictBox.innerHTML = `
    <span class="kv-chip ok">${forced} forced by the knockout</span>
    <span class="kv-chip amb">${amb} within alternate optima</span>
    <span class="kv-note">Of the ${verdicts.size} largest changes, ${forced} are real:
      the reaction's feasible flux range moved so far that no solution of the wild type
      and none of the knockout could share a value. The other ${amb} sit in ranges that
      still overlap, so the difference drawn is one valid choice among several and is
      shown hollow. The check is conservative: it explores every FBA optimum, a wider
      space than ${esc((ctx.method || 'pFBA').toUpperCase())} actually picks from.</span>`;

  plotScatter(P, ctx, moved, verdicts);
  plotTornado(P, ctx, moved, verdicts);
}

/* ── 1. Titration: how much of the reaction did the cell need? ─────────────── */
async function plotTitration(P, ctx) {
  const div = document.getElementById('ko-plot-titr');
  const cap = document.getElementById('ko-cap-titr');
  const carried = ctx.kos.filter(id => Math.abs(ctx.fW[id] || 0) > TOL);
  if (!carried.length) {
    div.innerHTML = `<div class="ko-empty">The knocked-out reaction carries no flux in
      the wild type, so there is nothing to titrate. Any growth change therefore comes
      from the network losing an option it was not using in this particular solution.</div>`;
    cap.textContent = '';
    return;
  }
  div.innerHTML = '<div class="ko-empty">Titrating…</div>';
  let res;
  try {
    res = await knockdownCurve(ctx.model, ctx.media, ctx.kos, { wtFluxes: ctx.fW, steps: 21 });
  } catch { div.innerHTML = '<div class="ko-empty">Titration failed.</div>'; return; }
  if (!res.optimal) { div.innerHTML = '<div class="ko-empty">Titration failed.</div>'; return; }

  const x = res.points.map(p => 100 * p.f), y = res.points.map(p => p.growth);
  const gW = ctx.gW;

  // where does growth first reach 99% of wild type? That is all the cell needed.
  let need = 100;
  for (const p of res.points) { if (p.growth >= 0.99 * gW) { need = 100 * p.f; break; } }

  /* Zoom to the curve. Anchoring the axis at zero would render a 3% growth change
     as a flat line, which is precisely the shape this plot exists to show. */
  const lo = Math.min(...y), hi = Math.max(...y);
  const span = Math.max(hi - lo, hi * 0.02, 1e-6);
  const yr = [lo - 0.18 * span, hi + 0.22 * span];

  P.newPlot(div, [
    { x, y, type: 'scatter', mode: 'lines', line: { color: C.ko, width: 2.5, shape: 'spline' },
      hovertemplate: 'capacity kept %{x:.0f}%<br>growth %{y:.4f} h⁻¹<extra></extra>' },
    { x: [0], y: [ctx.gK], type: 'scatter', mode: 'markers',
      marker: { color: C.off, size: 12, line: { color: '#fff', width: 1.5 } },
      hovertemplate: 'full knockout<br>growth %{y:.4f} h⁻¹<extra></extra>' },
    { x: [100], y: [gW], type: 'scatter', mode: 'markers',
      marker: { color: C.same, size: 12, line: { color: '#fff', width: 1.5 } },
      hovertemplate: 'wild type<br>growth %{y:.4f} h⁻¹<extra></extra>' },
  ], layout('Flux capacity kept, % of wild type', 'Growth (h⁻¹)', {
    margin: { l: 70, r: 18, t: 20, b: 52 },
    yaxis: { title: { text: 'Growth (h⁻¹)', font: { size: 12 } }, range: yr, automargin: true,
             gridcolor: C.grid, zerolinecolor: C.grid, tickfont: { size: 11 },
             tickformat: '.3f' },
    xaxis: { title: { text: 'Flux capacity kept, % of wild type', font: { size: 12 } },
             range: [-4, 104], gridcolor: C.grid, zerolinecolor: C.grid, tickfont: { size: 11 } },
    shapes: [
      { type: 'line', x0: -4, x1: 104, y0: gW, y1: gW,                    // the wild-type ceiling
        line: { color: C.same, width: 1.2, dash: 'dash' } },
      ...(need > 0.01 && need < 99 ? [{ type: 'line', x0: need, x1: need, y0: yr[0], y1: gW,
        line: { color: C.warn, width: 1.3, dash: 'dot' } }] : []),
    ],
    annotations: [
      { x: 2, y: gW, text: 'wild type', showarrow: false, xanchor: 'left', yshift: 9,
        font: { size: 10, color: C.axis } },
      // sits under the curve, so it cannot collide with the wild-type label above it
      ...(need > 0.01 && need < 99 ? [{ x: need, y: yr[0], text: `needs ${need.toFixed(0)}%`,
        showarrow: false, yshift: 12, xanchor: need > 55 ? 'right' : 'left', xshift: need > 55 ? -4 : 4,
        font: { size: 11, color: C.warn } }] : []),
    ],
  }), CFG);

  /* What the curve is worth knowing for: how much has to be taken away before anything
     happens. A reaction the cell over-provisions is a poor target, because inhibiting
     half of it does nothing at all. So always state the spare capacity, rather than
     sorting the answer into buckets whose edges would decide the wording. */
  const kept = ctx.gK / gW;
  const spare = 100 - need;
  if (kept >= 0.99) {
    cap.innerHTML = `Flat. Growth is unchanged even with the reaction removed entirely, so
      the cell has a route around it and no amount of inhibition would slow it down. The
      flux it carried in the wild type was a preference, not a requirement.`;
  } else if (need >= 99) {
    cap.innerHTML = `Growth falls the moment any capacity is taken away: no spare at all, and
      every percent of inhibition costs a percent of growth. Removed entirely, the cell still
      manages <b>${(100 * kept).toFixed(0)}%</b> of wild-type growth on another route.`;
  } else {
    cap.innerHTML = `The cell needs only <b>${need.toFixed(0)}%</b> of the flux it was
      carrying, so <b>${spare.toFixed(0)}%</b> of this reaction's activity is spare. An
      inhibitor would have to remove more than that before growth moved at all${
      spare >= 30 ? ', which is what separates a drug target from a dead end' : ''}.
      Removed entirely, the cell still manages <b>${(100 * kept).toFixed(0)}%</b> of
      wild-type growth.`;
  }
}

/* ── 2. Where in metabolism did it land? ──────────────────────────────────── */
function plotPathways(P, ctx, moved, S) {
  const div = document.getElementById('ko-plot-subsys');
  const cap = document.getElementById('ko-cap-subsys');
  const src = pickSource(ctx.model);

  /* Group case-insensitively: the reconstructions are inconsistent about it and would
     otherwise show "Transport, extracellular" and "transport, extracellular" as two
     different pathways. Display whichever spelling is the more common. */
  const agg = new Map(); let unassigned = 0;
  for (const m of moved) {
    const s = subOf(m.r, S, src);
    if (!s) { unassigned++; continue; }
    const k = s.toLowerCase();
    const a = agg.get(k) || { flux: 0, n: 0, forms: new Map() };
    a.flux += Math.abs(m.d); a.n++;
    a.forms.set(s, (a.forms.get(s) || 0) + 1);
    agg.set(k, a);
  }
  const rows = [...agg.values()].map(a => ({
    s: [...a.forms.entries()].sort((x, y) => y[1] - x[1])[0][0], flux: a.flux, n: a.n,
  })).sort((a, b) => b.flux - a.flux).slice(0, 11).reverse();
  if (!rows.length) {
    div.innerHTML = '<div class="ko-empty">No pathway annotation for the reactions that changed.</div>';
    cap.textContent = ''; return;
  }
  const short = (s) => s.length > 38 ? s.slice(0, 37) + '…' : s;

  // lollipop: the stem carries the magnitude, the head carries the count
  const traces = [];
  rows.forEach((r, i) => traces.push({
    x: [0, r.flux], y: [i, i], type: 'scatter', mode: 'lines',
    line: { color: C.grid, width: 2 }, hoverinfo: 'skip',
  }));
  traces.push({
    x: rows.map(r => r.flux), y: rows.map((_, i) => i),
    text: rows.map(r => String(r.n)),
    type: 'scatter', mode: 'markers+text', textposition: 'middle right',
    textfont: { size: 10, color: C.axis },
    marker: {
      color: rows.map(r => r.flux), colorscale: [[0, '#BFDBFE'], [1, C.ko]],
      size: rows.map(r => Math.min(26, 9 + 2.2 * Math.sqrt(r.n))),
      line: { color: '#fff', width: 1.5 },
    },
    customdata: rows.map(r => [r.n]),
    hovertemplate: '%{y}<br>%{customdata[0]} reactions changed<br>total |Δ| %{x:.2f} mmol gDW⁻¹ h⁻¹<extra></extra>',
  });
  P.newPlot(div, traces, layout('Total flux moved (mmol gDW⁻¹ h⁻¹)', '', {
    margin: { l: 8, r: 34, t: 12, b: 52 },
    yaxis: { tickmode: 'array', tickvals: rows.map((_, i) => i),
             ticktext: rows.map(r => short(r.s)), automargin: true, gridcolor: 'rgba(0,0,0,0)',
             tickfont: { size: 11 } },
  }), CFG);

  cap.innerHTML = `Each pathway's bar is the total flux that moved inside it; the number
    on the head is how many of its reactions changed.
    ${src === 'bigg'
      ? `This reconstruction carries no pathway annotation, so the pathways come from the
         BiGG identifiers its reactions share with the E.&nbsp;coli models.`
      : ''}
    ${unassigned
      ? `<b>${unassigned}</b> changed reactions could not be placed in a pathway and are not shown.`
      : ''}`;
}

/* ── 3. The whole redistribution, in one view ─────────────────────────────── */
function plotScatter(P, ctx, moved, verdicts) {
  const div = document.getElementById('ko-plot-scatter');
  const cap = document.getElementById('ko-cap-scatter');
  const { model, fW, fK } = ctx;
  const koSet = new Set(ctx.kos);

  const pts = [];
  for (const r of model.reactions) {
    if (koSet.has(r.id)) continue;
    const a = fW[r.id] || 0, b = fK[r.id] || 0;
    if (Math.abs(a) <= TOL && Math.abs(b) <= TOL) continue;
    pts.push({ r, wt: a, ko: b, cls: classify(a, b) });
  }
  const maxAbs = Math.max(1, ...pts.map(p => Math.max(Math.abs(p.wt), Math.abs(p.ko))));
  const ax = lmAxis(maxAbs);
  const lim = L(maxAbs) * 1.08;

  const order = ['same', 'up', 'down', 'off', 'on', 'reversed'];
  const traces = [{
    x: [-lim, lim], y: [-lim, lim], type: 'scatter', mode: 'lines',
    line: { color: C.grid, width: 1.5, dash: 'dash' }, hoverinfo: 'skip',
  }];
  for (const cls of order) {
    const g = pts.filter(p => p.cls === cls);
    if (!g.length) continue;
    const forced = (p) => verdicts && verdicts.has(p.r.id) ? verdicts.get(p.r.id).forced : null;
    traces.push({
      x: g.map(p => L(p.wt)), y: g.map(p => L(p.ko)),
      type: 'scatter', mode: 'markers', name: cls,
      marker: {
        color: g.map(p => forced(p) === false ? '#FFFFFF' : C[cls]),
        size: g.map(p => cls === 'same' ? 5 : (forced(p) === true ? 10 : 8)),
        opacity: cls === 'same' ? 0.5 : 0.9,
        line: { color: g.map(p => C[cls]), width: g.map(p => forced(p) === false ? 1.8 : 0.6) },
      },
      customdata: g.map(p => {
        const v = verdicts && verdicts.get(p.r.id);
        return [label(p.r, 46), p.wt.toFixed(3), p.ko.toFixed(3), CLASS_LABEL[p.cls],
                v ? (v.forced ? 'forced by the knockout' : 'within alternate optima') : ''];
      }),
      hovertemplate: '<b>%{customdata[0]}</b><br>wild type %{customdata[1]}  →  knockout %{customdata[2]}'
                   + '<br>%{customdata[3]}<br><i>%{customdata[4]}</i><extra></extra>',
    });
  }
  P.newPlot(div, traces, layout('Wild-type flux (mmol gDW⁻¹ h⁻¹)', 'Knockout flux (mmol gDW⁻¹ h⁻¹)', {
    margin: { l: 70, r: 18, t: 12, b: 56 },
    xaxis: Object.assign({ title: { text: 'Wild-type flux (mmol gDW⁻¹ h⁻¹)', font: { size: 12 } },
      range: [-lim, lim], gridcolor: C.grid, zerolinecolor: C.axis, tickfont: { size: 11 } }, ax),
    yaxis: Object.assign({ title: { text: 'Knockout flux (mmol gDW⁻¹ h⁻¹)', font: { size: 12 } },
      range: [-lim, lim], gridcolor: C.grid, zerolinecolor: C.axis, tickfont: { size: 11 } }, ax),
  }), CFG);

  const n = (c) => pts.filter(p => p.cls === c).length;
  cap.innerHTML = `Every reaction carrying flux, wild type against knockout, on a signed
    log scale so that zeros and reversals are both visible. Points on the dashed diagonal
    did not change. <span class="lg" style="--c:${C.off}">${n('off')} switched off</span>
    <span class="lg" style="--c:${C.on}">${n('on')} switched on</span>
    <span class="lg" style="--c:${C.reversed}">${n('reversed')} reversed direction</span>
    <span class="lg" style="--c:${C.up}">${n('up')} up</span>
    <span class="lg" style="--c:${C.down}">${n('down')} down</span>
    ${ctx.gK <= 1e-9
      ? 'The knockout is lethal, so everything the cell was running now sits on the zero line.'
      : (verdicts ? '<b>Hollow</b> points are changes that alternate optima could also explain.' : '')}`;
}

/* ── 4. Which reactions, by name ──────────────────────────────────────────── */
function plotTornado(P, ctx, moved, verdicts) {
  const div = document.getElementById('ko-plot-tornado');
  const cap = document.getElementById('ko-cap-tornado');
  const top = moved.slice(0, 18).reverse();
  if (!top.length) { div.innerHTML = '<div class="ko-empty">Nothing re-routed.</div>'; cap.textContent = ''; return; }

  const forced = (m) => verdicts && verdicts.has(m.id) ? verdicts.get(m.id).forced : null;
  P.newPlot(div, [{
    x: top.map(m => m.d), y: top.map((_, i) => i),
    type: 'bar', orientation: 'h',
    marker: {
      color: top.map(m => forced(m) === false ? 'rgba(255,255,255,0)' : C[m.cls]),
      line: { color: top.map(m => C[m.cls]), width: top.map(m => forced(m) === false ? 2 : 0) },
    },
    customdata: top.map(m => {
      const v = verdicts && verdicts.get(m.id);
      return [label(m.r, 46), m.wt.toFixed(3), m.ko.toFixed(3), CLASS_LABEL[m.cls],
              v ? (v.forced ? 'forced by the knockout' : 'within alternate optima') : ''];
    }),
    hovertemplate: '<b>%{customdata[0]}</b><br>wild type %{customdata[1]}  →  knockout %{customdata[2]}'
                 + '<br>%{customdata[3]}<br><i>%{customdata[4]}</i><extra></extra>',
  }], layout('Change in flux, knockout − wild type', '', {
    margin: { l: 8, r: 24, t: 12, b: 52 },
    yaxis: { tickmode: 'array', tickvals: top.map((_, i) => i),
             ticktext: top.map(m => label(m.r, 30)), automargin: true,
             gridcolor: 'rgba(0,0,0,0)', tickfont: { size: 11 } },
    bargap: 0.35,
  }), CFG);

  cap.innerHTML = `The ${top.length} reactions whose flux moved most, by enzyme name.
    Bars right of zero carry more flux after the knockout, bars left carry less.
    ${verdicts ? '<b>Outlined</b> bars are changes that alternate optima could also explain.' : ''}`;
}

/* ── 5. What the cell eats and excretes ───────────────────────────────────── */
function plotExchanges(P, ctx) {
  const div = document.getElementById('ko-plot-exch');
  const cap = document.getElementById('ko-cap-exch');
  const { model, fW, fK } = ctx;

  const all = [];
  for (const r of model.reactions) {
    if (!r.id.startsWith('EX_')) continue;
    const a = fW[r.id] || 0, b = fK[r.id] || 0;
    if (Math.abs(a) <= TOL && Math.abs(b) <= TOL) continue;
    if (Math.abs(b - a) <= TOL) continue;                  // unchanged boundary: not news
    all.push({ r, wt: a, ko: b, d: b - a });
  }
  if (!all.length) {
    div.innerHTML = `<div class="ko-empty">The cell eats and excretes exactly the same
      things, at the same rates. The knockout was absorbed entirely inside the network.</div>`;
    cap.textContent = ''; return;
  }
  all.sort((a, b) => Math.abs(b.d) - Math.abs(a.d));

  /* Every uptake that is coupled to biomass moves a little whenever growth moves, so a
     knockout that costs 3% growth nudges every salt in the medium by about 3%. Those
     rows are real but they are not the story, and on a shared axis they crowd out the
     compounds that actually changed. Keep only shifts worth a reader's attention. */
  const maxD = Math.abs(all[0].d);
  const rows = all.filter(m => Math.abs(m.d) >= Math.max(1e-4, 0.02 * maxD));
  const hidden = all.length - rows.length;
  const top = rows.slice(0, 14).reverse();

  const traces = [];
  top.forEach((m, i) => traces.push({
    x: [m.wt, m.ko], y: [i, i], type: 'scatter', mode: 'lines',
    line: { color: m.d > 0 ? C.up : C.down, width: 2.5 }, hoverinfo: 'skip',
  }));
  traces.push({
    x: top.map(m => m.wt), y: top.map((_, i) => i), type: 'scatter', mode: 'markers',
    marker: { color: '#FFFFFF', size: 10, line: { color: C.wt, width: 2 } },
    customdata: top.map(m => [label(m.r, 40)]),
    hovertemplate: '<b>%{customdata[0]}</b><br>wild type %{x:.3f}<extra></extra>',
  });
  traces.push({
    x: top.map(m => m.ko), y: top.map((_, i) => i), type: 'scatter', mode: 'markers',
    marker: { color: C.ko, size: 10, line: { color: '#fff', width: 1.5 } },
    customdata: top.map(m => [label(m.r, 40)]),
    hovertemplate: '<b>%{customdata[0]}</b><br>knockout %{x:.3f}<extra></extra>',
  });

  /* The dots can sit close together when a shift is small next to the axis span, so
     the size of the move is written out rather than left to be measured in pixels. */
  const lo = Math.min(0, ...top.map(m => Math.min(m.wt, m.ko)));
  const hi = Math.max(0, ...top.map(m => Math.max(m.wt, m.ko)));
  const pad = 0.16 * Math.max(hi - lo, 1e-6);

  P.newPlot(div, traces, layout('Exchange flux  (negative = taken up, positive = secreted)', '', {
    margin: { l: 8, r: 24, t: 12, b: 56 },
    xaxis: { title: { text: 'Exchange flux  (negative = taken up, positive = secreted)', font: { size: 12 } },
             range: [lo - pad, hi + pad * 1.5], automargin: true,
             gridcolor: C.grid, zerolinecolor: C.axis, tickfont: { size: 11 } },
    yaxis: { tickmode: 'array', tickvals: top.map((_, i) => i),
             ticktext: top.map(m => label(m.r, 26)), automargin: true,
             gridcolor: 'rgba(0,0,0,0)', tickfont: { size: 11 },
             range: [-0.7, top.length - 0.3] },
    shapes: [{ type: 'line', x0: 0, x1: 0, yref: 'paper', y0: 0, y1: 1,
               line: { color: C.axis, width: 1 } }],
    annotations: top.map((m, i) => ({
      x: Math.max(m.wt, m.ko), y: i, xanchor: 'left', xshift: 9,
      text: (m.d > 0 ? '+' : '') + m.d.toFixed(2),
      showarrow: false, font: { size: 10, color: m.d > 0 ? C.up : C.down },
    })),
  }), CFG);

  const crossed = top.filter(m => Math.sign(m.wt) !== Math.sign(m.ko) && Math.abs(m.wt) > TOL && Math.abs(m.ko) > TOL);
  cap.innerHTML = `The organism's boundary, before and after. Hollow dot is the wild type,
    filled dot the knockout, and the bar between them is the shift. Left of the line the
    cell is feeding, right of it excreting.
    ${crossed.length
      ? `<b>${crossed.map(m => label(m.r, 18)).join(', ')}</b> cross the line: the cell used
         to take ${crossed.length > 1 ? 'them' : 'it'} up and now secretes
         ${crossed.length > 1 ? 'them' : 'it'}, or the reverse.`
      : `Nothing crosses the line, so the cell eats and excretes the same set of compounds
         and only the rates moved.`}
    ${rows.length > top.length ? `Showing the ${top.length} largest of ${rows.length}. ` : ''}
    ${hidden ? `${hidden} further exchange${hidden > 1 ? 's' : ''} shifted by under 2% of the
      largest, which is the whole medium drifting with biomass rather than a real change,
      and ${hidden > 1 ? 'are' : 'is'} not shown.` : ''}`;
}
