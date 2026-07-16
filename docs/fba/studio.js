// EcopanGEM Flux Analysis Studio — FVA, Dynamic FBA & Multi-model tabs.
// (The Explore/Compare tab is handled by fba_ui.js.) All client-side.
import { runFBA, runPFBA, runFVA, runDFBA, productionEnvelope, phasePlane, essentialityScan, exchangeReport, listExchanges, bindMedium,
         singleGeneDeletion, doubleDeletion, fseof, sampleFluxes, findBlockedReactions, modelQC, listGenes } from './fba_engine.js';

// Glucose is spelled differently across BiGG generations. Pick whichever this model has.
import { mountMediaPicker, rebind as rebindMedia, adoptDefaultMedium, presets as mediaPresets } from './media_ui.js';

const GLC = ['EX_glc__D_e', 'EX_glc_D_e', 'EX_glc_e'];
const pickGlc = (ids) => GLC.find(g => ids.includes(g)) || null;

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (x, n = 3) => (x == null || isNaN(x)) ? '—' : Number(x).toFixed(n);
let PRESETS = null;

const PHYLO_COLORS = { A: '#4e79a7', B1: '#59a14f', B2: '#e15759', C: '#f28e2b', D: '#76b7b2',
  E: '#edc948', F: '#b07aa1', G: '#ff9da7', clade: '#9c755f', Unknown: '#bab0ac', '': '#bab0ac' };
const phyloColor = (p) => PHYLO_COLORS[p] || '#8888aa';

// ── shared data helpers ───────────────────────────────────────────────────────
async function presets() { if (!PRESETS) PRESETS = await (await fetch('fba/media_presets.json')).json(); return PRESETS; }
function fillMedia(sel, descEl) {
  sel.innerHTML = '';
  for (const [k, v] of Object.entries(PRESETS)) { const o = document.createElement('option'); o.value = k; o.textContent = v.label; sel.appendChild(o); }
  const upd = () => { if (descEl) { const p = PRESETS[sel.value]; descEl.textContent = `${p.desc} · ${Object.keys(p.bounds).length} exchanges · ${p.source}`; } };
  sel.addEventListener('change', upd); upd();
}
/* Bind a preset medium to THIS model. Never pass raw preset bounds to the solver:
   LactoPanGEM predates the BiGG double-underscore convention (EX_glc_D_e vs
   EX_glc__D_e), buildLP closes every exchange the medium does not name, and the
   result is a silently carbon-free medium and zero growth. See bindMedium. */
function mediaFor(model, key) {
  const p = PRESETS[key];
  if (!model || !p) return (p ? p.bounds : {});
  const comps = Object.entries(p.bounds).map(([exchange, lower_bound]) => ({ exchange, lower_bound }));
  return bindMedium(model, comps).bounds;
}

/* The medium a tab currently has selected, already bound to its model. Falls back
   to the preset dropdown if the picker has not mounted yet. */
function mediaOf(state, selId) {
  if (state.mediaBounds && Object.keys(state.mediaBounds).length) return state.mediaBounds;
  return mediaFor(state.model, $(selId).value);
}

/* Mount the shared picker on a tab and keep it bound as the model changes. */
function attachMedia(selId, state) {
  mountMediaPicker($(selId), state);
  const orig = state.onModelLoaded;
  state.onModelLoaded = () => { adoptDefaultMedium(state, $(selId)); if (orig) orig(); };
}

const metaByFile = new Map();
function meta(gemFile) { if (!metaByFile.size && window.gemMetadata) window.gemMetadata.forEach(m => metaByFile.set(m.gem_file, m)); return metaByFile.get(gemFile) || {}; }

const modelCache = new Map();
async function loadModel(gemFile) {
  if (modelCache.has(gemFile)) return modelCache.get(gemFile);
  const bn = window.gemBatchMap && window.gemBatchMap[gemFile];
  if (!bn) throw new Error(`Unknown model "${gemFile}".`);
  const resp = await fetch(window.BATCH_URL_BASE + String(bn).padStart(2, '0') + '.zip');
  if (!resp.ok) throw new Error('Batch download failed.');
  const zip = await window.JSZip.loadAsync(await resp.blob());
  const entry = zip.file(gemFile); if (!entry) throw new Error('Model not in batch.');
  const model = JSON.parse(await entry.async('text')); modelCache.set(gemFile, model); return model;
}

// reusable model combobox. onPick(gemFile). If multi, input keeps focus.
/* ── Deep links ───────────────────────────────────────────────────────────────
   The model browsers hand a strain straight to an analysis here, so a reader who has
   just found an interesting genome does not have to come and search for it again:

     ?model=<gem_file>[&tab=<analysis>][&medium=<preset>][&ko=<rxn,rxn>]
     ?models=<a,b,c>                       -> multi-model comparison

   Every model picker goes through makeCombo, so registering the loader here catches
   all eleven of them without touching each call site. */
const TAB_BY_PREFIX = { fva: 'fva', dfba: 'dfba', mm: 'multi', env: 'envelope', pp: 'phaseplane',
  ess: 'essential', genes: 'genes', syn: 'synlethal', fseof: 'design', samp: 'sampling', qc: 'qc' };
const MEDIA_SELECT = { fva: 'fva-media', dfba: 'dfba-media', multi: 'mm-media', cohort: 'cohort-media',
  envelope: 'env-media', phaseplane: 'pp-media', essential: 'ess-media', genes: 'genes-media',
  synlethal: 'syn-media', design: 'fseof-media', sampling: 'samp-media', qc: 'qc-media' };
const MODEL_LOADERS = {};

function makeCombo(input, menu, onPick) {
  const tab = TAB_BY_PREFIX[(input.id || '').replace(/-model-input$/, '')];
  if (tab) MODEL_LOADERS[tab] = onPick;

  let items = [], active = -1;
  const render = () => {
    const q = input.value.trim().toLowerCase();
    const md = window.gemMetadata || [];
    items = (q ? md.filter(m => (m.gem_file || '').toLowerCase().includes(q) || (m.genome_name || '').toLowerCase().includes(q) || (m.strain || '').toLowerCase().includes(q)) : md).slice(0, 40);
    menu.innerHTML = items.length ? items.map((m, i) => `<div class="fba-combo-item${i === active ? ' active' : ''}" data-i="${i}"><div class="nm">${esc(m.genome_name || m.strain || m.gem_file)}</div><div class="id">${esc(m.gem_file)}</div><div class="meta">${m.phylogroup ? 'Phylogroup ' + esc(m.phylogroup) : ''}${m.MLST ? ' · ST' + esc(m.MLST) : ''}${m.isolation_source ? ' · ' + esc(m.isolation_source) : ''}</div></div>`).join('') : `<div class="fba-combo-empty">No models match.</div>`;
    menu.classList.add('show');
  };
  input.addEventListener('focus', render); input.addEventListener('input', () => { active = -1; render(); });
  input.addEventListener('keydown', (e) => {
    if (!menu.classList.contains('show')) return;
    if (e.key === 'ArrowDown') { active = Math.min(active + 1, items.length - 1); render(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { active = Math.max(active - 1, 0); render(); e.preventDefault(); }
    else if (e.key === 'Enter') { const it = items[active] || items[0]; if (it) { onPick(it.gem_file); e.preventDefault(); } }
    else if (e.key === 'Escape') menu.classList.remove('show');
  });
  menu.addEventListener('mousedown', (e) => { const it = e.target.closest('.fba-combo-item'); if (!it) return; e.preventDefault(); onPick(items[+it.dataset.i].gem_file); });
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));
}

function modelCardHTML(gemFile, model) {
  const m = meta(gemFile);
  return `<div class="mc-name">${esc(m.genome_name || m.strain || gemFile)}</div><div class="mc-id">${esc(gemFile)}</div>
    <div class="mc-tags"><span class="fba-tag">${model.reactions.length} rxns</span><span class="fba-tag">${model.metabolites.length} mets</span>
    ${m.phylogroup ? `<span class="fba-tag">Phylogroup ${esc(m.phylogroup)}</span>` : ''}${m.MLST ? `<span class="fba-tag">ST ${esc(m.MLST)}</span>` : ''}${m.isolation_source ? `<span class="fba-tag">${esc(m.isolation_source)}</span>` : ''}</div>`;
}
function setStatus(id, msg, kind) { const e = $(id); e.textContent = msg; e.className = 'fba-status ' + (kind || ''); }
/* shared KPI + Plotly layout for the newer panels */
const PLOT_CFG = { responsive: true, displaylogo: false };
const PALETTE = ['#2c6fbb', '#c0392b', '#1a7f4b', '#e08a1e', '#7d3c98', '#16a085',
                 '#5b8ff9', '#e15759', '#59a14f', '#f28e2b', '#76b7b2', '#b07aa1'];

function kpi(v, l, color) {
  return `<div class="fba-kpi"><div class="v"${color ? ` style="color:${color}"` : ''}>${v}</div><div class="l">${l}</div></div>`;
}
function plotly(title, xt, yt) {
  return {
    title: title || '',
    margin: { l: 62, r: 20, t: title ? 40 : 18, b: 62 },
    xaxis: { title: xt || '', automargin: true },
    yaxis: { title: yt || '', automargin: true },
    paper_bgcolor: '#fff', plot_bgcolor: '#fff',
    font: { family: "'Fira Sans', system-ui, sans-serif", size: 11, color: '#334155' },
  };
}

function prog(wrapId, barId, frac) { const w = $(wrapId); w.style.display = frac == null ? 'none' : 'block'; if (frac != null) $(barId).style.width = Math.round(frac * 100) + '%'; }
function saveCSV(csv, base) { const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' })); a.download = base + '.csv'; a.click(); }
const bio = (id) => id.replace(/^EX_/, '').replace(/_e$/, '');

// ── Sidebar navigation ─────────────────────────────────────────────────────────
const TAB_META = {
  genes: ['Gene Essentiality', 'Delete every gene through its GPR rule &mdash; isozymes survive, complexes do not.'],
  synlethal: ['Synthetic Lethality', 'Double deletions: pairs where neither knockout hurts alone, but together they kill.'],
  design: ['Strain Design (FSEOF)', 'Force a product up and see which reactions must carry more flux.'],
  sampling: ['Flux Sampling', 'ACHR hit-and-run through the solution space &mdash; the distribution, not one optimum.'],
  qc: ['Model QC', 'Blocked reactions, mass and charge imbalance, dead ends, orphan genes.'],
  home: ['Flux Studio', 'A constraint-based modelling workbench for 4,659 strain-specific metabolic models — everything runs in your browser.'],
  explore: ['Explore & Compare', 'FBA / pFBA with editable media, knockouts and a live Escher flux map.'],
  dfba: ['Dynamic FBA', 'Batch-culture time course from Michaelis–Menten uptake kinetics.'],
  fva: ['Flux Variability Analysis', 'The min/max flux each reaction can carry at near-optimal growth.'],
  envelope: ['Production Envelope', 'The biomass-vs-product trade-off frontier for strain design.'],
  phaseplane: ['Phenotype Phase Plane', 'The growth surface over two uptake capacities.'],
  essential: ['Essentiality Screen', 'Single-reaction knockout scan across the network.'],
  multi: ['Multi-model analytics', 'FBA across many strains — scatter, PCA and heatmaps.'],
  cohort: ['Group comparison', 'Test which metabolic traits differ between metadata-defined groups.'],
};
function switchTab(tab) {
  if (!TAB_META[tab]) tab = 'home';
  document.querySelectorAll('#studio-nav .nav-item').forEach(b => b.classList.toggle('active', b.dataset.tab === tab));
  document.querySelectorAll('.studio-panel').forEach(p => p.classList.toggle('active', p.id === 'panel-' + tab));
  const [t, s] = TAB_META[tab]; $('page-title').textContent = t; $('page-sub').innerHTML = s;
  document.getElementById('sidebar').classList.remove('open');
  window.scrollTo({ top: 0, behavior: 'smooth' });
  setTimeout(() => document.querySelectorAll('#panel-' + tab + ' .js-plotly-plot').forEach(p => window.Plotly && window.Plotly.Plots.resize(p)), 40);
}
function initNav() {
  document.querySelectorAll('#studio-nav .nav-item, #tool-gallery .tool-card').forEach(b => b.addEventListener('click', () => switchTab(b.dataset.tab)));
  const mt = $('menu-toggle'); if (mt) mt.addEventListener('click', () => document.getElementById('sidebar').classList.toggle('open'));
}

// ── Plot downloads ──────────────────────────────────────────────────────────────
function wireDownloads() {
  document.addEventListener('click', e => {
    const btn = e.target.closest('.viz-dl'); if (!btn) return;
    const el = document.getElementById(btn.dataset.plot); if (!el) return;
    const type = btn.dataset.type || 'plotly';
    if (type === 'plotly' && window.Plotly) window.Plotly.downloadImage(el, { format: 'png', scale: 2, filename: btn.dataset.plot, width: el.clientWidth || 900, height: Math.max(360, el.clientHeight || 460) });
    else if (type === 'chartjs') { const a = document.createElement('a'); a.href = el.toDataURL('image/png'); a.download = btn.dataset.plot + '.png'; a.click(); }
    else if (type === 'escher') { const svg = el.querySelector('svg'); if (svg) { const a = document.createElement('a'); a.href = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(new XMLSerializer().serializeToString(svg)); a.download = btn.dataset.plot + '.svg'; a.click(); } }
  });
}
// add a ⬇ button to every static Plotly card that doesn't already have one
function decorateStaticPlots() {
  document.querySelectorAll('.fba-chart-card').forEach(card => {
    const h = card.querySelector('h6'); if (!h || h.querySelector('.viz-dl')) return;
    const pd = card.querySelector('.plot, .plot-tall'); if (!pd || !pd.id) return;
    const b = document.createElement('button'); b.className = 'viz-dl'; b.dataset.plot = pd.id; b.dataset.type = 'plotly';
    b.textContent = '⬇'; b.style.cssText = 'float:right;margin-top:-1px'; h.appendChild(b);
  });
}

// ── FVA tab ───────────────────────────────────────────────────────────────────
const fvaState = { model: null, file: null };
function initFVA() {
  fillMedia($('fva-media'), $('fva-media-desc'));
  attachMedia('fva-media', fvaState);
  makeCombo($('fva-model-input'), $('fva-model-menu'), async (gemFile) => {
    $('fva-model-input').value = '';
    const mc = $('fva-modelcard'); mc.style.display = 'block'; mc.innerHTML = '<span class="fba-hint-inline">Loading…</span>';
    try { const model = await loadModel(gemFile); fvaState.model = model; fvaState.file = gemFile; mc.innerHTML = modelCardHTML(gemFile, model); }
    catch (e) { mc.innerHTML = `<span style="color:var(--bad)">${esc(e.message)}</span>`; }
  });
  $('fva-frac').addEventListener('input', () => $('fva-frac-val').textContent = (+$('fva-frac').value).toFixed(2));
  $('fva-run').addEventListener('click', runFVAtab);
}
async function runFVAtab() {
  if (!fvaState.model) return setStatus('fva-status', 'Choose a model first.', 'err');
  const media = mediaOf(fvaState, 'fva-media');
  const fraction = +$('fva-frac').value;
  const exIds = listExchanges(fvaState.model).map(e => e.id);
  $('fva-run').disabled = true; $('fva-results').style.display = 'none';
  setStatus('fva-status', `Running FVA on ${exIds.length} exchanges…`, 'busy'); prog('fva-prog', 'fva-prog-bar', 0);
  try {
    const t0 = performance.now();
    const res = await runFVA(fvaState.model, media, exIds, { fraction, onProgress: (d, t) => prog('fva-prog', 'fva-prog-bar', d / t) });
    prog('fva-prog', 'fva-prog-bar', null);
    if (!res.optimal) { setStatus('fva-status', 'Model does not grow on this medium.', 'err'); $('fva-run').disabled = false; return; }
    renderFVA(res, fvaState.model);
    setStatus('fva-status', `Done — ${exIds.length} exchanges in ${((performance.now() - t0) / 1000).toFixed(1)} s. Optimum growth ${fmt(res.z)} h⁻¹.`, 'ok');
  } catch (e) { setStatus('fva-status', 'Error: ' + e.message, 'err'); console.error(e); prog('fva-prog', 'fva-prog-bar', null); }
  finally { $('fva-run').disabled = false; }
}
function renderFVA(res, model) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  // keep exchanges with non-trivial range or nonzero, sort by span
  let rows = Object.entries(res.ranges).map(([id, r]) => ({ id, name: nameById[id], min: r.min, max: r.max, span: r.max - r.min }))
    .filter(r => Math.abs(r.max) > 1e-6 || Math.abs(r.min) > 1e-6);
  rows.sort((a, b) => b.span - a.span);
  const top = rows.slice(0, 30).reverse();
  $('fva-results').style.display = 'block';
  window.Plotly.newPlot('fva-plot', [{
    type: 'scatter', mode: 'markers', x: top.map(r => r.min), y: top.map(bioName), name: 'min',
    marker: { color: '#2c6fbb', size: 8, symbol: 'line-ns-open' }, hovertemplate: '%{y}<br>min %{x:.3f}<extra></extra>'
  }, {
    type: 'scatter', mode: 'markers', x: top.map(r => r.max), y: top.map(bioName), name: 'max',
    marker: { color: '#c0392b', size: 8 }, hovertemplate: '%{y}<br>max %{x:.3f}<extra></extra>'
  }, {
    type: 'scatter', mode: 'lines', x: top.flatMap(r => [r.min, r.max, null]), y: top.flatMap(r => [bioName(r), bioName(r), null]),
    line: { color: '#9db8d6', width: 6 }, hoverinfo: 'skip', showlegend: false
  }], {
    margin: { l: 110, r: 20, t: 10, b: 40 }, height: Math.max(420, top.length * 16),
    xaxis: { title: 'Flux range (mmol·gDW⁻¹·h⁻¹) — biomass ≥ ' + res.fraction + '× optimum', zeroline: true, zerolinecolor: '#ccc' },
    yaxis: { automargin: true }, legend: { orientation: 'h', y: 1.05 }, font: { size: 11 }
  }, { responsive: true, displaylogo: false });
  function bioName(r) { return bio(r.id); }
  // table
  const tr = rows.map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td><td class="num">${fmt(r.min)}</td><td class="num">${fmt(r.max)}</td><td class="num">${fmt(r.span)}</td></tr>`).join('');
  $('fva-table').innerHTML = `<div class="fba-tablewrap" style="max-height:320px"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Min</th><th>Max</th><th>Span</th></tr></thead><tbody>${tr}</tbody></table></div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="fva-csv">⬇ Download FVA (CSV)</button>`;
  $('fva-csv').addEventListener('click', () => { let c = 'exchange_id,name,min,max,span\n'; rows.forEach(r => c += `${r.id},"${(r.name || '').replace(/"/g, '""')}",${r.min},${r.max},${r.span}\n`); saveCSV(c, `FVA_${fvaState.file.replace(/\.json.*/, '')}_${$('fva-media').value}`); });
}

// ── Dynamic FBA tab ────────────────────────────────────────────────────────────
const dfbaState = { model: null, file: null, series: null };
function initDFBA() {
  fillMedia($('dfba-media'), null);
  attachMedia('dfba-media', dfbaState);
  makeCombo($('dfba-model-input'), $('dfba-model-menu'), async (gemFile) => {
    $('dfba-model-input').value = '';
    const mc = $('dfba-modelcard'); mc.style.display = 'block'; mc.innerHTML = '<span class="fba-hint-inline">Loading…</span>';
    try {
      const model = await loadModel(gemFile); dfbaState.model = model; dfbaState.file = gemFile; mc.innerHTML = modelCardHTML(gemFile, model);
      // populate substrate options from carbon exchanges
      const sub = $('dfba-substrate'); const cur = sub.value;
      const exs = listExchanges(model).filter(e => /glc|glucose|fru|gal|lac|ac_e|succ|glyc|sucr|malt|xyl|arab/i.test(e.id));
      sub.innerHTML = exs.map(e => `<option value="${e.id}">${esc((e.name || e.id).replace(/ exchange$/i, ''))} (${e.id})</option>`).join('');
      const g = pickGlc([...sub.options].map(o => o.value));
      if (g) sub.value = g; else if (cur) sub.value = cur;
    } catch (e) { mc.innerHTML = `<span style="color:var(--bad)">${esc(e.message)}</span>`; }
  });
  $('dfba-run').addEventListener('click', runDFBAtab);
}
async function runDFBAtab() {
  if (!dfbaState.model) return setStatus('dfba-status', 'Choose a model first.', 'err');
  const media = mediaOf(dfbaState, 'dfba-media');
  const substrateEx = $('dfba-substrate').value;
  const opts = {
    substrateEx, substrate0: +$('dfba-s0').value, biomass0: +$('dfba-x0').value,
    vmax: +$('dfba-vmax').value, km: +$('dfba-km').value, dt: +$('dfba-dt').value, tmax: +$('dfba-tmax').value,
    trackEx: ['EX_ac_e', 'EX_for_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e'].filter(e => dfbaState.model.reactions.some(r => r.id === e)),
    onProgress: (i, n) => prog('dfba-prog', 'dfba-prog-bar', i / n),
  };
  $('dfba-run').disabled = true; $('dfba-results').style.display = 'none';
  setStatus('dfba-status', 'Simulating batch culture…', 'busy'); prog('dfba-prog', 'dfba-prog-bar', 0);
  try {
    const t0 = performance.now();
    const s = await runDFBA(dfbaState.model, media, opts);
    prog('dfba-prog', 'dfba-prog-bar', null);
    dfbaState.series = s;
    renderDFBA(s, substrateEx, dfbaState.model);
    const finalX = s.biomass[s.biomass.length - 1];
    setStatus('dfba-status', `Done — ${s.t.length} steps in ${((performance.now() - t0) / 1000).toFixed(1)} s. Final biomass ${fmt(finalX)} gDW/L.`, 'ok');
  } catch (e) { setStatus('dfba-status', 'Error: ' + e.message, 'err'); console.error(e); prog('dfba-prog', 'dfba-prog-bar', null); }
  finally { $('dfba-run').disabled = false; }
}
function renderDFBA(s, substrateEx, model) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  $('dfba-results').style.display = 'block';
  const traces = [{ x: s.t, y: s.biomass, name: 'Biomass (gDW/L)', yaxis: 'y2', line: { color: '#1a7f4b', width: 3 } }];
  const palette = ['#2c6fbb', '#c0392b', '#e08a1e', '#7d3c98', '#16a085'];
  let ci = 0;
  for (const [ex, arr] of Object.entries(s.conc)) {
    const isSub = ex === substrateEx;
    traces.push({ x: s.t, y: arr, name: (isSub ? '▸ ' : '') + bio(ex) + ' (mM)', line: { color: isSub ? '#333' : palette[ci++ % palette.length], width: isSub ? 3 : 2, dash: isSub ? 'solid' : 'solid' } });
  }
  window.Plotly.newPlot('dfba-plot', traces, {
    margin: { l: 55, r: 55, t: 10, b: 45 }, height: 460,
    xaxis: { title: 'Time (h)' },
    yaxis: { title: 'Concentration (mM)', rangemode: 'tozero' },
    yaxis2: { title: 'Biomass (gDW/L)', overlaying: 'y', side: 'right', rangemode: 'tozero', showgrid: false },
    legend: { orientation: 'h', y: 1.08 }, font: { size: 11 }, hovermode: 'x unified',
  }, { responsive: true, displaylogo: false });
  $('dfba-note').textContent = `Substrate ▸ ${bio(substrateEx)} depletes via Michaelis–Menten uptake; biomass (right axis, green) grows until the substrate runs out; fermentation products accumulate.`;
  $('dfba-csv').onclick = () => {
    const keys = Object.keys(s.conc);
    let csv = 'time_h,biomass_gDW_L,' + keys.map(k => bio(k) + '_mM').join(',') + '\n';
    for (let i = 0; i < s.t.length; i++) csv += `${s.t[i]},${s.biomass[i]},` + keys.map(k => s.conc[k][i]).join(',') + '\n';
    saveCSV(csv, `dFBA_${dfbaState.file.replace(/\.json.*/, '')}_${$('dfba-media').value}`);
  };
}

// ── Multi-model tab ────────────────────────────────────────────────────────────
const mm = { selected: [], results: null };
const mmState = { mediaSpec: null, openMinerals: false };
function initMulti() {
  fillMedia($('mm-media'), $('mm-media-desc'));
  mountMediaPicker($('mm-media'), mmState);
  makeCombo($('mm-model-input'), $('mm-model-menu'), (gemFile) => { addModel(gemFile); $('mm-model-input').value = ''; });
  document.querySelectorAll('.mm-quick button[data-add]').forEach(b => b.addEventListener('click', () => quickAdd(b.dataset.add)));
  $('mm-clear').addEventListener('click', () => { mm.selected = []; renderChips(); });
  $('mm-run').addEventListener('click', runMulti);
}
function addModel(gemFile) { if (!mm.selected.includes(gemFile) && window.gemBatchMap && window.gemBatchMap[gemFile]) mm.selected.push(gemFile); renderChips(); }
function quickAdd(kind) {
  const md = window.gemMetadata || [];
  if (kind === 'perphylo') {
    const byP = {}; md.forEach(m => { const p = m.phylogroup || 'Unknown'; (byP[p] = byP[p] || []).push(m.gem_file); });
    Object.values(byP).forEach(list => { for (let i = 0; i < Math.min(4, list.length); i++) addModel(list[Math.floor(i * list.length / 4)]); });
  } else {
    const n = +kind, pool = md.map(m => m.gem_file);
    for (let k = 0; k < n && mm.selected.length < 120; k++) addModel(pool[Math.floor(seededRand(mm.selected.length + k) * pool.length)]);
  }
}
let _seed = 12345; function seededRand(i) { _seed = (1103515245 * (_seed + i) + 12345) & 0x7fffffff; return _seed / 0x7fffffff; }
function renderChips() {
  $('mm-count').textContent = `${mm.selected.length} models selected`;
  $('mm-chips').innerHTML = mm.selected.map(f => { const m = meta(f); return `<span class="mm-chip" title="${esc(f)}">${esc((m.genome_name || m.strain || f).slice(0, 26))} <span class="x" data-f="${esc(f)}">✕</span></span>`; }).join('');
  $('mm-chips').querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { mm.selected = mm.selected.filter(f => f !== x.dataset.f); renderChips(); }));
}
async function runMulti() {
  if (mm.selected.length < 2) return setStatus('mm-status', 'Add at least 2 models.', 'err');
  const mmSpec = mmState.mediaSpec;   // may be a preset or any Media DB record; bound per model
  const meth = document.querySelector('input[name="mm-method"]:checked').value;
  $('mm-run').disabled = true; $('mm-results').style.display = 'none';
  setStatus('mm-status', `Running ${meth.toUpperCase()} on ${mm.selected.length} models…`, 'busy'); prog('mm-prog', 'mm-prog-bar', 0);
  try {
    const t0 = performance.now();
    const rows = [];
    for (let i = 0; i < mm.selected.length; i++) {
      const f = mm.selected[i];
      try {
        const model = await loadModel(f);
        const media = bindMedium(model, mmSpec.components, { openMinerals: mmState.openMinerals }).bounds;
        const fba = await runFBA(model, media);
        let res = fba;
        /* Parsimony is forced on a non-growing model whatever the method, because with
           the objective at zero every feasible vector is optimal and FBA returns an
           arbitrary one. Reading an exchange fingerprint off that would invent a
           secretion profile for a strain that is not growing at all, and on a minimal
           medium plenty of the Lactobacillaceae do not. */
        if (fba.optimal && (meth === 'pfba' || fba.growth <= 1e-9)) res = await runPFBA(model, media, fba);
        const ex = {}; for (const [id, v] of Object.entries(res.fluxes)) if (id.startsWith('EX_') && Math.abs(v) > 1e-6) ex[id] = v;
        rows.push({ file: f, meta: meta(f), growth: fba.optimal ? fba.growth : 0, ex });
      } catch (e) { rows.push({ file: f, meta: meta(f), growth: 0, ex: {}, error: e.message }); }
      setStatus('mm-status', `Solved ${i + 1}/${mm.selected.length}…`, 'busy'); prog('mm-prog', 'mm-prog-bar', (i + 1) / mm.selected.length);
    }
    prog('mm-prog', 'mm-prog-bar', null);
    mm.results = rows;
    renderMulti(rows);
    setStatus('mm-status', `Done — ${rows.length} models in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('mm-status', 'Error: ' + e.message, 'err'); console.error(e); prog('mm-prog', 'mm-prog-bar', null); }
  finally { $('mm-run').disabled = false; }
}
function renderMulti(rows) {
  $('mm-results').style.display = 'block';
  const feas = rows.filter(r => r.growth > 1e-9);
  const meanG = feas.length ? feas.reduce((s, r) => s + r.growth, 0) / feas.length : 0;
  $('mm-kpis').innerHTML =
    `<div class="fba-kpi"><div class="v">${rows.length}</div><div class="l">Models run</div></div>
     <div class="fba-kpi"><div class="v">${feas.length}</div><div class="l">Feasible (grow)</div></div>
     <div class="fba-kpi"><div class="v">${fmt(meanG)}</div><div class="l">Mean growth (h⁻¹)</div></div>
     <div class="fba-kpi"><div class="v">${new Set(rows.map(r => r.meta.phylogroup || 'Unknown')).size}</div><div class="l">Phylogroups</div></div>`;

  growthByPhylo(rows);
  buildAxisSelectors(rows);
  drawScatter(rows);
  drawPCA(rows);
  drawHeatmap(rows);
  $('mm-csv').onclick = () => downloadMultiCSV(rows);
  $('mm-detail').innerHTML = '';
}
function growthByPhylo(rows) {
  const groups = {};
  rows.forEach(r => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(r); });
  const traces = Object.entries(groups).map(([p, rs]) => ({
    type: 'box', boxpoints: 'all', jitter: 0.5, pointpos: 0, name: p, y: rs.map(r => r.growth),
    text: rs.map(r => r.meta.genome_name || r.file), customdata: rs.map(r => r.file),
    marker: { color: phyloColor(p), size: 6 }, line: { color: phyloColor(p) },
    hovertemplate: '%{text}<br>growth %{y:.3f}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-growth', traces, { margin: { l: 45, r: 15, t: 10, b: 35 }, height: 340, yaxis: { title: 'Growth (h⁻¹)', rangemode: 'tozero' }, showlegend: false, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}
function commonExchanges(rows, minCount = 2) {
  const count = {}; rows.forEach(r => Object.keys(r.ex).forEach(id => count[id] = (count[id] || 0) + 1));
  return Object.entries(count).filter(([, c]) => c >= minCount).map(([id]) => id).sort();
}
function buildAxisSelectors(rows) {
  const exs = commonExchanges(rows);
  const opts = [`<option value="__growth">Growth rate (h⁻¹)</option>`].concat(exs.map(id => `<option value="${id}">${esc(bio(id))} flux</option>`)).join('');
  $('mm-x').innerHTML = opts; $('mm-y').innerHTML = opts;
  $('mm-x').value = '__growth';
  $('mm-y').value = exs.includes('EX_ac_e') ? 'EX_ac_e' : (exs[0] || '__growth');
  $('mm-x').onchange = $('mm-y').onchange = () => drawScatter(rows);
}
function axisVal(r, key) { return key === '__growth' ? r.growth : (r.ex[key] || 0); }
function axisLabel(key) { return key === '__growth' ? 'Growth rate (h⁻¹)' : bio(key) + ' flux'; }
function drawScatter(rows) {
  const xk = $('mm-x').value, yk = $('mm-y').value;
  const groups = {};
  rows.forEach(r => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(r); });
  const traces = Object.entries(groups).map(([p, rs]) => ({
    type: 'scatter', mode: 'markers', name: p, x: rs.map(r => axisVal(r, xk)), y: rs.map(r => axisVal(r, yk)),
    text: rs.map(r => r.meta.genome_name || r.file), customdata: rs.map(r => r.file),
    marker: { color: phyloColor(p), size: 9, line: { color: '#fff', width: 0.5 } },
    hovertemplate: '%{text}<br>' + esc(axisLabel(xk)) + ' %{x:.3f}<br>' + esc(axisLabel(yk)) + ' %{y:.3f}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-scatter', traces, { margin: { l: 50, r: 15, t: 10, b: 45 }, height: 340, xaxis: { title: axisLabel(xk) }, yaxis: { title: axisLabel(yk) }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}

// PCA via Gram-matrix (n×n) eigendecomposition (Jacobi)
function drawPCA(rows) {
  const feats = commonExchanges(rows, 2);
  if (feats.length < 2 || rows.length < 3) { $('mm-plot-pca').innerHTML = '<div class="fba-hint-inline" style="padding:1rem">Need ≥3 models and ≥2 shared exchanges for PCA.</div>'; return; }
  const n = rows.length;
  // matrix X (n×p), center columns
  let X = rows.map(r => feats.map(f => r.ex[f] || 0));
  const mean = feats.map((_, j) => X.reduce((s, row) => s + row[j], 0) / n);
  X = X.map(row => row.map((v, j) => v - mean[j]));
  // Gram G = X Xᵀ (n×n)
  const G = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, k) => { let s = 0; for (let j = 0; j < feats.length; j++) s += X[i][j] * X[k][j]; return s; }));
  const { values, vectors } = jacobiEigen(G);
  const order = values.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
  const tot = values.reduce((s, v) => s + Math.max(0, v), 0) || 1;
  const pc = (rank) => { const [lam, idx] = order[rank]; const sq = Math.sqrt(Math.max(0, lam)); return { scores: vectors.map(row => row[idx] * sq), ev: Math.max(0, lam) / tot }; };
  const p1 = pc(0), p2 = pc(1);
  const groups = {};
  rows.forEach((r, i) => { const p = r.meta.phylogroup || 'Unknown'; (groups[p] = groups[p] || []).push(i); });
  const traces = Object.entries(groups).map(([p, idxs]) => ({
    type: 'scatter', mode: 'markers', name: p, x: idxs.map(i => p1.scores[i]), y: idxs.map(i => p2.scores[i]),
    text: idxs.map(i => rows[i].meta.genome_name || rows[i].file), customdata: idxs.map(i => rows[i].file),
    marker: { color: phyloColor(p), size: 9, line: { color: '#fff', width: 0.5 } },
    hovertemplate: '%{text}<extra>' + p + '</extra>',
  }));
  window.Plotly.newPlot('mm-plot-pca', traces, { margin: { l: 50, r: 15, t: 10, b: 45 }, height: 340,
    xaxis: { title: `PC1 (${(p1.ev * 100).toFixed(0)}%)` }, yaxis: { title: `PC2 (${(p2.ev * 100).toFixed(0)}%)` }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showDetail(ev.points[0].customdata)));
}
// Jacobi eigenvalue algorithm for symmetric matrices
function jacobiEigen(A) {
  const n = A.length; const a = A.map(r => r.slice());
  const V = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => i === j ? 1 : 0));
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
    if (off < 1e-12) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-14) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
      for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
      for (let k = 0; k < n; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
    }
  }
  return { values: a.map((r, i) => r[i]), vectors: V };
}
function drawHeatmap(rows) {
  const feats = commonExchanges(rows, 2);
  // most variable exchanges
  const varOf = (f) => { const vals = rows.map(r => r.ex[f] || 0); const m = vals.reduce((s, v) => s + v, 0) / vals.length; return vals.reduce((s, v) => s + (v - m) ** 2, 0); };
  const top = feats.map(f => [f, varOf(f)]).sort((a, b) => b[1] - a[1]).slice(0, 20).map(x => x[0]);
  const order = rows.map((r, i) => [r.growth, i]).sort((a, b) => b[0] - a[0]).map(x => x[1]);
  const z = order.map(i => top.map(f => rows[i].ex[f] || 0));
  window.Plotly.newPlot('mm-plot-heat', [{
    type: 'heatmap', z, x: top.map(bio), y: order.map(i => (rows[i].meta.genome_name || rows[i].file).slice(0, 22)),
    colorscale: 'RdBu', reversescale: true, zmid: 0, colorbar: { title: 'flux', thickness: 10 },
    hovertemplate: '%{y}<br>%{x}: %{z:.2f}<extra></extra>',
  }], { margin: { l: 150, r: 10, t: 10, b: 70 }, height: Math.max(300, order.length * 16 + 90), xaxis: { tickangle: -45, tickfont: { size: 9 } }, yaxis: { tickfont: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false });
}
function showDetail(file) {
  if (!file) return;
  const r = mm.results.find(x => x.file === file); if (!r) return;
  const m = r.meta; const rep = { sec: Object.entries(r.ex).filter(([, v]) => v > 1e-6).sort((a, b) => b[1] - a[1]).slice(0, 6) };
  $('mm-detail').innerHTML = `<div class="fba-note" style="background:var(--accent)">
    <strong>${esc(m.genome_name || m.strain || file)}</strong> <code>${esc(file)}</code>
    ${m.phylogroup ? '· Phylogroup ' + esc(m.phylogroup) : ''} · growth <strong>${fmt(r.growth)}</strong> h⁻¹
    · top secretions: ${rep.sec.map(([id, v]) => `${bio(id)} ${v.toFixed(1)}`).join(', ') || '—'}
    &nbsp;<a href="#" id="mm-open-explore">open in Explore →</a></div>`;
  $('mm-open-explore').onclick = (e) => { e.preventDefault(); if (window.fbaSetModel) { document.querySelector('#studio-tabs .studio-tab[data-tab="explore"]').click(); window.fbaSetModel('a', file); } };
}
function downloadMultiCSV(rows) {
  const feats = commonExchanges(rows, 1);
  let csv = 'gem_file,genome_name,phylogroup,growth_h-1,' + feats.map(bio).join(',') + '\n';
  rows.forEach(r => { csv += `${r.file},"${(r.meta.genome_name || '').replace(/"/g, '""')}",${r.meta.phylogroup || ''},${r.growth},` + feats.map(f => r.ex[f] || 0).join(',') + '\n'; });
  saveCSV(csv, `multimodel_${$('mm-media').value}`);
}

// ── Statistics (validated vs scipy) ───────────────────────────────────────────
function erf(x) { const t = 1 / (1 + 0.3275911 * Math.abs(x)); const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x); return x >= 0 ? y : -y; }
function normalCdf(z) { return 0.5 * (1 + erf(z / Math.SQRT2)); }
function mannWhitneyU(a, b) {
  const n1 = a.length, n2 = b.length, N = n1 + n2;
  if (!n1 || !n2) return { U: NaN, p: 1 };
  const all = a.map(v => [v, 0]).concat(b.map(v => [v, 1])).sort((x, y) => x[0] - y[0]);
  const ranks = new Array(N); let i = 0; const tie = [];
  while (i < N) { let j = i; while (j < N - 1 && all[j + 1][0] === all[i][0]) j++; const r = (i + j + 2) / 2; for (let k = i; k <= j; k++) ranks[k] = r; if (j > i) tie.push(j - i + 1); i = j + 1; }
  let R1 = 0; for (let k = 0; k < N; k++) if (all[k][1] === 0) R1 += ranks[k];
  const U1 = R1 - n1 * (n1 + 1) / 2, U = Math.min(U1, n1 * n2 - U1), mu = n1 * n2 / 2;
  const tieTerm = tie.reduce((s, t) => s + (t * t * t - t), 0);
  const sigma = Math.sqrt(n1 * n2 / 12 * ((N + 1) - tieTerm / (N * (N - 1))));
  if (!(sigma > 0)) return { U, p: 1 };
  return { U, p: Math.min(1, 2 * normalCdf((U - mu + 0.5) / sigma)) };
}
function benjaminiHochberg(ps) {
  const m = ps.length, idx = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]), q = new Array(m); let prev = 1;
  for (let k = m - 1; k >= 0; k--) { const [p, i] = idx[k]; const v = Math.min(prev, p * m / (k + 1)); q[i] = v; prev = v; }
  return q;
}

// ── Cohort comparison tab ──────────────────────────────────────────────────────
const COHORT_FIELDS = [['phylogroup', 'Phylogroup'], ['MLST', 'MLST (ST)'], ['pathovar', 'Pathotype'], ['isolation_source', 'Isolation source'], ['host_name', 'Host'], ['isolation_country', 'Country'], ['mash_cluster', 'MASH cluster'], ['serovar', 'Serovar'], ['oxygen_requirement', 'Oxygen requirement'], ['disease', 'Disease']];
const cohort = { a: { field: 'phylogroup', values: new Set() }, b: { field: 'phylogroup', values: new Set(), complement: false }, results: null };
const cohortState = { mediaSpec: null, openMinerals: false };
const clean = (v) => { const s = String(v == null ? '' : v).trim(); return (s === '' || s === 'nan' || s === 'None' || s === '-1') ? '' : s; };

function initCohort() {
  fillMedia($('cohort-media'), null);
  mountMediaPicker($('cohort-media'), cohortState);
  ['a', 'b'].forEach(c => setupCohortBuilder(c));
  $('cohort-b-complement').addEventListener('change', () => { $('cohort-b-manual').style.display = $('cohort-b-complement').checked ? 'none' : 'block'; updateCohortCount('b'); });
  $('cohort-run').addEventListener('click', runCohort);
}
function builderEl(c) { return document.querySelector(`.cohort-build[data-cohort="${c}"]`); }
function setupCohortBuilder(c) {
  const root = c === 'b' ? $('cohort-b-manual') : builderEl('a');
  const fieldSel = root.querySelector('.cohort-field');
  fieldSel.innerHTML = COHORT_FIELDS.map(([k, l]) => `<option value="${k}">${l}</option>`).join('');
  fieldSel.value = cohort[c].field;
  fieldSel.addEventListener('change', () => { cohort[c].field = fieldSel.value; cohort[c].values = new Set(); renderCohortValues(c); });
  root.querySelector('.cohort-valsearch').addEventListener('input', () => renderCohortValues(c));
  renderCohortValues(c);
}
function fieldCounts(field) {
  const c = {}; (window.gemMetadata || []).forEach(m => { const v = clean(m[field]); if (v) c[v] = (c[v] || 0) + 1; });
  return Object.entries(c).sort((a, b) => b[1] - a[1]);
}
function renderCohortValues(c) {
  const root = c === 'b' ? $('cohort-b-manual') : builderEl('a');
  const box = root.querySelector('.cohort-values');
  const q = root.querySelector('.cohort-valsearch').value.trim().toLowerCase();
  let vals = fieldCounts(cohort[c].field);
  if (q) vals = vals.filter(([v]) => v.toLowerCase().includes(q));
  box.innerHTML = vals.slice(0, 300).map(([v, n]) =>
    `<label><input type="checkbox" value="${esc(v)}" ${cohort[c].values.has(v) ? 'checked' : ''}> ${esc(v)} <span class="cnt">${n}</span></label>`).join('') || '<div class="fba-hint-inline" style="padding:4px">no values</div>';
  box.querySelectorAll('input').forEach(cb => cb.addEventListener('change', () => { cb.checked ? cohort[c].values.add(cb.value) : cohort[c].values.delete(cb.value); updateCohortCount(c); }));
  updateCohortCount(c);
}
function cohortFiles(c) {
  const md = window.gemMetadata || [];
  if (c === 'b' && $('cohort-b-complement').checked) { const aset = new Set(cohortFiles('a')); return md.filter(m => !aset.has(m.gem_file)).map(m => m.gem_file); }
  const vals = cohort[c].values;
  if (!vals.size) return [];
  return md.filter(m => vals.has(clean(m[cohort[c].field]))).map(m => m.gem_file);
}
function updateCohortCount(c) { $('cohort-' + c + '-count').textContent = cohortFiles(c).length + ' GEMs'; }

async function runCohort() {
  const aFiles = cohortFiles('a'), bFiles = cohortFiles('b');
  if (aFiles.length < 3 || bFiles.length < 3) return setStatus('cohort-status', 'Each cohort needs ≥3 GEMs. Broaden your selection.', 'err');
  const cap = Math.max(3, Math.min(80, +$('cohort-cap').value || 25));
  const cohortSpec = cohortState.mediaSpec;   // bound per model inside the loop
  const meth = document.querySelector('input[name="cohort-method"]:checked').value;
  const sampA = sample(aFiles, cap), sampB = sample(bFiles, cap);
  const overlap = new Set(sampA).size + new Set(sampB).size - new Set([...sampA, ...sampB]).size;
  const jobs = [...sampA.map(f => ['A', f]), ...sampB.map(f => ['B', f])];
  $('cohort-run').disabled = true; $('cohort-results').style.display = 'none';
  setStatus('cohort-status', `Running ${meth.toUpperCase()} on ${jobs.length} models…`, 'busy'); prog('cohort-prog', 'cohort-prog-bar', 0);
  try {
    const t0 = performance.now();
    const rows = [];
    for (let i = 0; i < jobs.length; i++) {
      const [grp, f] = jobs[i];
      try { const model = await loadModel(f); const media = bindMedium(model, cohortSpec.components, { openMinerals: cohortState.openMinerals }).bounds; const fba = await runFBA(model, media); let res = fba; if (fba.optimal && (meth === 'pfba' || fba.growth <= 1e-9)) res = await runPFBA(model, media, fba); const ex = {}; for (const [id, v] of Object.entries(res.fluxes)) if (id.startsWith('EX_') && Math.abs(v) > 1e-6) ex[id] = v; rows.push({ grp, file: f, meta: meta(f), growth: fba.optimal ? fba.growth : 0, ex }); }
      catch (e) { rows.push({ grp, file: f, meta: meta(f), growth: 0, ex: {} }); }
      prog('cohort-prog', 'cohort-prog-bar', (i + 1) / jobs.length); setStatus('cohort-status', `Solved ${i + 1}/${jobs.length}…`, 'busy');
    }
    prog('cohort-prog', 'cohort-prog-bar', null);
    cohort.results = { rows, aFiles, bFiles, sampA, sampB, overlap, meth };
    renderCohort(cohort.results);
    setStatus('cohort-status', `Done — ${jobs.length} models in ${((performance.now() - t0) / 1000).toFixed(1)} s.` + (overlap ? ` (${overlap} strain(s) in both cohorts)` : ''), 'ok');
  } catch (e) { setStatus('cohort-status', 'Error: ' + e.message, 'err'); console.error(e); prog('cohort-prog', 'cohort-prog-bar', null); }
  finally { $('cohort-run').disabled = false; }
}
function sample(arr, n) { if (arr.length <= n) return arr.slice(); const a = arr.slice(), out = []; for (let k = 0; k < n; k++) { const j = k + Math.floor(seededRand(k) * (a.length - k)); [a[k], a[j]] = [a[j], a[k]]; out.push(a[k]); } return out; }
function labA() { return cohortLabel('a'); }
function labB() { return $('cohort-b-complement').checked ? 'B (complement)' : cohortLabel('b'); }
function cohortLabel(c) { const f = COHORT_FIELDS.find(x => x[0] === cohort[c].field)[1]; return `${f}=${[...cohort[c].values].join('/') || '—'}`; }

function renderCohort(R) {
  const A = R.rows.filter(r => r.grp === 'A'), B = R.rows.filter(r => r.grp === 'B');
  const gA = A.map(r => r.growth), gB = B.map(r => r.growth);
  const feasA = gA.filter(g => g > 1e-9), feasB = gB.filter(g => g > 1e-9);
  const meanA = feasA.length ? feasA.reduce((s, v) => s + v, 0) / feasA.length : 0;
  const meanB = feasB.length ? feasB.reduce((s, v) => s + v, 0) / feasB.length : 0;
  const gTest = mannWhitneyU(gA, gB);
  $('cohort-results').style.display = 'block';
  $('cohort-kpis').innerHTML =
    `<div class="fba-kpi"><div class="v" style="color:var(--primary)">${A.length}</div><div class="l">Group A models</div></div>
     <div class="fba-kpi"><div class="v" style="color:var(--bad)">${B.length}</div><div class="l">Group B models</div></div>
     <div class="fba-kpi"><div class="v">${fmt(meanA)} / ${fmt(meanB)}</div><div class="l">Mean growth A / B</div></div>
     <div class="fba-kpi"><div class="v" style="color:${gTest.p < 0.05 ? 'var(--ok)' : 'var(--mute)'}">${gTest.p < 1e-4 ? gTest.p.toExponential(1) : gTest.p.toFixed(3)}</div><div class="l">Growth MWU p-value</div></div>`;

  // growth box
  window.Plotly.newPlot('cohort-plot-growth', [
    { type: 'box', boxpoints: 'all', jitter: 0.5, name: 'A: ' + labA(), y: gA, marker: { color: '#2c6fbb', size: 5 }, line: { color: '#2c6fbb' }, text: A.map(r => r.meta.genome_name || r.file) },
    { type: 'box', boxpoints: 'all', jitter: 0.5, name: 'B: ' + labB(), y: gB, marker: { color: '#c0392b', size: 5 }, line: { color: '#c0392b' }, text: B.map(r => r.meta.genome_name || r.file) },
  ], { margin: { l: 45, r: 10, t: 10, b: 30 }, height: 340, yaxis: { title: 'Growth (h⁻¹)', rangemode: 'tozero' }, showlegend: true, legend: { font: { size: 9 }, orientation: 'h', y: 1.12 }, font: { size: 11 } }, { responsive: true, displaylogo: false });

  // differential exchange fluxes
  const exIds = new Set(); R.rows.forEach(r => Object.keys(r.ex).forEach(id => exIds.add(id)));
  const diff = [];
  for (const id of exIds) {
    const va = A.map(r => r.ex[id] || 0), vb = B.map(r => r.ex[id] || 0);
    const nzA = va.filter(v => Math.abs(v) > 1e-9).length, nzB = vb.filter(v => Math.abs(v) > 1e-9).length;
    if (nzA + nzB < 2) continue;
    const mA = va.reduce((s, v) => s + v, 0) / va.length, mB = vb.reduce((s, v) => s + v, 0) / vb.length;
    const t = mannWhitneyU(va, vb);
    diff.push({ id, meanA: mA, meanB: mB, delta: mA - mB, p: t.p });
  }
  const qs = benjaminiHochberg(diff.map(d => d.p));
  diff.forEach((d, i) => d.q = qs[i]);
  diff.sort((a, b) => a.p - b.p);
  _cohortDiff = diff;
  renderVolcano(diff);
  renderCohortBars(diff);
  renderCohortPCA(R.rows);
  renderCohortTable(diff);
  $('cohort-detail').innerHTML = '';
  $('cohort-csv').onclick = () => { let c = 'exchange_id,metabolite,mean_flux_A,mean_flux_B,delta_A_minus_B,mwu_p,bh_q\n'; diff.forEach(d => { c += `${d.id},${bio(d.id)},${d.meanA},${d.meanB},${d.delta},${d.p},${d.q}\n`; }); saveCSV(c, `cohort_${$('cohort-media').value}`); };
}
function renderVolcano(diff) {
  const sig = diff.filter(d => d.q < 0.05), ns = diff.filter(d => !(d.q < 0.05));
  const mk = (arr, color) => ({ type: 'scatter', mode: 'markers', x: arr.map(d => d.delta), y: arr.map(d => -Math.log10(Math.max(d.p, 1e-300))), text: arr.map(d => bio(d.id)), customdata: arr.map(d => d.id), marker: { color, size: 8, line: { color: '#fff', width: 0.4 } }, hovertemplate: '%{text}<br>Δ(A−B) %{x:.3f}<br>-log10 p %{y:.2f}<extra></extra>' });
  window.Plotly.newPlot('cohort-plot-volcano', [
    Object.assign(mk(ns, '#b8c2cf'), { name: 'ns' }),
    Object.assign(mk(sig, '#7d3c98'), { name: 'FDR<0.05' }),
  ], { margin: { l: 45, r: 10, t: 10, b: 40 }, height: 340, xaxis: { title: 'Δ mean flux (A − B)', zeroline: true }, yaxis: { title: '−log₁₀ p' }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false })
    .then(gd => gd.on('plotly_click', ev => showCohortDetail(ev.points[0].customdata)));
}
function renderCohortBars(diff) {
  const top = diff.slice().sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta)).slice(0, 12).reverse();
  window.Plotly.newPlot('cohort-plot-bars', [
    { type: 'bar', orientation: 'h', name: 'A', y: top.map(d => bio(d.id)), x: top.map(d => d.meanA), marker: { color: '#2c6fbb' } },
    { type: 'bar', orientation: 'h', name: 'B', y: top.map(d => bio(d.id)), x: top.map(d => d.meanB), marker: { color: '#c0392b' } },
  ], { barmode: 'group', margin: { l: 70, r: 10, t: 10, b: 35 }, height: 340, xaxis: { title: 'Mean flux' }, legend: { font: { size: 9 } }, font: { size: 10 } }, { responsive: true, displaylogo: false });
}
function renderCohortPCA(rows) {
  const feats = commonExchanges(rows, 2);
  if (feats.length < 2 || rows.length < 4) { $('cohort-plot-pca').innerHTML = '<div class="fba-hint-inline" style="padding:1rem">Not enough data for PCA.</div>'; return; }
  const n = rows.length;
  let X = rows.map(r => feats.map(f => r.ex[f] || 0));
  const mean = feats.map((_, j) => X.reduce((s, row) => s + row[j], 0) / n);
  X = X.map(row => row.map((v, j) => v - mean[j]));
  const G = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, k) => { let s = 0; for (let j = 0; j < feats.length; j++) s += X[i][j] * X[k][j]; return s; }));
  const { values, vectors } = jacobiEigen(G);
  const order = values.map((v, i) => [v, i]).sort((a, b) => b[0] - a[0]);
  const tot = values.reduce((s, v) => s + Math.max(0, v), 0) || 1;
  const pc = (rank) => { const [lam, idx] = order[rank]; const sq = Math.sqrt(Math.max(0, lam)); return { s: vectors.map(row => row[idx] * sq), ev: Math.max(0, lam) / tot }; };
  const p1 = pc(0), p2 = pc(1);
  const grp = (g, color, name) => { const idxs = rows.map((r, i) => [r, i]).filter(([r]) => r.grp === g).map(([, i]) => i); return { type: 'scatter', mode: 'markers', name, x: idxs.map(i => p1.s[i]), y: idxs.map(i => p2.s[i]), text: idxs.map(i => rows[i].meta.genome_name || rows[i].file), customdata: idxs.map(i => rows[i].file), marker: { color, size: 9, line: { color: '#fff', width: 0.5 } }, hovertemplate: '%{text}<extra>' + name + '</extra>' }; };
  window.Plotly.newPlot('cohort-plot-pca', [grp('A', '#2c6fbb', 'Group A'), grp('B', '#c0392b', 'Group B')],
    { margin: { l: 45, r: 10, t: 10, b: 40 }, height: 340, xaxis: { title: `PC1 (${(p1.ev * 100).toFixed(0)}%)` }, yaxis: { title: `PC2 (${(p2.ev * 100).toFixed(0)}%)` }, legend: { font: { size: 9 } }, font: { size: 11 } }, { responsive: true, displaylogo: false });
}
function renderCohortTable(diff) {
  const rows = diff.slice(0, 40).map(d => `<tr><td><code>${esc(d.id)}</code></td><td>${esc(bio(d.id))}</td><td class="num">${fmt(d.meanA)}</td><td class="num">${fmt(d.meanB)}</td><td class="num" style="color:${d.delta >= 0 ? 'var(--primary)' : 'var(--bad)'}">${(d.delta >= 0 ? '+' : '') + fmt(d.delta)}</td><td class="num">${d.p < 1e-4 ? d.p.toExponential(1) : d.p.toFixed(4)}</td><td class="num" style="color:${d.q < 0.05 ? 'var(--ok)' : 'var(--mute)'}">${d.q < 1e-4 ? d.q.toExponential(1) : d.q.toFixed(4)}</td></tr>`).join('');
  $('cohort-table').innerHTML = `<div class="fba-tablewrap" style="max-height:320px"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Metabolite</th><th>Mean A</th><th>Mean B</th><th>Δ(A−B)</th><th>MWU p</th><th>BH q</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}
function showCohortDetail(exId) {
  const d = cohort.results && cohortDiffLookup(exId); if (!d) return;
  $('cohort-detail').innerHTML = `<div class="fba-note" style="background:var(--accent)"><strong>${esc(bio(exId))}</strong> <code>${esc(exId)}</code> — mean flux A ${fmt(d.meanA)}, B ${fmt(d.meanB)}, Δ ${fmt(d.delta)} · MWU p ${d.p.toExponential(2)} · BH q ${d.q.toExponential(2)} ${d.q < 0.05 ? '<strong style="color:var(--ok)">(significant)</strong>' : ''}</div>`;
}
let _cohortDiff = null;
function cohortDiffLookup(id) { return (_cohortDiff || []).find(d => d.id === id); }

// ── single-model picker helper (Envelope / PhasePlane / Essentiality) ──────────
function wireModelPicker(inputId, menuId, cardId, state, onLoad) {
  makeCombo($(inputId), $(menuId), async (gemFile) => {
    $(inputId).value = '';
    const mc = $(cardId); mc.style.display = 'block'; mc.innerHTML = '<span class="fba-hint-inline">Loading…</span>';
    try {
      const model = await loadModel(gemFile);
      state.model = model; state.file = gemFile; state.meta = meta(gemFile);
      mc.innerHTML = modelCardHTML(gemFile, model);
      if (state.onModelLoaded) state.onModelLoaded();   // rebind the medium to THIS model
      if (onLoad) onLoad(model);
    }
    catch (e) { mc.innerHTML = `<span style="color:var(--bad)">${esc(e.message)}</span>`; }
  });
}

// ── Production Envelope ─────────────────────────────────────────────────────────
const envState = { model: null, file: null };
function initEnvelope() {
  fillMedia($('env-media'), null);
  attachMedia('env-media', envState);
  wireModelPicker('env-model-input', 'env-model-menu', 'env-modelcard', envState, (model) => {
    const exs = listExchanges(model), has = id => exs.some(e => e.id === id);
    const pri = ['EX_ac_e', 'EX_etoh_e', 'EX_lac__D_e', 'EX_succ_e', 'EX_for_e', 'EX_ala__L_e', 'EX_pyr_e', 'EX_akg_e', 'EX_co2_e', 'EX_h2_e'].filter(has);
    const rest = exs.map(e => e.id).filter(id => !pri.includes(id) && id !== 'EX_h2o_e' && id !== 'EX_h_e');
    $('env-product').innerHTML = pri.concat(rest).map(id => `<option value="${id}">${bio(id)} (${id})</option>`).join('');
    if (has('EX_ac_e')) $('env-product').value = 'EX_ac_e';
  });
  $('env-run').addEventListener('click', runEnvelope);
}
async function runEnvelope() {
  if (!envState.model) return setStatus('env-status', 'Choose a model first.', 'err');
  const media = mediaOf(envState, 'env-media'), product = $('env-product').value;
  $('env-run').disabled = true; $('env-results').style.display = 'none';
  setStatus('env-status', 'Computing envelope…', 'busy'); prog('env-prog', 'env-prog-bar', 0);
  try {
    const t0 = performance.now();
    const res = await productionEnvelope(envState.model, media, product, { points: 24, onProgress: (d, t) => prog('env-prog', 'env-prog-bar', d / t) });
    prog('env-prog', 'env-prog-bar', null);
    $('env-results').style.display = 'block';
    const x = res.points.map(p => p.product), ymax = res.points.map(p => p.growthMax), ymin = res.points.map(p => p.growthMin);
    const peak = res.points.reduce((a, b) => b.growthMax > a.growthMax ? b : a, res.points[0]);
    window.Plotly.newPlot('env-plot', [
      { x, y: ymin, mode: 'lines', name: 'min biomass', line: { color: '#c0392b', width: 2 } },
      { x, y: ymax, mode: 'lines', name: 'max biomass', line: { color: '#2c6fbb', width: 3 }, fill: 'tonexty', fillcolor: 'rgba(44,111,187,0.12)' },
    ], {
      margin: { l: 55, r: 20, t: 30, b: 50 }, height: 460,
      xaxis: { title: `${bio(product)} production flux (mmol·gDW⁻¹·h⁻¹)` }, yaxis: { title: 'Growth rate (h⁻¹)', rangemode: 'tozero' },
      legend: { orientation: 'h', y: 1.12 }, font: { size: 11 },
      annotations: [{ x: peak.product, y: peak.growthMax, text: `peak ${fmt(peak.growthMax, 3)} @ ${fmt(peak.product, 1)}`, showarrow: true, arrowhead: 2, ax: 0, ay: -30, font: { size: 10 } }],
    }, { responsive: true, displaylogo: false });
    $('env-note').textContent = `The shaded region is the feasible biomass range at each production level. Max ${bio(product)} production is ${fmt(res.prodMax, 2)} mmol·gDW⁻¹·h⁻¹ (at zero growth).`;
    setStatus('env-status', `Done in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('env-status', 'Error: ' + e.message, 'err'); console.error(e); prog('env-prog', 'env-prog-bar', null); }
  finally { $('env-run').disabled = false; }
}

// ── Phenotype Phase Plane ───────────────────────────────────────────────────────
const ppState = { model: null, file: null };
function initPhasePlane() {
  fillMedia($('pp-media'), null);
  attachMedia('pp-media', ppState);
  wireModelPicker('pp-model-input', 'pp-model-menu', 'pp-modelcard', ppState, (model) => {
    const exs = listExchanges(model), opts = exs.map(e => `<option value="${e.id}">${bio(e.id)} (${e.id})</option>`).join('');
    $('pp-x').innerHTML = opts; $('pp-y').innerHTML = opts;
    const gx = pickGlc(exs.map(e => e.id));
    if (gx) $('pp-x').value = gx;
    if (exs.some(e => e.id === 'EX_o2_e')) $('pp-y').value = 'EX_o2_e';
  });
  $('pp-run').addEventListener('click', runPP);
}
async function runPP() {
  if (!ppState.model) return setStatus('pp-status', 'Choose a model first.', 'err');
  const media = mediaOf(ppState, 'pp-media'), xId = $('pp-x').value, yId = $('pp-y').value, n = Math.max(6, Math.min(40, +$('pp-n').value || 20));
  if (xId === yId) return setStatus('pp-status', 'Choose two different axes.', 'err');
  $('pp-run').disabled = true; $('pp-results').style.display = 'none';
  setStatus('pp-status', `Solving ${n}×${n} grid…`, 'busy'); prog('pp-prog', 'pp-prog-bar', 0);
  try {
    const t0 = performance.now();
    const res = await phasePlane(ppState.model, media, xId, yId, { n, xMax: 20, yMax: 20, onProgress: (d, t) => prog('pp-prog', 'pp-prog-bar', d / t) });
    prog('pp-prog', 'pp-prog-bar', null);
    $('pp-results').style.display = 'block';
    window.Plotly.newPlot('pp-plot', [{ type: 'contour', z: res.Z, x: res.xs, y: res.ys, colorscale: 'Viridis', contours: { coloring: 'heatmap' }, colorbar: { title: 'growth h⁻¹', thickness: 12 }, hovertemplate: `${bio(xId)} %{x:.1f}<br>${bio(yId)} %{y:.1f}<br>growth %{z:.3f}<extra></extra>` }],
      { margin: { l: 60, r: 10, t: 15, b: 55 }, height: 480, xaxis: { title: `${bio(xId)} uptake capacity` }, yaxis: { title: `${bio(yId)} uptake capacity` }, font: { size: 11 } }, { responsive: true, displaylogo: false });
    setStatus('pp-status', `Done — ${n * n} solves in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('pp-status', 'Error: ' + e.message, 'err'); console.error(e); prog('pp-prog', 'pp-prog-bar', null); }
  finally { $('pp-run').disabled = false; }
}

// ── Essentiality Screen ─────────────────────────────────────────────────────────
const essState = { model: null, file: null };
function initEssential() {
  fillMedia($('ess-media'), null);
  attachMedia('ess-media', essState);
  wireModelPicker('ess-model-input', 'ess-model-menu', 'ess-modelcard', essState);
  $('ess-run').addEventListener('click', runEss);
}
function essClass(ratio) { return ratio < 0.01 ? 'Essential' : ratio < 0.5 ? 'Severe' : ratio < 0.95 ? 'Mild' : 'Neutral'; }
const ESS_COLORS = { Essential: '#c0392b', Severe: '#e08a1e', Mild: '#5b8ff9', Neutral: '#cfd8e3' };
async function runEss() {
  if (!essState.model) return setStatus('ess-status', 'Choose a model first.', 'err');
  const media = mediaOf(essState, 'ess-media'), scope = $('ess-scope').value;
  const all = essState.model.reactions.map(r => r.id);
  const ids = scope === 'exchange' ? all.filter(id => id.startsWith('EX_')) : scope === 'metabolic' ? all.filter(id => !id.startsWith('EX_')) : all;
  $('ess-run').disabled = true; $('ess-results').style.display = 'none';
  setStatus('ess-status', `Knocking out ${ids.length} reactions… (this can take a minute)`, 'busy'); prog('ess-prog', 'ess-prog-bar', 0);
  try {
    const t0 = performance.now();
    const res = await essentialityScan(essState.model, media, ids, { onProgress: (d, t) => { if (d % 20 === 0 || d === t) { prog('ess-prog', 'ess-prog-bar', d / t); setStatus('ess-status', `Tested ${d}/${t}…`, 'busy'); } } });
    prog('ess-prog', 'ess-prog-bar', null);
    renderEss(res, ids.length);
    setStatus('ess-status', `Done — ${ids.length} knockouts in ${((performance.now() - t0) / 1000).toFixed(1)} s. WT growth ${fmt(res.wtGrowth)} h⁻¹.`, 'ok');
  } catch (e) { setStatus('ess-status', 'Error: ' + e.message, 'err'); console.error(e); prog('ess-prog', 'ess-prog-bar', null); }
  finally { $('ess-run').disabled = false; }
}
function renderEss(res, n) {
  const rows = res.results.map(r => ({ ...r, cls: essClass(r.ratio) })).sort((a, b) => a.ratio - b.ratio);
  const counts = { Essential: 0, Severe: 0, Mild: 0, Neutral: 0 };
  rows.forEach(r => counts[r.cls]++);
  $('ess-results').style.display = 'block';
  $('ess-kpis').innerHTML =
    `<div class="fba-kpi"><div class="v" style="color:var(--bad)">${counts.Essential}</div><div class="l">Essential (&lt;1%)</div></div>
     <div class="fba-kpi"><div class="v" style="color:var(--warn)">${counts.Severe}</div><div class="l">Severe (&lt;50%)</div></div>
     <div class="fba-kpi"><div class="v" style="color:var(--primary)">${counts.Mild}</div><div class="l">Mild (&lt;95%)</div></div>
     <div class="fba-kpi"><div class="v">${counts.Neutral}</div><div class="l">Dispensable</div></div>`;
  window.Plotly.newPlot('ess-plot-pie', [{ type: 'pie', labels: Object.keys(counts), values: Object.values(counts), marker: { colors: Object.keys(counts).map(k => ESS_COLORS[k]) }, textinfo: 'label+percent', hole: 0.45 }],
    { margin: { l: 10, r: 10, t: 10, b: 10 }, height: 320, font: { size: 11 }, showlegend: false }, { responsive: true, displaylogo: false });
  window.Plotly.newPlot('ess-plot-hist', [{ type: 'histogram', x: rows.map(r => r.ratio), nbinsx: 25, marker: { color: '#2c6fbb' } }],
    { margin: { l: 45, r: 10, t: 10, b: 40 }, height: 320, xaxis: { title: 'Knockout growth / WT growth' }, yaxis: { title: 'reactions' }, font: { size: 11 } }, { responsive: true, displaylogo: false });
  const disp = rows.filter(r => r.ratio < 0.95).slice(0, 60);
  $('ess-table').innerHTML = `<div class="fba-hint-inline" style="margin-bottom:4px">Reactions with a growth impact (ratio &lt; 0.95), most severe first — top ${disp.length}:</div>
    <div class="fba-tablewrap" style="max-height:320px"><table class="fba-flux"><thead><tr><th>Reaction</th><th>Name</th><th>Subsystem</th><th>KO growth</th><th>Ratio</th><th>Class</th></tr></thead>
    <tbody>${disp.map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td><td>${esc(r.subsystem)}</td><td class="num">${fmt(r.growth)}</td><td class="num" style="color:${ESS_COLORS[r.cls]}">${fmt(r.ratio, 3)}</td><td>${r.cls}</td></tr>`).join('')}</tbody></table></div>`;
  $('ess-csv').onclick = () => { let c = 'reaction_id,name,subsystem,ko_growth,ratio,class\n'; rows.forEach(r => c += `${r.id},"${(r.name || '').replace(/"/g, '""')}","${(r.subsystem || '').replace(/"/g, '""')}",${r.growth},${r.ratio},${r.cls}\n`); saveCSV(c, `essentiality_${essState.file.replace(/\.json.*/, '')}_${$('ess-media').value}`); };
}

// ── init ──────────────────────────────────────────────────────────────────────
async function init() {
  await presets();
  await mediaPresets();      // the shared picker keeps its own copy; load it before mounting
  initNav(); initFVA(); initDFBA(); initMulti(); initCohort(); initEnvelope(); initPhasePlane(); initEssential();
  initGenes(); initSyn(); initFseof(); initSamp(); initQC();
  wireDownloads(); decorateStaticPlots();

  const q = new URLSearchParams(location.search);
  if (q.get('model') || q.get('models') || q.get('tab')) {
    const apply = () => applyDeepLink(q);
    // the pickers cannot resolve a gem_file until the model index has landed
    if (window.gemBatchMap) apply(); else window.addEventListener('gem-data-ready', apply, { once: true });
  }
}

async function applyDeepLink(q) {
  const models = q.get('models'), model = q.get('model'), medium = q.get('medium');
  const ko = (q.get('ko') || '').split(',').map(s => s.trim()).filter(Boolean);
  let tab = q.get('tab');

  if (models) {
    models.split(',').map(s => s.trim()).filter(Boolean).forEach(addModel);
    switchTab('multi'); setMediumOn('multi', medium);
    return;
  }
  if (!model) { if (tab) switchTab(tab); return; }

  // a knockout only means anything in Explore, which is where the knockout study lives
  if (ko.length) tab = 'explore';
  if (!tab || tab === 'explore' || !MODEL_LOADERS[tab]) tab = 'explore';
  switchTab(tab);

  if (tab === 'explore') {
    if (window.fbaSetModel) await window.fbaSetModel('a', model, { medium, ko });
    return;
  }
  /* Set the medium BEFORE the model arrives. The picker swaps in a collection default
     when a model lands unless the user has chosen one, and a link IS a choice. */
  setMediumOn(tab, medium);
  await MODEL_LOADERS[tab](model);
}

function setMediumOn(tab, key) {
  if (!key) return;
  const sel = $(MEDIA_SELECT[tab]);
  if (!sel || !PRESETS || !PRESETS[key]) return;      // unknown medium: keep the default
  sel.value = key;
  sel.dispatchEvent(new Event('change'));             // the shared picker rebinds on change
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();


/* ════════════════════════ Gene essentiality ═══════════════════════════════ */
const genesState = { model: null, file: null, res: null };
function initGenes() {
  fillMedia($('genes-media'), $('genes-media-desc'));
  attachMedia('genes-media', genesState);
  wireModelPicker('genes-model-input', 'genes-model-menu', 'genes-modelcard', genesState);
  $('genes-run').addEventListener('click', runGenes);
}
const geneClass = (r) => r < 0.01 ? 'Essential' : r < 0.5 ? 'Severe' : r < 0.95 ? 'Mild' : 'Dispensable';
const CLS_COL = { Essential: '#c0392b', Severe: '#e08a1e', Mild: '#5b8ff9', Dispensable: '#9aa7b8' };

async function runGenes() {
  if (!genesState.model) return setStatus('genes-status', 'Choose a model first.', 'err');
  const media = mediaOf(genesState, 'genes-media');
  const cap = +$('genes-scope').value;
  const all = listGenes(genesState.model).map(g => g.id);
  const genes = cap ? all.slice(0, cap) : all;
  const btn = $('genes-run'); btn.disabled = true;
  setStatus('genes-status', `Deleting ${genes.length} genes through their GPRs…`, 'busy');
  try {
    const t0 = performance.now();
    const res = await singleGeneDeletion(genesState.model, media, genes,
      { onProgress: (d, t) => prog('genes-prog', 'genes-prog-bar', d / t) });
    prog('genes-prog', 'genes-prog-bar', null);
    genesState.res = res;
    renderGenes(res);
    setStatus('genes-status', `Done — ${genes.length} single-gene deletions in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('genes-status', 'Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function renderGenes(res) {
  const rows = res.results.map(r => ({ ...r, cls: geneClass(r.ratio) }));
  const counts = {}; rows.forEach(r => counts[r.cls] = (counts[r.cls] || 0) + 1);
  $('genes-results').style.display = 'block';
  $('genes-kpis').innerHTML =
    kpi(fmt(res.wtGrowth), 'Wild-type growth (h⁻¹)') +
    kpi(counts.Essential || 0, 'Essential genes', '#c0392b') +
    kpi(counts.Severe || 0, 'Severe', '#e08a1e') +
    kpi(counts.Mild || 0, 'Mild', '#5b8ff9') +
    kpi(counts.Dispensable || 0, 'Dispensable', '#9aa7b8') +
    kpi(rows.length, 'Genes tested');
  const order = ['Essential', 'Severe', 'Mild', 'Dispensable'];
  Plotly.newPlot('genes-plot-pie', [{
    type: 'pie', hole: 0.55, labels: order, values: order.map(c => counts[c] || 0),
    marker: { colors: order.map(c => CLS_COL[c]) }, textinfo: 'label+percent',
  }], plotly('', '', ''), PLOT_CFG);
  Plotly.newPlot('genes-plot-hist', [{
    type: 'histogram', x: rows.map(r => r.ratio), nbinsx: 30, marker: { color: '#2c6fbb' },
  }], plotly('', 'Growth ratio (KO / wild type)', 'Genes'), PLOT_CFG);

  const ess = rows.filter(r => r.cls !== 'Dispensable').sort((a, b) => a.ratio - b.ratio);
  $('genes-table').innerHTML = `
    <h6>Genes that matter (${ess.length})</h6>
    <div class="fba-tablewrap"><table class="fba-flux">
      <thead><tr><th>Gene</th><th>Name</th><th>Reactions switched off</th><th>Growth</th><th>Ratio</th><th>Class</th></tr></thead>
      <tbody>${ess.slice(0, 300).map(r => `
        <tr><td><code>${esc(r.gene)}</code></td><td>${esc(r.name)}</td>
        <td>${r.nOff ? `<code>${r.rxns.slice(0, 3).map(esc).join('</code> <code>')}</code>${r.rxns.length > 3 ? ` +${r.rxns.length - 3}` : ''}` : '<span style="color:var(--mute)">none</span>'}</td>
        <td class="num">${fmt(r.growth)}</td><td class="num">${(100 * r.ratio).toFixed(1)}%</td>
        <td><span class="fba-badge" style="background:${CLS_COL[r.cls]}">${r.cls}</span></td></tr>`).join('')}</tbody>
    </table></div>`;
  $('genes-csv').onclick = () => {
    let c = 'gene,name,n_reactions_off,reactions_off,growth,ratio,class\n';
    rows.forEach(r => c += `${r.gene},"${(r.name || '').replace(/"/g, '""')}",${r.nOff},"${r.rxns.join(';')}",${r.growth},${r.ratio},${r.cls}\n`);
    saveCSV(c, `gene_essentiality_${genesState.file.replace(/\.json.*/, '')}`);
  };
}

/* ════════════════════════ Synthetic lethality ═════════════════════════════ */
const synState = { model: null, file: null, res: null };
function initSyn() {
  fillMedia($('syn-media'), $('syn-media-desc'));
  attachMedia('syn-media', synState);
  wireModelPicker('syn-model-input', 'syn-model-menu', 'syn-modelcard', synState);
  $('syn-run').addEventListener('click', runSyn);
}

async function runSyn() {
  if (!synState.model) return setStatus('syn-status', 'Choose a model first.', 'err');
  const model = synState.model;
  const media = mediaOf(synState, 'syn-media');
  const kind = $('syn-kind').value, n = +$('syn-n').value;
  const btn = $('syn-run'); btn.disabled = true;
  setStatus('syn-status', 'Ranking candidates by single-deletion impact…', 'busy');
  try {
    const t0 = performance.now();
    // candidates: the reactions/genes that hurt most on their own but are not already lethal
    let ids;
    if (kind === 'gene') {
      const gs = listGenes(model).map(g => g.id).slice(0, 250);
      const sg = await singleGeneDeletion(model, media, gs, { onProgress: (d, t) => prog('syn-prog', 'syn-prog-bar', 0.5 * d / t) });
      ids = sg.results.filter(r => r.ratio > 0.05).sort((a, b) => a.ratio - b.ratio).slice(0, n).map(r => r.gene);
    } else {
      const rx = model.reactions.filter(r => !r.id.startsWith('EX_') && !/BIOMASS/i.test(r.id)).map(r => r.id).slice(0, 250);
      const sr = await essentialityScan(model, media, rx, { onProgress: (d, t) => prog('syn-prog', 'syn-prog-bar', 0.5 * d / t) });
      ids = sr.results.filter(r => r.ratio > 0.05).sort((a, b) => a.ratio - b.ratio).slice(0, n).map(r => r.id);
    }
    setStatus('syn-status', `Testing ${(ids.length * (ids.length - 1)) / 2} pairs…`, 'busy');
    const res = await doubleDeletion(model, media, ids, kind,
      { onProgress: (d, t) => prog('syn-prog', 'syn-prog-bar', 0.5 + 0.5 * d / t) });
    prog('syn-prog', 'syn-prog-bar', null);
    synState.res = res;
    renderSyn(res);
    const sl = res.pairs.filter(p => p.synthetic).length;
    setStatus('syn-status', `Done — ${res.pairs.length} pairs in ${((performance.now() - t0) / 1000).toFixed(1)} s. ${sl} synthetic lethal.`, 'ok');
  } catch (e) { setStatus('syn-status', 'Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function renderSyn(res) {
  const ids = res.ids, n = ids.length;
  const Z = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) =>
    i === j ? res.singles[ids[i]] : null));
  res.pairs.forEach(p => {
    const i = ids.indexOf(p.a), j = ids.indexOf(p.b);
    Z[i][j] = p.ratio; Z[j][i] = p.ratio;
  });
  const sl = res.pairs.filter(p => p.synthetic);
  $('syn-results').style.display = 'block';
  $('syn-kpis').innerHTML =
    kpi(fmt(res.wtGrowth), 'Wild-type growth (h⁻¹)') +
    kpi(res.pairs.length, 'Pairs tested') +
    kpi(sl.length, 'Synthetic lethal', sl.length ? '#c0392b' : '#1a7f4b') +
    kpi(res.pairs.filter(p => p.ratio < 0.5 && !p.synthetic).length, 'Severe pairs', '#e08a1e') +
    kpi(res.kind === 'gene' ? 'genes' : 'reactions', 'Deleting');
  Plotly.newPlot('syn-plot', [{
    type: 'heatmap', z: Z, x: ids, y: ids, zmin: 0, zmax: 1,
    colorscale: [[0, '#c0392b'], [0.5, '#f0c674'], [1, '#1a7f4b']],
    colorbar: { title: 'growth / WT' },
    hovertemplate: '%{y} + %{x}<br>growth ratio %{z:.3f}<extra></extra>',
  }], plotly('', '', ''), PLOT_CFG);
  $('syn-table').innerHTML = `
    <h6>Synthetic lethal pairs (${sl.length})</h6>
    ${sl.length ? `<div class="fba-tablewrap"><table class="fba-flux">
      <thead><tr><th>A</th><th>B</th><th>A alone</th><th>B alone</th><th>Together</th></tr></thead>
      <tbody>${sl.map(p => `<tr>
        <td><code>${esc(p.a)}</code></td><td><code>${esc(p.b)}</code></td>
        <td class="num">${(100 * p.ra).toFixed(0)}%</td><td class="num">${(100 * p.rb).toFixed(0)}%</td>
        <td class="num" style="color:var(--bad)">${(100 * p.ratio).toFixed(1)}%</td></tr>`).join('')}</tbody>
    </table></div>` : `<div class="fba-note">No synthetic lethal pair in this candidate set: every pair that kills has a member that already kills on its own.</div>`}`;
  $('syn-csv').onclick = () => {
    let c = 'a,b,growth_a_alone,growth_b_alone,growth_together,synthetic_lethal\n';
    res.pairs.forEach(p => c += `${p.a},${p.b},${p.ra},${p.rb},${p.ratio},${p.synthetic}\n`);
    saveCSV(c, `synthetic_lethality_${synState.file.replace(/\.json.*/, '')}`);
  };
}

/* ═══════════════════════════ FSEOF strain design ══════════════════════════ */
const fseofState = { model: null, file: null, res: null };
function initFseof() {
  fillMedia($('fseof-media'), $('fseof-media-desc'));
  attachMedia('fseof-media', fseofState);
  wireModelPicker('fseof-model-input', 'fseof-model-menu', 'fseof-modelcard', fseofState, (model) => {
    const sel = $('fseof-product'); sel.innerHTML = '';
    listExchanges(model).forEach(e => {
      const o = document.createElement('option'); o.value = e.id;
      o.textContent = `${e.id.replace(/^EX_/, '').replace(/_e$/, '')} — ${e.name.replace(/ exchange$/i, '')}`;
      sel.appendChild(o);
    });
    const pref = ['EX_succ_e', 'EX_ac_e', 'EX_lac__D_e', 'EX_lac_D_e', 'EX_etoh_e', 'EX_for_e'];
    const hit = pref.find(x => [...sel.options].some(o => o.value === x));
    if (hit) sel.value = hit;
  });
  $('fseof-run').addEventListener('click', runFseof);
}

async function runFseof() {
  if (!fseofState.model) return setStatus('fseof-status', 'Choose a model first.', 'err');
  const media = mediaOf(fseofState, 'fseof-media');
  const product = $('fseof-product').value;
  const steps = Math.max(4, Math.min(20, +$('fseof-steps').value || 10));
  const btn = $('fseof-run'); btn.disabled = true;
  setStatus('fseof-status', `Scanning ${steps} enforced levels of ${product}…`, 'busy');
  try {
    const t0 = performance.now();
    const res = await fseof(fseofState.model, media, product,
      { steps, onProgress: (d, t) => prog('fseof-prog', 'fseof-prog-bar', d / t) });
    prog('fseof-prog', 'fseof-prog-bar', null);
    if (!res.prodMax) { setStatus('fseof-status', 'This model cannot secrete that product on this medium.', 'err'); btn.disabled = false; return; }
    fseofState.res = res;
    renderFseof(res);
    setStatus('fseof-status', `Done — ${res.targets.length} amplification targets in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('fseof-status', 'Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function renderFseof(res) {
  const top = res.targets.slice(0, 12);
  $('fseof-results').style.display = 'block';
  $('fseof-kpis').innerHTML =
    kpi(fmt(res.prodMax, 2), 'Max product flux') +
    kpi(res.targets.length, 'Amplification targets', '#1a7f4b') +
    kpi(res.levels.length, 'Enforced levels') +
    kpi(res.levels.length ? fmt(res.levels[0].growth) : '—', 'Growth at lowest enforcement') +
    kpi(res.levels.length ? fmt(res.levels[res.levels.length - 1].growth) : '—', 'Growth at highest');
  Plotly.newPlot('fseof-plot-tradeoff', [{
    x: res.levels.map(l => l.enforced), y: res.levels.map(l => l.growth),
    mode: 'lines+markers', line: { color: '#2c6fbb', width: 3 },
  }], plotly('', `Enforced ${res.productEx} flux`, 'Growth (h⁻¹)'), PLOT_CFG);
  Plotly.newPlot('fseof-plot-targets',
    top.map((t, i) => ({
      x: res.levels.map(l => l.enforced), y: t.series, mode: 'lines', name: t.id,
      line: { width: 2, color: PALETTE[i % PALETTE.length] },
    })), plotly('', `Enforced ${res.productEx} flux`, 'Reaction flux'), PLOT_CFG);
  $('fseof-table').innerHTML = `
    <h6>Over-expression targets (${res.targets.length})</h6>
    <div class="fba-tablewrap"><table class="fba-flux">
      <thead><tr><th>Reaction</th><th>Name</th><th>Subsystem</th><th>Wild type</th><th>At max enforcement</th><th>Fold</th></tr></thead>
      <tbody>${res.targets.slice(0, 200).map(t => `<tr>
        <td><code>${esc(t.id)}</code></td><td>${esc(t.name)}</td><td>${esc(t.subsystem)}</td>
        <td class="num">${fmt(t.first, 3)}</td><td class="num">${fmt(t.last, 3)}</td>
        <td class="num" style="color:var(--ok)">${isFinite(t.fold) ? '×' + t.fold.toFixed(1) : 'de novo'}</td></tr>`).join('')}</tbody>
    </table></div>`;
  $('fseof-csv').onclick = () => {
    let c = 'reaction,name,subsystem,wild_type_flux,flux_at_max_enforcement,fold\n';
    res.targets.forEach(t => c += `${t.id},"${(t.name || '').replace(/"/g, '""')}","${t.subsystem}",${t.first},${t.last},${t.fold}\n`);
    saveCSV(c, `FSEOF_${fseofState.file.replace(/\.json.*/, '')}_${res.productEx}`);
  };
}

/* ═════════════════════════════ Flux sampling ══════════════════════════════ */
const sampState = { model: null, file: null, res: null };
function initSamp() {
  fillMedia($('samp-media'), $('samp-media-desc'));
  attachMedia('samp-media', sampState);
  wireModelPicker('samp-model-input', 'samp-model-menu', 'samp-modelcard', sampState);
  $('samp-run').addEventListener('click', runSamp);
}

async function runSamp() {
  if (!sampState.model) return setStatus('samp-status', 'Choose a model first.', 'err');
  const media = mediaOf(sampState, 'samp-media');
  const n = +$('samp-n').value, frac = +$('samp-frac').value;
  const btn = $('samp-run'); btn.disabled = true;
  setStatus('samp-status', 'Building warm-up points, then walking the polytope…', 'busy');
  try {
    const t0 = performance.now();
    const res = await sampleFluxes(sampState.model, media,
      { samples: n, warmup: 120, fraction: frac, onProgress: (d, t) => prog('samp-prog', 'samp-prog-bar', d / t) });
    prog('samp-prog', 'samp-prog-bar', null);
    if (!res.ok) { setStatus('samp-status', 'Cannot sample: ' + res.reason, 'err'); btn.disabled = false; return; }
    sampState.res = res;
    renderSamp(res);
    setStatus('samp-status', `Done — ${res.samples.length} samples from ${res.warmup} warm-up points in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('samp-status', 'Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function renderSamp(res) {
  const active = res.track.filter(id => res.samples.some(s => Math.abs(s[id]) > 1e-6)).slice(0, 14);
  $('samp-results').style.display = 'block';
  $('samp-kpis').innerHTML =
    kpi(res.samples.length, 'Samples') +
    kpi(res.warmup, 'Warm-up points') +
    kpi((100 * res.fraction).toFixed(0) + '%', 'Growth held at') +
    kpi(active.length, 'Exchanges with spread');
  Plotly.newPlot('samp-plot', active.map((id, i) => ({
    type: 'violin', y: res.samples.map(s => s[id]), name: id.replace(/^EX_/, '').replace(/_e$/, ''),
    box: { visible: true }, meanline: { visible: true },
    line: { color: PALETTE[i % PALETTE.length] },
  })), plotly('', '', 'Flux (mmol gDW⁻¹ h⁻¹)'), PLOT_CFG);
  $('samp-csv').onclick = () => {
    let c = res.track.join(',') + '\n';
    res.samples.forEach(s => c += res.track.map(id => s[id]).join(',') + '\n');
    saveCSV(c, `flux_samples_${sampState.file.replace(/\.json.*/, '')}`);
  };
}

/* ════════════════════════════════ Model QC ════════════════════════════════ */
const qcState = { model: null, file: null, res: null };
function initQC() {
  fillMedia($('qc-media'), $('qc-media-desc'));
  attachMedia('qc-media', qcState);
  wireModelPicker('qc-model-input', 'qc-model-menu', 'qc-modelcard', qcState);
  $('qc-run').addEventListener('click', runQC);
}

async function runQC() {
  if (!qcState.model) return setStatus('qc-status', 'Choose a model first.', 'err');
  const model = qcState.model;
  const media = mediaOf(qcState, 'qc-media');
  const cap = +$('qc-scope').value;
  const btn = $('qc-run'); btn.disabled = true;
  setStatus('qc-status', 'Checking mass and charge balance…', 'busy');
  try {
    const t0 = performance.now();
    const q = modelQC(model);                                   // structural, no LP
    const ids = model.reactions.map(r => r.id);
    const scan = cap ? ids.slice(0, cap) : ids;
    setStatus('qc-status', `Scanning ${scan.length} reactions for blocked flux…`, 'busy');
    const bl = await findBlockedReactions(model, media, { reactionIds: scan,
      onProgress: (d, t) => prog('qc-prog', 'qc-prog-bar', d / t) });
    prog('qc-prog', 'qc-prog-bar', null);
    qcState.res = { q, bl };
    renderQC(q, bl);
    setStatus('qc-status', `Done in ${((performance.now() - t0) / 1000).toFixed(1)} s.`, 'ok');
  } catch (e) { setStatus('qc-status', 'Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function renderQC(q, bl) {
  $('qc-results').style.display = 'block';
  $('qc-kpis').innerHTML =
    kpi(bl.blocked.length, `Blocked (of ${bl.tested})`, bl.blocked.length ? '#e08a1e' : '#1a7f4b') +
    kpi(q.massImbalanced.length, 'Mass imbalanced', q.massImbalanced.length ? '#c0392b' : '#1a7f4b') +
    kpi(q.chargeImbalanced.length, 'Charge imbalanced', q.chargeImbalanced.length ? '#c0392b' : '#1a7f4b') +
    kpi(q.deadEnds.length, 'Dead-end metabolites', '#e08a1e') +
    kpi(q.noFormula.length, 'Metabolites w/o formula') +
    kpi(q.orphanGenes.length, 'Orphan genes');
  const tbl = (title, rows, head, body) => `
    <div><h6>${title} (${rows.length})</h6>
      ${rows.length ? `<div class="fba-tablewrap"><table class="fba-flux"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`
                    : `<div class="fba-note" style="background:var(--ok-bg);border-color:var(--ok)">None found.</div>`}</div>`;
  $('qc-tables').innerHTML = `
    <div class="fba-two">
      ${tbl('Mass-imbalanced reactions', q.massImbalanced,
        '<th>Reaction</th><th>Name</th><th>Element imbalance</th>',
        q.massImbalanced.slice(0, 200).map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td>
          <td><code>${Object.entries(r.delta).map(([e, v]) => `${e}${v > 0 ? '+' : ''}${(+v.toFixed(2))}`).join(' ')}</code></td></tr>`).join(''))}
      ${tbl('Charge-imbalanced reactions', q.chargeImbalanced,
        '<th>Reaction</th><th>Name</th><th>Δ charge</th>',
        q.chargeImbalanced.slice(0, 200).map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td>
          <td class="num">${r.charge > 0 ? '+' : ''}${+r.charge.toFixed(2)}</td></tr>`).join(''))}
    </div>
    <div class="fba-two" style="margin-top:1rem">
      ${tbl('Blocked reactions (no flux possible on this medium)', bl.blocked.map(id => ({ id })),
        '<th>Reaction</th>', bl.blocked.slice(0, 300).map(r => `<tr><td><code>${esc(r.id)}</code></td></tr>`).join(''))}
      ${tbl('Dead-end metabolites', q.deadEnds,
        '<th>Metabolite</th><th>Name</th><th>Problem</th>',
        q.deadEnds.slice(0, 300).map(m => `<tr><td><code>${esc(m.id)}</code></td><td>${esc(m.name)}</td>
          <td>${m.onlyProduced ? 'only produced, never consumed' : 'only consumed, never produced'}</td></tr>`).join(''))}
    </div>`;
  $('qc-csv').onclick = () => {
    let c = 'category,id,detail\n';
    bl.blocked.forEach(r => c += `blocked,${r},\n`);
    q.massImbalanced.forEach(r => c += `mass_imbalanced,${r.id},"${Object.entries(r.delta).map(([e, v]) => e + v).join(' ')}"\n`);
    q.chargeImbalanced.forEach(r => c += `charge_imbalanced,${r.id},${r.charge}\n`);
    q.deadEnds.forEach(m => c += `dead_end,${m.id},${m.onlyProduced ? 'only_produced' : 'only_consumed'}\n`);
    q.orphanGenes.forEach(g => c += `orphan_gene,${g},\n`);
    saveCSV(c, `model_QC_${qcState.file.replace(/\.json.*/, '')}`);
  };
}
