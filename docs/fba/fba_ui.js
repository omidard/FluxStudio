// EcopanGEM — Flux Balance Analysis UI controller (ES module).
// Single & comparative FBA/pFBA with editable media, knockouts, a rich model
// picker, charts (Chart.js) and Escher flux maps. All client-side.
// Depends on page globals: gemBatchMap, BATCH_URL_BASE, JSZip, Chart, escher.
import { runFBA, runPFBA, runLMOMA, looplessSolution, exchangeReport, listExchanges, listReactions,
         bindMedium, exchangeIndex, resolveExchange } from './fba_engine.js';
import * as MediaDB from './media.js';
import { renderKOPlots } from './ko_plots.js';

const $ = (id) => document.getElementById(id);
const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls; if (html != null) e.innerHTML = html; return e; };
const fmt = (x, n = 4) => (x == null || isNaN(x)) ? '—' : Number(x).toFixed(n);
const pct = (x) => (100 * x).toFixed(0) + '%';

const S = {
  presets: null,
  mode: 'single',      // single | compare | ko
  charts: {},
  conditions: {},
  lastRun: null,
};

// Which pangenome a model belongs to. Drives the badge and the model filter.
const COLLECTION = {
  EcopanGEM:   { key: 'eco',   short: 'E. coli',           tag: 'ECO' },
  LactoPanGEM: { key: 'lacto', short: 'Lactobacillaceae',  tag: 'LACTO' },
};
const collOf = (m) => COLLECTION[m && m.dataset] || { key: 'other', short: '', tag: '' };

/* Lactobacillaceae are fastidious: they need amino acids, nucleotides and
   vitamins, and do not grow on a minimal medium at all. Measured across a
   sample of LactoPanGEM models, 0/10 grow on M9 + glucose while 5/10 grow on
   MRS or BHI. Handing a Lactobacillus model the E. coli default would show zero
   growth on arrival and read as a broken tool, so each collection gets a
   default medium that suits its organism. */
const MINIMAL = new Set(['M9_glucose_aerobic', 'M9_glucose_anaerobic']);
/* CDM (Ardalani et al., mSystems 2024, Table 2) is a chemically defined medium
   built for this family: 16 of 16 LactoPanGEM strains tested grow on it, against
   5 of 10 on MRS or BHI and 0 of 10 on M9. It is the right default here. */
const DEFAULT_MEDIUM = { EcopanGEM: 'M9_glucose_aerobic', LactoPanGEM: 'CDM' };
const organismOf = (m) => m.gtdb_species || m.organism || m.genome_name || m.strain || '';

// ── Media presets ─────────────────────────────────────────────────────────────
async function loadPresets() {
  if (!S.presets) S.presets = await (await fetch('fba/media_presets.json')).json();
  return S.presets;
}

/* ── Media binding ────────────────────────────────────────────────────────────
   A medium is a list of components keyed by BiGG exchange. It only becomes real
   once bound to a specific model, because the solver closes every exchange the
   medium does not name. EcopanGEM and LactoPanGEM were built a BiGG generation
   apart (EX_glc__D_e vs EX_glc_D_e), so binding has to resolve across both
   spellings. Nothing is applied silently: the coverage line below the picker
   states exactly how many compounds landed and which did not. */
function rebindMedia(cond) {
  if (!cond.model || !cond.mediaSpec) return;
  const r = bindMedium(cond.model, cond.mediaSpec.components, {
    openMinerals: cond.openMinerals,
  });
  cond.mediaBounds = r.bounds;
  cond.mediaReport = r;
  cond.mediaEdited = false;
  renderMediaReport(cond);
}

function specFromPreset(key) {
  const p = S.presets[key];
  return { kind: 'preset', id: key, label: p.label, desc: p.desc, source: p.source,
           components: MediaDB.presetToComponents(p) };
}

function specFromDB(m) {
  return { kind: 'db', id: m.id, label: m.name, source: (m.provenance && m.provenance.source_type) || m.namespace || 'Media DB',
           desc: m.description || '', citation: (m.provenance && m.provenance.citation) || '',
           components: m.components };
}

// ── Model loading (batch zip, cached) ────────────────────────────────────────
const modelCache = new Map();
async function loadModel(gemFile) {
  if (modelCache.has(gemFile)) return modelCache.get(gemFile);
  const batchNum = window.gemBatchMap && window.gemBatchMap[gemFile];
  if (!batchNum) throw new Error(`Unknown model "${gemFile}".`);
  const url = window.BATCH_URL_BASE + String(batchNum).padStart(2, '0') + '.zip';
  const resp = await fetch(url);
  if (!resp.ok) throw new Error('Failed to download model batch.');
  const zip = await window.JSZip.loadAsync(await resp.blob());
  const entry = zip.file(gemFile);
  if (!entry) throw new Error('Model not found inside batch.');
  const model = JSON.parse(await entry.async('text'));
  modelCache.set(gemFile, model);
  return model;
}
const metaByFile = new Map();
function indexMeta() {
  if (metaByFile.size || !window.gemMetadata) return;
  window.gemMetadata.forEach(m => metaByFile.set(m.gem_file, m));
}

// ── Condition card ────────────────────────────────────────────────────────────
function newCondition(slot) {
  return { slot, modelFile: null, model: null, meta: null, exchanges: null, reactions: null,
           mediaSpec: null, mediaBounds: {}, mediaReport: null, mediaEdited: false,
           openMinerals: false, knockouts: new Set(), card: null };
}

function buildConditionCard(slot) {
  const cond = newCondition(slot);
  cond.mediaSpec = specFromPreset('M9_glucose_aerobic');
  S.conditions[slot] = cond;
  const card = el('div', 'fba-cond'); card.dataset.slot = slot;
  card.innerHTML = `
    <div class="fba-cond-head"><span class="fba-cond-badge">${slot.toUpperCase()}</span>
      <span class="title">Condition ${slot.toUpperCase()}</span></div>

    <div class="fba-field">
      <label>Strain model (GEM)</label>
      <div class="coll-filter" role="group" aria-label="Filter models by pangenome">
        <button type="button" class="coll-chip active" data-coll="">All 4,659</button>
        <button type="button" class="coll-chip" data-coll="EcopanGEM"><i>E. coli</i> 2,313</button>
        <button type="button" class="coll-chip" data-coll="LactoPanGEM">Lactobacillaceae 2,346</button>
      </div>
      <div class="fba-combo">
        <input class="form-control form-control-sm fba-combo-input" placeholder="Search species, strain or accession…" autocomplete="off">
        <div class="fba-combo-menu"></div>
      </div>
      <div class="fba-modelcard" style="display:none"></div>
    </div>

    <div class="fba-field">
      <label>Growth medium
        <button type="button" class="fba-linkbtn me-toggle" disabled>edit compounds</button>
        <button type="button" class="fba-linkbtn me-reset" disabled>reset</button></label>
      <div class="media-quick"></div>
      <button type="button" class="media-browse">
        <svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" aria-hidden="true"><circle cx="11" cy="11" r="7"/><path d="M20 20l-4.3-4.3"/></svg>
        Browse <span class="mc">all</span> media
      </button>
      <div class="media-report"></div>
      <label class="media-min"><input type="checkbox" class="me-minerals"> open essential inorganic ions and water</label>
      <div class="fba-media-editor" style="display:none"></div>
    </div>

    <div class="fba-field">
      <label>Reaction knockouts <span class="fba-hint-inline">(optional — set flux to 0)</span></label>
      <div class="fba-ko">
        <input class="form-control form-control-sm fba-ko-input" placeholder="Search reaction ID or name to knock out…" autocomplete="off" disabled>
        <div class="fba-ko-menu"></div>
        <div class="fba-ko-chips"></div>
      </div>
    </div>`;
  cond.card = card;

  // quick media chips (the curated presets)
  const quick = card.querySelector('.media-quick');
  quick.innerHTML = Object.entries(S.presets).map(([k, v], i) =>
    `<button type="button" class="media-chip${i === 0 ? ' active' : ''}" data-key="${esc(k)}">${esc(v.label)}</button>`).join('');
  quick.addEventListener('click', (e) => {
    const b = e.target.closest('.media-chip'); if (!b) return;
    quick.querySelectorAll('.media-chip').forEach(x => x.classList.toggle('active', x === b));
    cond.mediaSpec = specFromPreset(b.dataset.key);
    cond.mediaTouched = true;
    rebindMedia(cond);
    if (card.querySelector('.fba-media-editor').style.display !== 'none') renderMediaEditor(cond);
  });

  MediaDB.fillCount(card.querySelector('.media-browse .mc'), n => `all ${n.toLocaleString()}`);
  card.querySelector('.media-browse').addEventListener('click', () => openMediaBrowser(cond));
  card.querySelector('.me-minerals').addEventListener('change', (e) => {
    cond.openMinerals = e.target.checked;
    rebindMedia(cond);
  });

  wireCombo(cond);
  wireMediaEditor(cond);
  wireKO(cond);
  renderMediaReport(cond);
  return card;
}

/* The coverage line. A medium that binds 32 of 56 compounds is not "mostly
   fine": the 24 that did not bind were silently removed from the medium. Say so. */
function renderMediaReport(cond) {
  const box = cond.card.querySelector('.media-report');
  const spec = cond.mediaSpec;
  if (!spec) { box.innerHTML = ''; return; }
  if (!cond.model) {
    box.innerHTML = `<div class="mr-line"><b>${esc(spec.label)}</b> · ${spec.components.length} compounds
      <span class="mr-hint">pick a model to bind it</span></div>`;
    return;
  }
  const r = cond.mediaReport || { mapped: [], missing: [], added: [], coverage: 0 };
  const n = spec.components.length;
  const cls = r.coverage >= 0.9 ? 'ok' : r.coverage >= 0.6 ? 'warn' : 'bad';
  const renamed = r.mapped.filter(m => m.renamed).length;
  box.innerHTML = `
    <div class="mr-line">
      <b>${esc(spec.label)}</b>
      <span class="mr-cov ${cls}">${r.mapped.length}/${n} bound · ${pct(r.coverage)}</span>
      ${cond.mediaEdited ? '<span class="mr-hint">edited</span>' : ''}
    </div>
    <div class="mr-sub">
      ${renamed ? `<span class="mr-fix" title="Resolved across BiGG naming generations, e.g. EX_glc__D_e to EX_glc_D_e">${renamed} id${renamed > 1 ? 's' : ''} resolved</span>` : ''}
      ${r.added.length ? `<span class="mr-hint">+${r.added.length} minerals opened</span>` : ''}
      ${r.missing.length ? `<button type="button" class="mr-missing">${r.missing.length} not in this model</button>` : ''}
      ${spec.source ? `<span class="mr-src">${esc(spec.source)}</span>` : ''}
    </div>
    <div class="mr-list" hidden>${r.missing.slice(0, 60).map(m =>
      `<code>${esc(m.exchange)}</code> <span>${esc(m.name || '')}</span>`).join('')}</div>`;
  const btn = box.querySelector('.mr-missing');
  if (btn) btn.addEventListener('click', () => {
    const l = box.querySelector('.mr-list'); l.hidden = !l.hidden;
  });
}

// ── Model combobox ────────────────────────────────────────────────────────────
function wireCombo(cond) {
  const input = cond.card.querySelector('.fba-combo-input');
  const menu = cond.card.querySelector('.fba-combo-menu');
  const chips = cond.card.querySelector('.coll-filter');
  let items = [], active = -1, collFilter = '';

  chips.addEventListener('click', (e) => {
    const b = e.target.closest('.coll-chip'); if (!b) return;
    chips.querySelectorAll('.coll-chip').forEach(x => x.classList.toggle('active', x === b));
    collFilter = b.dataset.coll || '';
    active = -1; render(); input.focus();
  });

  const render = () => {
    const q = input.value.trim().toLowerCase();
    let meta = window.gemMetadata || [];
    if (collFilter) meta = meta.filter(m => m.dataset === collFilter);
    items = (q ? meta.filter(m =>
        (m.gem_file || '').toLowerCase().includes(q) ||
        (m.gtdb_species || '').toLowerCase().includes(q) ||
        (m.organism || '').toLowerCase().includes(q) ||
        (m.genome_name || '').toLowerCase().includes(q) ||
        (m.strain || '').toLowerCase().includes(q))
      : meta).slice(0, 40);
    if (!items.length) { menu.innerHTML = `<div class="fba-combo-empty">No models match.</div>`; menu.classList.add('show'); return; }
    menu.innerHTML = items.map((m, i) => {
      const c = collOf(m);
      const facts = [
        m.strain ? esc(m.strain) : '',
        m.phylogroup ? 'Phylogroup ' + esc(m.phylogroup) : '',
        m.MLST ? 'ST' + esc(m.MLST) : '',
        m.isolation_source ? esc(m.isolation_source) : '',
        m.country ? esc(m.country) : '',
      ].filter(Boolean).slice(0, 3).join(' · ');
      const size = [m.n_reactions ? `${m.n_reactions} rxns` : '', m.n_genes ? `${m.n_genes} genes` : '']
        .filter(Boolean).join(' · ');
      return `
      <div class="fba-combo-item mdl${i === active ? ' active' : ''}" data-i="${i}">
        <span class="mdl-badge ${c.key}">${c.tag}</span>
        <span class="mdl-body">
          <span class="nm"><i>${esc(organismOf(m))}</i></span>
          <span class="meta">${facts}</span>
          <span class="id">${esc(m.assembly_accession || m.gem_file)}${size ? ' · ' + size : ''}</span>
        </span>
      </div>`;
    }).join('');
    menu.classList.add('show');
  };
  input.addEventListener('focus', render);
  input.addEventListener('input', () => { active = -1; render(); });
  input.addEventListener('keydown', (e) => {
    if (!menu.classList.contains('show')) return;
    if (e.key === 'ArrowDown') { active = Math.min(active + 1, items.length - 1); render(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { active = Math.max(active - 1, 0); render(); e.preventDefault(); }
    else if (e.key === 'Enter') { if (items[active] || items[0]) { pick((items[active] || items[0]).gem_file); e.preventDefault(); } }
    else if (e.key === 'Escape') menu.classList.remove('show');
  });
  menu.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-combo-item'); if (!it) return;
    e.preventDefault(); pick(items[+it.dataset.i].gem_file);
  });
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));

  const pick = async (gemFile) => {
    menu.classList.remove('show');
    input.value = '';
    await selectModel(cond, gemFile);
  };
}

async function selectModel(cond, gemFile) {
  indexMeta();
  cond.modelFile = gemFile;
  cond.meta = metaByFile.get(gemFile) || {};
  const mc = cond.card.querySelector('.fba-modelcard');
  mc.style.display = 'block';
  mc.innerHTML = `<span class="fba-hint-inline">Loading model…</span>`;
  try {
    const model = await loadModel(gemFile);
    cond.model = model;
    cond.exchanges = listExchanges(model);
    cond.reactions = listReactions(model);
    const m = cond.meta, c = collOf(m);
    const obj = (model.reactions.find(r => r.objective_coefficient) || {}).id || '—';
    mc.innerHTML = `
      <div class="mc-head">
        <span class="mdl-badge ${c.key}">${c.tag}</span>
        <span class="mc-name"><i>${esc(organismOf(m))}</i></span>
      </div>
      <div class="mc-id">${esc(m.assembly_accession || gemFile)}</div>
      <div class="mc-tags">
        <span class="fba-tag">${model.reactions.length} rxns</span>
        <span class="fba-tag">${model.metabolites.length} mets</span>
        <span class="fba-tag">${model.genes.length} genes</span>
        <span class="fba-tag">${cond.exchanges.length} exchanges</span>
        <span class="fba-tag" title="Objective (biomass) reaction">obj ${esc(obj)}</span>
        ${m.strain ? `<span class="fba-tag">${esc(m.strain)}</span>` : ''}
        ${m.phylogroup ? `<span class="fba-tag">Phylogroup ${esc(m.phylogroup)}</span>` : ''}
        ${m.MLST ? `<span class="fba-tag">ST ${esc(m.MLST)}</span>` : ''}
        ${m.isolation_source ? `<span class="fba-tag">${esc(m.isolation_source)}</span>` : ''}
        ${m.country ? `<span class="fba-tag">${esc(m.country)}</span>` : ''}
      </div>`;
    // enable editors
    cond.card.querySelector('.me-toggle').disabled = false;
    cond.card.querySelector('.me-reset').disabled = false;
    cond.card.querySelector('.fba-ko-input').disabled = false;

    // Give the organism a medium it can actually grow on, unless the user chose one.
    const want = DEFAULT_MEDIUM[m.dataset];
    let swapped = null;
    if (want && !cond.mediaTouched && cond.mediaSpec.kind === 'preset'
        && MINIMAL.has(cond.mediaSpec.id) && want !== cond.mediaSpec.id) {
      cond.mediaSpec = specFromPreset(want);
      cond.card.querySelectorAll('.media-chip').forEach(x =>
        x.classList.toggle('active', x.dataset.key === want));
      swapped = S.presets[want].label;
    }
    rebindMedia(cond);            // the medium is only real once bound to THIS model
    if (swapped) {
      cond.card.querySelector('.media-report').insertAdjacentHTML('beforeend',
        `<div class="mr-swap">Switched to <b>${esc(swapped)}</b>: Lactobacillaceae do not grow on minimal media.</div>`);
    }
  } catch (e) {
    mc.innerHTML = `<span style="color:var(--bad)">${esc(e.message)}</span>`;
  }
}

/* ── Media browser: search the whole Media DB ─────────────────────────────── */
async function openMediaBrowser(cond) {
  const ov = el('div', 'dlg-overlay');
  ov.innerHTML = `
    <div class="dlg" role="dialog" aria-modal="true" aria-label="Browse growth media">
      <div class="dlg-head">
        <div>
          <h5>Growth media</h5>
          <p class="dlg-sub"><b class="md-count">Curated</b> media, keyed to BiGG exchanges. DSMZ MediaDive, HMDB biospecimens, USDA foods, and media from GEM papers.</p>
        </div>
        <button class="dlg-x" aria-label="Close">&times;</button>
      </div>
      <div class="dlg-tools">
        <input class="form-control form-control-sm md-q" placeholder="Search media, e.g. CDM, MRS, LB, BHI, blood…" autocomplete="off">
        <select class="form-select form-select-sm md-cat">
          <option value="">All categories</option>
          <option value="laboratory">Laboratory</option>
          <option value="food">Food</option>
          <option value="biospecimen">Biospecimen</option>
        </select>
      </div>
      <div class="md-list"><div class="md-empty">Loading catalog…</div></div>
    </div>`;
  document.body.appendChild(ov);
  MediaDB.fillCategoryCounts(ov.querySelector('.md-cat'));
  const close = () => ov.remove();
  ov.querySelector('.dlg-x').addEventListener('click', close);
  ov.addEventListener('click', (e) => { if (e.target === ov) close(); });
  document.addEventListener('keydown', function esc2(e) {
    if (e.key === 'Escape') { close(); document.removeEventListener('keydown', esc2); }
  });

  const q = ov.querySelector('.md-q'), cat = ov.querySelector('.md-cat'), list = ov.querySelector('.md-list');
  let cat_ = null;
  try { cat_ = await MediaDB.catalog(); }
  catch (e) { list.innerHTML = `<div class="md-empty">Media DB unavailable: ${esc(e.message)}</div>`; return; }
  ov.querySelector('.md-count').textContent = cat_.length.toLocaleString();   // the catalog itself, not a guess

  const idx = cond.model ? exchangeIndex(cond.model) : null;
  const draw = () => {
    const res = MediaDB.search(cat_, q.value, { category: cat.value }, 80);
    if (!res.length) { list.innerHTML = `<div class="md-empty">No media match.</div>`; return; }
    list.innerHTML = res.map(m => `
      <button type="button" class="md-item" data-id="${esc(m.id)}">
        <span class="md-name">${esc(m.name)}</span>
        <span class="md-meta">
          <span class="md-tag ${esc(m.category)}">${esc(m.category)}</span>
          <span>${m.n_components} compounds</span>
          ${m.aerobic === true ? '<span>aerobic</span>' : m.aerobic === false ? '<span>anaerobic</span>' : ''}
          ${m.source_db ? `<span class="md-src">${esc(m.source_db)}</span>` : ''}
        </span>
      </button>`).join('');
  };
  q.addEventListener('input', draw);
  cat.addEventListener('change', draw);
  draw();
  q.focus();

  list.addEventListener('click', async (e) => {
    const b = e.target.closest('.md-item'); if (!b) return;
    b.classList.add('busy');
    try {
      const m = await MediaDB.medium(b.dataset.id);
      cond.mediaSpec = specFromDB(m);
      cond.mediaTouched = true;
      cond.card.querySelectorAll('.media-chip').forEach(x => x.classList.remove('active'));
      rebindMedia(cond);
      if (cond.card.querySelector('.fba-media-editor').style.display !== 'none') renderMediaEditor(cond);
      close();
    } catch (err) {
      b.classList.remove('busy');
      list.insertAdjacentHTML('afterbegin', `<div class="md-empty">Could not load: ${esc(err.message)}</div>`);
    }
  });
}

// ── Media editor ──────────────────────────────────────────────────────────────
function wireMediaEditor(cond) {
  const toggle = cond.card.querySelector('.me-toggle');
  const reset = cond.card.querySelector('.me-reset');
  const box = cond.card.querySelector('.fba-media-editor');
  toggle.addEventListener('click', () => {
    if (box.style.display === 'none') { box.style.display = 'block'; renderMediaEditor(cond); toggle.textContent = 'hide compounds'; }
    else { box.style.display = 'none'; toggle.textContent = 'edit compounds'; }
  });
  reset.addEventListener('click', () => {
    rebindMedia(cond);
    if (box.style.display !== 'none') renderMediaEditor(cond);
  });
}

function renderMediaEditor(cond) {
  const box = cond.card.querySelector('.fba-media-editor');
  const nameById = {}; (cond.exchanges || []).forEach(e => nameById[e.id] = e.name);
  const rows = Object.entries(cond.mediaBounds)
    .sort((a, b) => a[0].localeCompare(b[0]))
    .map(([id, lb]) => `
      <tr data-id="${id}">
        <td class="me-cmpd">${esc((nameById[id] || '').replace(/ exchange$/i, ''))}</td>
        <td class="me-id">${esc(id)}</td>
        <td><input type="number" step="0.1" class="me-rate" value="${lb}"></td>
        <td><span class="fba-me-rm" title="remove">✕</span></td>
      </tr>`).join('');
  box.innerHTML = `
    <div class="me-head"><b>Uptake bounds (mmol·gDW⁻¹·h⁻¹) — negative = uptake</b>
      <span class="fba-hint-inline">${Object.keys(cond.mediaBounds).length} compounds</span></div>
    <div class="fba-me-list"><table class="fba-me"><tbody>${rows}</tbody></table></div>
    <div class="fba-me-add">
      <input class="form-control form-control-sm me-add-input" placeholder="+ add compound (search exchange)…" autocomplete="off">
      <div class="fba-combo-menu me-add-menu"></div>
    </div>`;
  // rate edits
  box.querySelectorAll('.me-rate').forEach(inp => inp.addEventListener('change', () => {
    const id = inp.closest('tr').dataset.id;
    const v = parseFloat(inp.value);
    if (!isNaN(v)) { cond.mediaBounds[id] = v; cond.mediaEdited = true; renderMediaReport(cond); }
  }));
  // remove
  box.querySelectorAll('.fba-me-rm').forEach(x => x.addEventListener('click', () => {
    const id = x.closest('tr').dataset.id;
    delete cond.mediaBounds[id]; cond.mediaEdited = true; renderMediaEditor(cond); renderMediaReport(cond);
  }));
  // add-compound combobox
  const ai = box.querySelector('.me-add-input'), am = box.querySelector('.me-add-menu');
  let addItems = [];
  const renderAdd = () => {
    const q = ai.value.trim().toLowerCase();
    const avail = (cond.exchanges || []).filter(e => !(e.id in cond.mediaBounds));
    addItems = (q ? avail.filter(e => e.id.toLowerCase().includes(q) || (e.name || '').toLowerCase().includes(q)) : avail).slice(0, 30);
    am.innerHTML = addItems.length ? addItems.map((e, i) =>
      `<div class="fba-combo-item" data-i="${i}"><div class="nm">${esc((e.name || '').replace(/ exchange$/i, ''))}</div><div class="id">${esc(e.id)}</div></div>`).join('')
      : `<div class="fba-combo-empty">No more exchanges.</div>`;
    am.classList.add('show');
  };
  ai.addEventListener('focus', renderAdd); ai.addEventListener('input', renderAdd);
  ai.addEventListener('blur', () => setTimeout(() => am.classList.remove('show'), 150));
  am.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-combo-item'); if (!it) return; e.preventDefault();
    const ex = addItems[+it.dataset.i];
    cond.mediaBounds[ex.id] = -0.5; cond.mediaEdited = true; renderMediaEditor(cond); renderMediaReport(cond);
  });
}

// ── Knockout editor ───────────────────────────────────────────────────────────
function wireKO(cond) {
  const input = cond.card.querySelector('.fba-ko-input');
  const menu = cond.card.querySelector('.fba-ko-menu');
  let items = [];
  const render = () => {
    const q = input.value.trim().toLowerCase();
    if (!q || !cond.reactions) { menu.classList.remove('show'); return; }
    items = cond.reactions.filter(r => !cond.knockouts.has(r.id) &&
      (r.id.toLowerCase().includes(q) || (r.name || '').toLowerCase().includes(q))).slice(0, 25);
    menu.innerHTML = items.length ? items.map((r, i) =>
      `<div class="fba-ko-item" data-i="${i}"><span class="id">${esc(r.id)}</span> <span class="nm">${esc(r.name || '')}</span></div>`).join('')
      : `<div class="fba-combo-empty">No match.</div>`;
    menu.classList.add('show');
  };
  input.addEventListener('input', render); input.addEventListener('focus', render);
  input.addEventListener('blur', () => setTimeout(() => menu.classList.remove('show'), 150));
  menu.addEventListener('mousedown', (e) => {
    const it = e.target.closest('.fba-ko-item'); if (!it) return; e.preventDefault();
    cond.knockouts.add(items[+it.dataset.i].id); input.value = ''; menu.classList.remove('show'); renderChips(cond);
  });
}
function renderChips(cond) {
  const wrap = cond.card.querySelector('.fba-ko-chips');
  wrap.innerHTML = [...cond.knockouts].map(id =>
    `<span class="fba-ko-chip">${esc(id)} <span class="x" data-id="${esc(id)}">✕</span></span>`).join('');
  wrap.querySelectorAll('.x').forEach(x => x.addEventListener('click', () => { cond.knockouts.delete(x.dataset.id); renderChips(cond); }));
}

// ── Run ───────────────────────────────────────────────────────────────────────
function method() { return document.querySelector('input[name="fba-method"]:checked').value; }

async function runAnalysis(cond, meth) {
  return runWith(cond, meth, [...cond.knockouts]);
}

async function runWith(cond, meth, knockouts, ref) {
  const opts = { knockouts };
  const fba = await runFBA(cond.model, cond.mediaBounds, opts);
  let result = fba;
  const lethal = fba.optimal && fba.growth <= 1e-9;

  /* A dead cell's flux vector is not whatever the solver hands back. With the
     objective at zero, EVERY feasible vector is optimal, so FBA returns an arbitrary
     one: typically loops spinning at the +/-1000 bounds. Plotted, that reads as a
     corpse re-routing its metabolism, which is nonsense. So when the knockout is
     lethal we always minimise total flux, whatever method was asked for, and report
     the parsimonious state: the cell stops. */
  if (lethal) {
    const p = await runPFBA(cond.model, cond.mediaBounds, fba, opts);
    if (p.optimal) result = { ...p, growth: 0, lethal: true };
  } else if (meth === 'pfba' && fba.optimal) {
    result = await runPFBA(cond.model, cond.mediaBounds, fba, opts);
  } else if (meth === 'moma' && ref && fba.optimal) {
    // MOMA asks a different question: not "what maximises growth after the
    // knockout" but "what is the smallest change from how the cell was running".
    result = await runLMOMA(cond.model, cond.mediaBounds, ref, opts);
  }
  // CycleFreeFlux: same growth, same exchanges, no thermodynamically impossible loops
  if (looplessOn() && result.optimal && result.growth > 1e-9) {
    const ll = await looplessSolution(cond.model, cond.mediaBounds, result.fluxes, opts);
    if (ll.optimal) { result = { ...result, fluxes: ll.fluxes, loopless: true, loopsRemoved: ll.removed.length }; }
  }
  return { cond, fba, result, meth, knockouts };
}
function looplessOn() { const c = $('fba-loopless'); return !!(c && c.checked); }

/* ── Knockout study: the same model and medium, solved with and without the
   knockouts, so the effect is read as a difference rather than an absolute. ── */
function koVerdict(ratio) {
  if (ratio < 0.01) return { label: 'ESSENTIAL', color: 'var(--bad)' };
  if (ratio < 0.50) return { label: 'SEVERE', color: 'var(--warn)' };
  if (ratio < 0.95) return { label: 'MILD', color: 'var(--warn)' };
  return { label: 'NO EFFECT', color: 'var(--ok)' };
}

function renderKO(RW, RK) {
  const cond = RW.cond;
  const fW = RW.result.fluxes || {}, fK = RK.result.fluxes || {};
  const gW = RW.fba.optimal ? (RW.result.growth || 0) : 0;
  const gK = RK.fba.optimal ? (RK.result.growth || 0) : 0;
  const ratio = gW > 1e-9 ? gK / gW : 0;
  const v = koVerdict(ratio);
  const kos = [...cond.knockouts];

  // every reaction whose flux moved, largest change first. The knocked-out
  // reactions themselves are excluded from the re-routing count: they did not
  // re-route, they were removed, and counting the intervention as its own effect
  // inflates the number.
  const koSet = new Set(kos);
  const deltas = [];
  for (const r of cond.model.reactions) {
    const a = fW[r.id] || 0, b = fK[r.id] || 0, d = b - a;
    if (Math.abs(d) > 1e-6 && !koSet.has(r.id)) deltas.push({ id: r.id, name: r.name || '', wt: a, ko: b, d });
  }
  deltas.sort((x, y) => Math.abs(y.d) - Math.abs(x.d));

  const repW = exchangeReport(cond.model, fW), repK = exchangeReport(cond.model, fK);
  const setOf = (arr) => new Set(arr.map(x => x.id));
  const upW = setOf(repW.uptake), upK = setOf(repK.uptake);
  const seW = setOf(repW.secretion), seK = setOf(repK.secretion);
  const gained = [...seK].filter(x => !seW.has(x));
  const lost = [...seW].filter(x => !seK.has(x));

  const box = $('fba-results'); box.style.display = 'block';
  box.innerHTML = `
    <hr>
    <div class="ko-head">
      <code>${esc(cond.modelFile)}</code> on <strong>${esc(cond.mediaSpec.label)}</strong>
      <span class="fba-badge">${RW.meth.toUpperCase()}</span>
      <span class="ko-chipline">${kos.map(k => `<span class="fba-ko-chip">${esc(k)}</span>`).join('')}</span>
    </div>
    <div class="fba-kpis ko-kpis">
      ${kpi(fmt(gW), 'Growth, wild type (h⁻¹)', 'var(--primary)')}
      ${kpi(gK <= 1e-9 ? '0' : fmt(gK), 'Growth, knockout (h⁻¹)', gK <= 1e-9 ? 'var(--bad)' : 'var(--ink)')}
      ${kpi((gK - gW >= 0 ? '+' : '') + fmt(gK - gW), 'Δ growth', (gK - gW) < -1e-9 ? 'var(--bad)' : 'var(--ok)')}
      ${kpi(pct(ratio), 'of wild type', v.color)}
      ${kpi(v.label, 'Verdict', v.color)}
      ${kpi(deltas.length, gK <= 1e-9 ? 'Reactions shut down' : 'Reactions re-routed')}
    </div>

    <div class="ko-verdict" id="ko-verdict"></div>

    <div class="ko-viz">
      <div class="ko-viz-row">
        <figure class="ko-fig">
          <h6>How much of it did the cell actually need?</h6>
          <div id="ko-plot-titr" class="ko-plot"></div>
          <figcaption id="ko-cap-titr"></figcaption>
        </figure>
        <figure class="ko-fig">
          <h6>Where in metabolism did it land?</h6>
          <div id="ko-plot-subsys" class="ko-plot"></div>
          <figcaption id="ko-cap-subsys"></figcaption>
        </figure>
      </div>

      <figure class="ko-fig wide">
        <h6>The whole flux distribution, wild type against knockout</h6>
        <div id="ko-plot-scatter" class="ko-plot tall"></div>
        <figcaption id="ko-cap-scatter"></figcaption>
      </figure>

      <div class="ko-viz-row">
        <figure class="ko-fig">
          <h6>Which reactions moved most</h6>
          <div id="ko-plot-tornado" class="ko-plot tall"></div>
          <figcaption id="ko-cap-tornado"></figcaption>
        </figure>
        <figure class="ko-fig">
          <h6>What the cell eats and excretes</h6>
          <div id="ko-plot-exch" class="ko-plot tall"></div>
          <figcaption id="ko-cap-exch"></figcaption>
        </figure>
      </div>
    </div>

    <div class="ko-maps">
      <div class="ko-map">
        <h6>Wild type<button class="viz-dl" data-plot="ko-map-wt" data-type="escher">⬇ SVG</button></h6>
        <div id="ko-map-wt" class="ko-map-box"></div>
      </div>
      <div class="ko-map">
        <h6>Knockout <span class="ko-dim">(${kos.length} reaction${kos.length > 1 ? 's' : ''} removed)</span><button class="viz-dl" data-plot="ko-map-ko" data-type="escher">⬇ SVG</button></h6>
        <div id="ko-map-ko" class="ko-map-box"></div>
      </div>
    </div>
    <div class="fba-map-note" id="ko-map-note"></div>

    <div class="fba-two" style="margin-top:1.2rem">
      <div>
        <h6>Largest flux changes (${deltas.length})</h6>
        <div class="fba-tablewrap">
          <table class="fba-flux">
            <thead><tr><th>Reaction</th><th>Wild type</th><th>Knockout</th><th>Δ</th></tr></thead>
            <tbody>${deltas.slice(0, 200).map(x => `
              <tr title="${esc(x.name)}">
                <td><code>${esc(x.id)}</code></td>
                <td class="num">${fmt(x.wt, 3)}</td>
                <td class="num">${fmt(x.ko, 3)}</td>
                <td class="num" style="color:${x.d < 0 ? 'var(--bad)' : 'var(--ok)'}">${(x.d >= 0 ? '+' : '') + fmt(x.d, 3)}</td>
              </tr>`).join('')}</tbody>
          </table>
        </div>
      </div>
      <div>
        <h6>What the cell does differently</h6>
        <div class="ko-shift">
          <div class="ko-shift-row"><span class="l">Nutrients taken up</span><span class="v">${repW.uptake.length} <span class="ko-arrow">→</span> ${repK.uptake.length}</span></div>
          <div class="ko-shift-row"><span class="l">Products secreted</span><span class="v">${repW.secretion.length} <span class="ko-arrow">→</span> ${repK.secretion.length}</span></div>
          ${gained.length ? `<div class="ko-shift-list"><b>New products</b>${gained.slice(0, 14).map(x => `<code class="up">${esc(bio({ id: x }))}</code>`).join('')}</div>` : ''}
          ${lost.length ? `<div class="ko-shift-list"><b>Products lost</b>${lost.slice(0, 14).map(x => `<code class="dn">${esc(bio({ id: x }))}</code>`).join('')}</div>` : ''}
          ${!gained.length && !lost.length ? `<div class="ko-shift-list"><span class="fba-hint-inline">Same set of secreted products; only the magnitudes moved.</span></div>` : ''}
        </div>
        <button class="btn btn-sm btn-outline-secondary mt-2" id="ko-csv">⬇ Download Δ fluxes (CSV)</button>
      </div>
    </div>`;

  renderEscher('ko-map-wt', null, fW, cond.model);
  renderEscher('ko-map-ko', 'ko-map-note', fK, cond.model);

  renderKOPlots({
    model: cond.model, media: cond.mediaBounds,
    fW, fK, gW, gK, kos, method: RW.meth,
  }).catch(e => { $('ko-verdict').innerHTML = `<span class="kv-note">Plots failed: ${esc(e.message)}</span>`; });

  $('ko-csv').addEventListener('click', () => {
    const rows = [['reaction', 'name', 'wild_type_flux', 'knockout_flux', 'delta']]
      .concat(deltas.map(x => [x.id, (x.name || '').replace(/,/g, ';'), x.wt, x.ko, x.d]));
    saveCSV(rows.map(r => r.join(',')).join('\n'),
      `KO_${cond.modelFile.replace(/\.json.*/, '')}_${kos.join('-')}`);
  });
}

async function run() {
  const meth = method();
  const btn = $('fba-run'); btn.disabled = true;
  $('fba-results').style.display = 'none';
  try {
    if (S.mode === 'single') {
      const a = S.conditions.a;
      if (!a.model) return setStatus('Choose a model for Condition A first.', 'err');
      setStatus(`Solving ${meth.toUpperCase()} in your browser…`, 'busy');
      const t0 = performance.now();
      const RA = await runAnalysis(a, meth);
      // growth across media (with current knockouts)
      const across = {};
      for (const k of Object.keys(S.presets)) {
        if (a.mediaSpec.kind === 'preset' && k === a.mediaSpec.id && !a.mediaEdited) { across[k] = RA.fba.growth; continue; }
        // bind each preset to THIS model: a Lactobacillaceae GEM needs its own
        // exchange spelling, or every medium reads as zero growth.
        const b = bindMedium(a.model, MediaDB.presetToComponents(S.presets[k]), { openMinerals: a.openMinerals });
        across[k] = (await runFBA(a.model, b.bounds, { knockouts: [...a.knockouts] })).growth;
      }
      RA.across = across;
      S.lastRun = { mode: 'single', RA };
      renderSingle(RA);
      setStatus(`Done — solved in ${(performance.now() - t0).toFixed(0)} ms on your machine.`, 'ok');
    } else if (S.mode === 'ko') {
      const a = S.conditions.a;
      if (!a.model) return setStatus('Choose a model first.', 'err');
      if (!a.knockouts.size) return setStatus('Add at least one reaction knockout to compare against.', 'err');
      setStatus(`Solving wild type and knockout (${meth.toUpperCase()})…`, 'busy');
      const t0 = performance.now();
      // MOMA measures distance FROM the wild type, so the wild type is solved by
      // pFBA first and handed to the knockout as its reference state.
      const wtMeth = meth === 'moma' ? 'pfba' : meth;
      const RW = await runWith(a, wtMeth, []);
      const RK = await runWith(a, meth, [...a.knockouts], RW.result.fluxes);
      S.lastRun = { mode: 'ko', RA: RW, RB: RK };
      renderKO(RW, RK);
      setStatus(`Done — both solved in ${(performance.now() - t0).toFixed(0)} ms.`, 'ok');
    } else {
      const a = S.conditions.a, b = S.conditions.b;
      if (!a.model || !b.model) return setStatus('Choose a model for both Condition A and B.', 'err');
      setStatus(`Solving both conditions (${meth.toUpperCase()})…`, 'busy');
      const t0 = performance.now();
      const RA = await runAnalysis(a, meth);
      const RB = await runAnalysis(b, meth);
      S.lastRun = { mode: 'compare', RA, RB };
      renderCompare(RA, RB);
      setStatus(`Done — both solved in ${(performance.now() - t0).toFixed(0)} ms.`, 'ok');
    }
  } catch (e) { setStatus('Error: ' + e.message, 'err'); console.error(e); }
  finally { btn.disabled = false; }
}

function setStatus(msg, kind) { const e = $('fba-status'); e.textContent = msg; e.className = 'fba-status ' + (kind || ''); }
function kpi(v, l, color) { return `<div class="fba-kpi"><div class="v"${color ? ` style="color:${color}"` : ''}>${v}</div><div class="l">${l}</div></div>`; }
function condLabel(cond) { return `${cond.mediaSpec ? cond.mediaSpec.label : '—'}${cond.mediaEdited ? ' *' : ''}${cond.knockouts.size ? ` · KO:${[...cond.knockouts].join(',')}` : ''}`; }

// ── Single results ────────────────────────────────────────────────────────────
function renderSingle(R) {
  const cond = R.cond, fluxes = R.result.fluxes, rep = exchangeReport(cond.model, fluxes);
  const growth = R.result.growth || 0, infeasible = !R.fba.optimal || growth <= 1e-9;
  const box = $('fba-results'); box.style.display = 'block';
  box.innerHTML = `
    <hr>
    <div style="font-size:0.95rem;margin-bottom:0.4rem"><code>${esc(cond.modelFile)}</code> on <strong>${esc(cond.mediaSpec.label)}</strong>
      <span class="fba-badge">${R.meth.toUpperCase()}</span> ${cond.knockouts.size ? `<span class="fba-badge" style="background:var(--bad)">${cond.knockouts.size} KO</span>` : ''}</div>
    <div class="fba-kpis">
      ${kpi(infeasible ? '0' : fmt(growth), 'Growth rate (h⁻¹)', infeasible ? '#c0392b' : '#1a7f4b')}
      ${kpi(infeasible ? 'NO GROWTH' : 'FEASIBLE', 'Status', infeasible ? '#c0392b' : '#1a7f4b')}
      ${kpi(rep.uptake.length, 'Nutrients taken up')}
      ${kpi(rep.secretion.length, 'Products secreted')}
      ${kpi(R.result.pfba ? fmt(R.result.totalFlux, 0) : '—', 'Total flux Σ|v| (pFBA)')}
    </div>`;
  if (infeasible) {
    const isLacto = (cond.meta || {}).dataset === 'LactoPanGEM';
    const minimal = cond.mediaSpec.kind === 'preset' && MINIMAL.has(cond.mediaSpec.id);
    const why = cond.knockouts.size
      ? 'The knockouts remove a reaction the model cannot route around on this medium, so this reaction set is essential here.'
      : (isLacto && minimal)
        ? 'Lactobacillaceae are fastidious: they need amino acids, nucleotides and vitamins, and do not grow on a minimal medium. Try <b>MRS</b> or <b>BHI</b>.'
        : isLacto
          ? 'This medium does not cover everything this strain needs. Not every LactoPanGEM model is gap-filled for every medium, so a rich medium can still leave a strain unable to grow. Try MRS or BHI, or open the essential inorganic ions.'
          : 'No biomass on this medium. Try a richer medium, or open the essential inorganic ions.';
    box.insertAdjacentHTML('beforeend',
      `<div class="fba-note"><b>No growth (biomass ≈ 0).</b> ${why}</div>`);
    return;
  }
  box.insertAdjacentHTML('beforeend', `
    <div class="fba-charts">
      <div class="fba-chart-card"><h6>Growth across media${cond.knockouts.size ? ' (with knockouts)' : ''}<button class="viz-dl" data-plot="ch-across" data-type="chartjs" style="float:right">⬇</button></h6><div class="fba-chart-box"><canvas id="ch-across"></canvas></div></div>
      <div class="fba-chart-card"><h6>Top secreted end products<button class="viz-dl" data-plot="ch-sec" data-type="chartjs" style="float:right">⬇</button></h6><div class="fba-chart-box"><canvas id="ch-sec"></canvas></div></div>
      <div class="fba-chart-card"><h6>Top nutrient uptakes<button class="viz-dl" data-plot="ch-up" data-type="chartjs" style="float:right">⬇</button></h6><div class="fba-chart-box"><canvas id="ch-up"></canvas></div></div>
    </div>
    <div class="fba-map-wrap"><h6 style="font-size:0.85rem;font-weight:700;color:var(--ink);margin-bottom:0.4rem">Metabolic flux map<button class="viz-dl" data-plot="ch-map" data-type="escher" style="float:right">⬇ SVG</button></h6>
      <div id="ch-map" style="width:100%;height:560px;border:1px solid var(--canvas-line);border-radius:8px;background:#fff;overflow:hidden"></div>
      <div class="fba-map-note" id="ch-map-note"></div></div>
    <div id="ch-tables" style="margin-top:1.2rem"></div>`);

  barChart('ch-across', Object.keys(R.across).map(k => S.presets[k].label),
    Object.values(R.across).map(v => Math.max(0, v)),
    Object.keys(R.across).map(k => (cond.mediaSpec.kind === 'preset' && k === cond.mediaSpec.id) ? '#1a7f4b' : '#9db8d6'), 'Growth (h⁻¹)', false);
  barChart('ch-sec', rep.secretion.slice(0, 10).map(bio), rep.secretion.slice(0, 10).map(x => x.flux), '#c0392b', 'Flux', true);
  barChart('ch-up', rep.uptake.slice(0, 10).map(bio), rep.uptake.slice(0, 10).map(x => Math.abs(x.flux)), '#2c6fbb', 'Flux', true);
  renderEscher('ch-map', 'ch-map-note', fluxes, cond.model);
  singleTables(cond, R, rep);
}

function singleTables(cond, R, rep) {
  const rows = (arr, s) => arr.map(x => `<tr><td><code>${esc(x.id)}</code></td><td>${esc(x.name || '')}</td><td class="num" style="color:${s > 0 ? 'var(--bad)' : 'var(--primary)'}">${fmt(x.flux)}</td></tr>`).join('');
  $('ch-tables').innerHTML = `
    <div class="fba-two">
      <div><h6>Uptake (${rep.uptake.length})</h6><div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead><tbody>${rows(rep.uptake, -1)}</tbody></table></div></div>
      <div><h6>Secretion / end products (${rep.secretion.length})</h6><div class="fba-tablewrap"><table class="fba-flux"><thead><tr><th>Exchange</th><th>Name</th><th>Flux</th></tr></thead><tbody>${rows(rep.secretion, 1)}</tbody></table></div></div>
    </div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="ch-csv">⬇ Download all fluxes (CSV)</button>`;
  $('ch-csv').addEventListener('click', () => downloadFluxCSV(cond.model, R.result.fluxes, `FBA_${cond.modelFile.replace(/\.json.*/, '')}_${cond.mediaSpec.id}_${R.meth}`));
}

// ── Compare results ───────────────────────────────────────────────────────────
function renderCompare(RA, RB) {
  const gA = RA.result.growth || 0, gB = RB.result.growth || 0;
  const fA = RA.result.fluxes, fB = RB.result.fluxes;
  const infA = !RA.fba.optimal || gA <= 1e-9, infB = !RB.fba.optimal || gB <= 1e-9;
  const delta = gB - gA, pct = gA > 1e-9 ? (delta / gA * 100) : null;
  const box = $('fba-results'); box.style.display = 'block';
  box.innerHTML = `
    <hr>
    <div class="fba-cmp-heads">
      <div><span class="fba-cond-badge">A</span> <code>${esc(RA.cond.modelFile)}</code> — ${esc(condLabel(RA.cond))}</div>
      <div><span class="fba-cond-badge" style="background:var(--bad)">B</span> <code>${esc(RB.cond.modelFile)}</code> — ${esc(condLabel(RB.cond))}</div>
    </div>
    <div class="fba-kpis">
      ${kpi(infA ? '0' : fmt(gA), 'Growth A (h⁻¹)', '#2c6fbb')}
      ${kpi(infB ? '0' : fmt(gB), 'Growth B (h⁻¹)', '#c0392b')}
      ${kpi((delta >= 0 ? '+' : '') + fmt(delta), 'Δ Growth (B−A)', delta >= 0 ? '#1a7f4b' : '#c0392b')}
      ${kpi(pct == null ? '—' : (pct >= 0 ? '+' : '') + pct.toFixed(1) + '%', 'Δ Growth %')}
    </div>
    <div class="fba-charts" style="grid-template-columns:1fr 1fr">
      <div class="fba-chart-card"><h6>Growth: A vs B</h6><div class="fba-chart-box"><canvas id="cc-growth"></canvas></div></div>
      <div class="fba-chart-card"><h6>Secreted end products: A vs B</h6><div class="fba-chart-box"><canvas id="cc-sec"></canvas></div></div>
    </div>
    <div class="fba-two" style="margin-top:0.4rem">
      <div class="fba-map-wrap"><h6 style="font-size:0.82rem;font-weight:700;color:var(--primary)">Flux map — Condition A</h6>
        <div id="cc-map-a" style="width:100%;height:460px;border:1px solid var(--canvas-line);border-radius:8px;background:#fff;overflow:hidden"></div></div>
      <div class="fba-map-wrap"><h6 style="font-size:0.82rem;font-weight:700;color:var(--bad)">Flux map — Condition B</h6>
        <div id="cc-map-b" style="width:100%;height:460px;border:1px solid var(--canvas-line);border-radius:8px;background:#fff;overflow:hidden"></div></div>
    </div>
    <div class="fba-map-note">Each map is coloured by its own condition's |flux|. Reactions outside the e_coli_core central-metabolism map are not shown.</div>
    <div id="cc-diff" style="margin-top:1.2rem"></div>`;

  barChart('cc-growth', ['Condition A', 'Condition B'], [Math.max(0, gA), Math.max(0, gB)], ['#2c6fbb', '#c0392b'], 'Growth (h⁻¹)', false);

  // secretion comparison over union of top secretions
  const repA = exchangeReport(RA.cond.model, fA), repB = exchangeReport(RB.cond.model, fB);
  const topIds = [...new Set([...repA.secretion.slice(0, 8), ...repB.secretion.slice(0, 8)].map(x => x.id))];
  const mapA = Object.fromEntries(repA.secretion.map(x => [x.id, x.flux]));
  const mapB = Object.fromEntries(repB.secretion.map(x => [x.id, x.flux]));
  groupedChart('cc-sec', topIds.map(id => id.replace(/^EX_/, '').replace(/_e$/, '')),
    topIds.map(id => mapA[id] || 0), topIds.map(id => mapB[id] || 0), 'Secretion flux');

  renderEscher('cc-map-a', null, fA, RA.cond.model);
  renderEscher('cc-map-b', null, fB, RB.cond.model);
  diffTable(RA, RB, fA, fB);
}

function diffTable(RA, RB, fA, fB) {
  const ids = new Set([...Object.keys(fA), ...Object.keys(fB)]);
  const nameById = {}; RA.cond.model.reactions.forEach(r => nameById[r.id] = r.name || ''); RB.cond.model.reactions.forEach(r => { if (!nameById[r.id]) nameById[r.id] = r.name || ''; });
  const rows = [...ids].map(id => ({ id, a: fA[id] || 0, b: fB[id] || 0, d: (fB[id] || 0) - (fA[id] || 0), name: nameById[id] || '' }))
    .filter(r => Math.abs(r.d) > 1e-6).sort((x, y) => Math.abs(y.d) - Math.abs(x.d));
  const top = rows.slice(0, 60);
  $('cc-diff').innerHTML = `
    <h6>Largest flux differences (B − A) <span class="fba-hint-inline">${rows.length} reactions differ; top ${top.length} shown</span></h6>
    <div class="fba-tablewrap" style="max-height:360px"><table class="fba-flux"><thead><tr><th>Reaction</th><th>Name</th><th>Flux A</th><th>Flux B</th><th>Δ (B−A)</th></tr></thead>
      <tbody>${top.map(r => `<tr><td><code>${esc(r.id)}</code></td><td>${esc(r.name)}</td><td class="num">${fmt(r.a, 3)}</td><td class="num">${fmt(r.b, 3)}</td><td class="num" style="color:${r.d >= 0 ? 'var(--ok)' : 'var(--bad)'}">${(r.d >= 0 ? '+' : '') + fmt(r.d, 3)}</td></tr>`).join('')}</tbody></table></div>
    <button class="btn btn-sm btn-outline-secondary mt-2" id="cc-csv">⬇ Download flux comparison (CSV)</button>`;
  $('cc-csv').addEventListener('click', () => {
    let csv = 'reaction_id,name,flux_A,flux_B,delta_B_minus_A\n';
    rows.forEach(r => { csv += `${r.id},"${(r.name || '').replace(/"/g, '""')}",${r.a},${r.b},${r.d}\n`; });
    saveCSV(csv, `FBA_compare_${RA.cond.modelFile.replace(/\.json.*/, '')}_vs_${RB.cond.modelFile.replace(/\.json.*/, '')}`);
  });
}

// ── Charts ────────────────────────────────────────────────────────────────────
function destroyChart(id) { if (S.charts[id]) { S.charts[id].destroy(); delete S.charts[id]; } }
function barChart(canvasId, labels, data, colors, axisTitle, horizontal) {
  destroyChart(canvasId);
  S.charts[canvasId] = new Chart($(canvasId), {
    type: 'bar', data: { labels, datasets: [{ data, backgroundColor: colors, borderRadius: 4 }] },
    options: {
      indexAxis: horizontal ? 'y' : 'x', responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: false }, tooltip: { callbacks: { label: c => ' ' + c.parsed[horizontal ? 'x' : 'y'].toFixed(3) } } },
      scales: { [horizontal ? 'x' : 'y']: { title: { display: true, text: axisTitle, font: { size: 11 } }, ticks: { font: { size: 10 } } }, [horizontal ? 'y' : 'x']: { ticks: { font: { size: 10 } } } },
    },
  });
}
function groupedChart(canvasId, labels, dataA, dataB, axisTitle) {
  destroyChart(canvasId);
  S.charts[canvasId] = new Chart($(canvasId), {
    type: 'bar',
    data: { labels, datasets: [{ label: 'A', data: dataA, backgroundColor: '#2c6fbb', borderRadius: 3 }, { label: 'B', data: dataB, backgroundColor: '#c0392b', borderRadius: 3 }] },
    options: { indexAxis: 'y', responsive: true, maintainAspectRatio: false,
      plugins: { legend: { display: true, labels: { font: { size: 10 }, boxWidth: 12 } } },
      scales: { x: { title: { display: true, text: axisTitle, font: { size: 11 } }, ticks: { font: { size: 10 } } }, y: { ticks: { font: { size: 10 } } } } },
  });
}
function bio(x) { return x.id.replace(/^EX_/, '').replace(/_e$/, ''); }

// ── Escher ────────────────────────────────────────────────────────────────────
let coreMap = null;
/* Re-key a model's flux vector onto the reaction ids the Escher map actually uses.
   Two mismatches to bridge: (1) the objective/biomass reaction id differs per model
   (e.g. BIOMASS_Ec_iML1515_core_75p37M) from the map's (BIOMASS_Ecoli_core_w_GAM),
   so biomass never coloured; (2) BiGG id generations differ (stereo "__", "(e)"→"_e").
   Returns fluxes augmented with the map's own ids so they render. */
function remapFluxesToMap(fluxes, model, map) {
  const mapData = Array.isArray(map) ? map[1] : map;
  const mapIds = Object.values(mapData.reactions || {}).map(r => r.bigg_id);
  const out = Object.assign({}, fluxes);
  // (1) biomass: the model's objective flux -> every biomass reaction on the map
  const obj = model && (model.reactions.find(r => r.objective_coefficient) || {}).id;
  if (obj && fluxes[obj] != null) for (const id of mapIds) if (/biomass/i.test(id)) out[id] = fluxes[obj];
  // (2) normalized fallback for stereo/compartment id-form differences
  const norm = s => String(s).toLowerCase().replace(/[\[(]([a-z])[\])]$/, '_$1').replace(/_+/g, '_');
  const byNorm = {}; for (const [k, v] of Object.entries(fluxes)) byNorm[norm(k)] = v;
  for (const id of mapIds) if (out[id] == null) { const v = byNorm[norm(id)]; if (v != null) out[id] = v; }
  return out;
}
// The bundled Escher maps the user can switch between (E. coli reference maps —
// the pan-genome models share most BiGG ids with these).
const ESCHER_MAPS = [
  { id: 'core', label: 'Core metabolism (e_coli_core)', file: 'fba/maps/core.json' },
  { id: 'central', label: 'Central metabolism (iJO1366)', file: 'fba/maps/central.json' },
  { id: 'nucleotide', label: 'Nucleotide & histidine (iJO1366)', file: 'fba/maps/nucleotide.json' },
  { id: 'fa_biosynth', label: 'Fatty-acid biosynthesis (iJO1366)', file: 'fba/maps/fa_biosynth.json' },
  { id: 'fa_betaox', label: 'Fatty-acid β-oxidation (iJO1366)', file: 'fba/maps/fa_betaox.json' },
];
let escherMapId = 'core';
const _mapCache = {};                 // id -> parsed map json
const _escherState = {};              // containerId -> {noteId, fluxes, model}
async function _loadMap(id) {
  const m = ESCHER_MAPS.find(x => x.id === id) || ESCHER_MAPS[0];
  if (!_mapCache[m.id]) _mapCache[m.id] = await (await fetch(m.file)).json();
  return _mapCache[m.id];
}
window.setEscherMap = async function (id) {                 // called by the switcher
  escherMapId = id;
  for (const [cid, st] of Object.entries(_escherState)) await renderEscher(cid, st.noteId, st.fluxes, st.model);
};
function _mapSwitcher() {
  return `<select class="escher-mapsel" onchange="setEscherMap(this.value)" title="Switch metabolic map">`
    + ESCHER_MAPS.map(m => `<option value="${m.id}"${m.id === escherMapId ? ' selected' : ''}>${esc(m.label)}</option>`).join('')
    + `</select>`;
}
async function renderEscher(containerId, noteId, fluxes, model) {
  const note = noteId ? $(noteId) : null;
  const cont = $(containerId); if (!cont) return;
  _escherState[containerId] = { noteId, fluxes, model };
  if (!window.escher) { if (note) note.textContent = 'Escher unavailable.'; return; }
  try {
    const map = await _loadMap(escherMapId);
    // one map switcher per container, inserted just above it
    let sel = document.getElementById(containerId + '-mapsel');
    if (!sel) { const d = document.createElement('div'); d.id = containerId + '-mapsel'; d.className = 'escher-mapbar';
      d.innerHTML = `<span class="escher-maplbl">Map</span> ${_mapSwitcher()}`; cont.parentNode.insertBefore(d, cont); }
    else { const s = sel.querySelector('select'); if (s) s.value = escherMapId; }
    const mapped = remapFluxesToMap(fluxes, model, map);
    cont.innerHTML = '';
    window.escher.Builder(map, null, null, window.escher.libs.d3_select('#' + containerId), {
      menu: 'zoom', scroll_behavior: 'zoom', fill_screen: false, never_ask_before_quit: true,
      reaction_data: mapped, reaction_styles: ['color', 'size', 'text'],
      reaction_scale: [
        { type: 'min', color: '#d0d0d0', size: 6 }, { type: 'value', value: 0, color: '#d0d0d0', size: 6 },
        { type: 'median', color: '#5b8ff9', size: 14 }, { type: 'max', color: '#c0392b', size: 28 }],
    });
    const label = (ESCHER_MAPS.find(x => x.id === escherMapId) || {}).label || '';
    if (note) note.textContent = `${label}. Colour & thickness = |flux|; biomass is mapped from the model's objective. Hover a reaction for its value; reactions absent from this map are not shown.`;
  } catch (e) { if (note) note.textContent = 'Map error: ' + e.message; console.error(e); }
}

// ── CSV ───────────────────────────────────────────────────────────────────────
function downloadFluxCSV(model, fluxes, base) {
  const nameById = {}; model.reactions.forEach(r => nameById[r.id] = r.name || '');
  let csv = 'reaction_id,name,flux\n';
  for (const [id, v] of Object.entries(fluxes)) csv += `${id},"${(nameById[id] || '').replace(/"/g, '""')}",${v}\n`;
  saveCSV(csv, base);
}
function saveCSV(csv, base) {
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' }));
  a.download = base + '.csv'; a.click();
}
function esc(s) { return String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c])); }

// ── Mode toggle & init ────────────────────────────────────────────────────────
function setMode(mode) {
  S.mode = mode;
  document.querySelectorAll('#fba-modetoggle button').forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
  const wrap = $('fba-conditions');
  wrap.classList.toggle('compare', mode === 'compare');
  S.conditions.b.card.style.display = mode === 'compare' ? '' : 'none';
  const runBtn = $('fba-run');
  if (runBtn) runBtn.textContent = mode === 'ko' ? '▶ Run knockout study' : '▶ Run analysis';
  const koHint = $('fba-ko-hint');
  if (koHint) koHint.style.display = mode === 'ko' ? '' : 'none';
  const momaOpt = $('fba-moma-opt');
  if (momaOpt) {
    momaOpt.style.display = mode === 'ko' ? '' : 'none';
    const r = momaOpt.querySelector('input');
    if (mode !== 'ko' && r.checked) { r.checked = false; document.querySelector('input[name="fba-method"][value="pfba"]').checked = true; }
  }
}

async function init() {
  await loadPresets();
  const wrap = $('fba-conditions');
  wrap.appendChild(buildConditionCard('a'));
  wrap.appendChild(buildConditionCard('b'));
  S.conditions.b.card.style.display = 'none';
  document.querySelectorAll('#fba-modetoggle button').forEach(b => b.addEventListener('click', () => setMode(b.dataset.mode)));
  $('fba-run').addEventListener('click', run);

  /* Handoff from a GEM detail modal, from the model table, or from a deep link out of
     any of the model browsers. opts: {medium: presetKey, ko: [reactionId]}. */
  window.fbaSetModel = async (slot, gemFile, opts = {}) => {
    slot = (slot === 'b') ? 'b' : 'a';
    if (slot === 'b') setMode('compare');
    const cond = S.conditions[slot];

    /* A medium named in the link is an explicit choice, so mark it as touched: the
       collection default would otherwise swap itself in when the model lands. */
    if (opts.medium && S.presets && S.presets[opts.medium]) {
      cond.mediaSpec = specFromPreset(opts.medium);
      cond.mediaTouched = true;
      cond.card.querySelectorAll('.media-chip').forEach(x =>
        x.classList.toggle('active', x.dataset.key === opts.medium));
    }
    const ko = (opts.ko || []).filter(Boolean);
    if (ko.length) setMode('ko');

    await selectModel(cond, gemFile);

    /* Knockouts go on only once the model is in, so a reaction this strain does not
       carry is dropped rather than left sitting in a chip that knocks out nothing. */
    if (ko.length && cond.model) {
      const have = new Set(cond.model.reactions.map(r => r.id));
      ko.filter(id => have.has(id)).forEach(id => cond.knockouts.add(id));
      renderChips(cond);
    }
    const anchor = $('panel-explore') || $('fba-conditions');
    if (anchor) anchor.scrollIntoView({ behavior: 'smooth' });
  };
  window.fbaSelectModel = (gemFile) => window.fbaSetModel('a', gemFile); // back-compat

  // URL handoff from the browser page: analysis.html?model=<gem>&slot=a
  const q = new URLSearchParams(location.search);
  if (q.get('model')) {
    const slot = q.get('slot') === 'b' ? 'b' : 'a';
    if (slot === 'b') setMode('compare');
    selectModel(S.conditions[slot], q.get('model'));
  }
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
else init();
