/* ── Shared media picker ──────────────────────────────────────────────────────
   One widget, mounted by every analysis, so the whole Media DB is available
   everywhere and not just in Explore. The media count is read from the database
   rather than written here, because a hardcoded one goes stale.

   A medium only becomes real once bound to a model, because the solver closes
   every exchange the medium does not name. So the widget owns the binding and
   always shows the coverage: how many compounds landed, how many ids had to be
   resolved across BiGG naming generations, and what is simply not in this model.
   ────────────────────────────────────────────────────────────────────────── */
import * as MediaDB from './media.js';
import { bindMedium } from './fba_engine.js';

const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c =>
  ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const pct = (x) => (100 * x).toFixed(0) + '%';

/* Presets, loaded once and shared. */
let PRESETS = null;
export async function presets() {
  if (!PRESETS) PRESETS = await (await fetch('fba/media_presets.json')).json();
  return PRESETS;
}
export function presetSpec(key) {
  const p = PRESETS[key];
  return { kind: 'preset', id: key, label: p.label, desc: p.desc, source: p.source,
           components: MediaDB.presetToComponents(p) };
}
export function dbSpec(m) {
  return { kind: 'db', id: m.id, label: m.name,
           source: (m.provenance && m.provenance.citation) || (m.provenance && m.provenance.source_type) || 'Media DB',
           desc: m.description || '', components: m.components };
}

/* Which medium suits which pangenome by default. Lactobacillaceae are fastidious
   and do not grow on a minimal medium at all: 0 of 10 LactoPanGEM models grow on
   M9 + glucose, 5 of 10 on MRS or BHI, and 16 of 16 on CDM. So CDM is the default
   for that collection. */
export const MINIMAL = new Set(['M9_glucose_aerobic', 'M9_glucose_anaerobic']);
export const DEFAULT_MEDIUM = { EcopanGEM: 'M9_glucose_aerobic', LactoPanGEM: 'CDM' };

/* ── The browse-all dialog ─────────────────────────────────────────────────── */
export async function openMediaBrowser(onPick) {
  const ov = document.createElement('div');
  ov.className = 'dlg-overlay';
  ov.innerHTML = `
    <div class="dlg" role="dialog" aria-modal="true" aria-label="Browse growth media">
      <div class="dlg-head">
        <div>
          <h5>Growth media</h5>
          <p class="dlg-sub"><b class="md-count">Curated</b> media, keyed to BiGG exchange reactions. DSMZ MediaDive, HMDB biospecimens, USDA foods, and media from GEM papers. Bounds bind to whichever model you have loaded.</p>
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
  const close = () => { ov.remove(); document.removeEventListener('keydown', onEsc); };
  const onEsc = (e) => { if (e.key === 'Escape') close(); };
  ov.querySelector('.dlg-x').addEventListener('click', close);
  ov.addEventListener('click', (e) => { if (e.target === ov) close(); });
  document.addEventListener('keydown', onEsc);

  const q = ov.querySelector('.md-q'), cat = ov.querySelector('.md-cat'), list = ov.querySelector('.md-list');
  let catalog;
  try { catalog = await MediaDB.catalog(); }
  catch (e) { list.innerHTML = `<div class="md-empty">Media DB unavailable: ${esc(e.message)}</div>`; return; }
  ov.querySelector('.md-count').textContent = catalog.length.toLocaleString();   // the catalog itself, not a guess

  const draw = () => {
    const res = MediaDB.search(catalog, q.value, { category: cat.value }, 80);
    list.innerHTML = res.length ? res.map(m => `
      <button type="button" class="md-item" data-id="${esc(m.id)}">
        <span class="md-name">${esc(m.name)}</span>
        <span class="md-meta">
          <span class="md-tag ${esc(m.category)}">${esc(m.category)}</span>
          <span>${m.n_components} compounds</span>
          ${m.aerobic === true ? '<span>aerobic</span>' : m.aerobic === false ? '<span>anaerobic</span>' : ''}
          ${m.source_db ? `<span class="md-src">${esc(m.source_db)}</span>` : ''}
        </span>
      </button>`).join('') : `<div class="md-empty">No media match.</div>`;
  };
  q.addEventListener('input', draw);
  cat.addEventListener('change', draw);
  draw(); q.focus();

  list.addEventListener('click', async (e) => {
    const b = e.target.closest('.md-item'); if (!b) return;
    b.classList.add('busy');
    try {
      const m = await MediaDB.medium(b.dataset.id);
      close();
      onPick(dbSpec(m));
    } catch (err) {
      b.classList.remove('busy');
      list.insertAdjacentHTML('afterbegin', `<div class="md-empty">Could not load: ${esc(err.message)}</div>`);
    }
  });
}

/* ── The picker widget ────────────────────────────────────────────────────────
   Mount next to any <select> that fillMedia() already populated. `state` is the
   tab's state object; the widget writes state.mediaSpec / .mediaBounds /
   .mediaReport onto it, and rebinds whenever the model or the medium changes.
   Call rebind(state) after a model loads. */
export function mountMediaPicker(selectEl, state, opts = {}) {
  const host = document.createElement('div');
  host.className = 'media-mount';
  host.innerHTML = `
    <button type="button" class="media-browse">
      <svg viewBox="0 0 24 24" width="14" height="14" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" aria-hidden="true"><circle cx="11" cy="11" r="7"/><path d="M20 20l-4.3-4.3"/></svg>
      Browse <span class="mc">all</span> media
    </button>
    <div class="media-report"></div>
    <label class="media-min"><input type="checkbox" class="me-minerals"> open essential inorganic ions and water</label>`;
  const anchor = selectEl.nextElementSibling && selectEl.nextElementSibling.classList.contains('fba-media-desc')
    ? selectEl.nextElementSibling : selectEl;
  anchor.parentNode.insertBefore(host, anchor.nextSibling);
  state._mediaHost = host;
  MediaDB.fillCount(host.querySelector('.mc'), n => `all ${n.toLocaleString()}`);

  state.mediaSpec = presetSpec(selectEl.value);
  state.openMinerals = false;

  selectEl.addEventListener('change', () => {
    state.mediaSpec = presetSpec(selectEl.value);
    state.mediaTouched = true;
    rebind(state);
  });
  host.querySelector('.media-browse').addEventListener('click', () =>
    openMediaBrowser((spec) => {
      state.mediaSpec = spec;
      state.mediaTouched = true;
      selectEl.selectedIndex = -1;
      rebind(state);
    }));
  host.querySelector('.me-minerals').addEventListener('change', (e) => {
    state.openMinerals = e.target.checked;
    rebind(state);
  });
  renderReport(state);
  return host;
}

/** Bind state.mediaSpec to state.model and refresh the coverage line. */
export function rebind(state) {
  if (!state.mediaSpec) return;
  if (state.model) {
    const r = bindMedium(state.model, state.mediaSpec.components, { openMinerals: state.openMinerals });
    state.mediaBounds = r.bounds;
    state.mediaReport = r;
  }
  renderReport(state);
}

/** Called when a model finishes loading: pick a medium the organism can grow on. */
export function adoptDefaultMedium(state, selectEl) {
  const ds = (state.meta || {}).dataset;
  const want = DEFAULT_MEDIUM[ds];
  if (want && PRESETS[want] && !state.mediaTouched && state.mediaSpec
      && state.mediaSpec.kind === 'preset' && MINIMAL.has(state.mediaSpec.id)
      && want !== state.mediaSpec.id) {
    state.mediaSpec = presetSpec(want);
    if (selectEl) selectEl.value = want;
    state._swapped = PRESETS[want].label;
  }
  rebind(state);
}

function renderReport(state) {
  const host = state._mediaHost; if (!host) return;
  const box = host.querySelector('.media-report');
  const spec = state.mediaSpec;
  if (!spec) { box.innerHTML = ''; return; }
  if (!state.model) {
    box.innerHTML = `<div class="mr-line"><b>${esc(spec.label)}</b> · ${spec.components.length} compounds
      <span class="mr-hint">pick a model to bind it</span></div>`;
    return;
  }
  const r = state.mediaReport || { mapped: [], missing: [], added: [], coverage: 0 };
  const n = spec.components.length;
  const cls = r.coverage >= 0.9 ? 'ok' : r.coverage >= 0.6 ? 'warn' : 'bad';
  const renamed = r.mapped.filter(m => m.renamed).length;
  box.innerHTML = `
    <div class="mr-line">
      <b>${esc(spec.label)}</b>
      <span class="mr-cov ${cls}">${r.mapped.length}/${n} bound · ${pct(r.coverage)}</span>
    </div>
    <div class="mr-sub">
      ${renamed ? `<span class="mr-fix" title="Resolved across BiGG naming generations, e.g. EX_glc_D_e to EX_glc__D_e">${renamed} id${renamed > 1 ? 's' : ''} resolved</span>` : ''}
      ${r.added.length ? `<span class="mr-hint">+${r.added.length} minerals opened</span>` : ''}
      ${r.missing.length ? `<button type="button" class="mr-missing">${r.missing.length} not in this model</button>` : ''}
      ${spec.source ? `<span class="mr-src">${esc(String(spec.source).slice(0, 90))}</span>` : ''}
    </div>
    <div class="mr-list" hidden>${r.missing.slice(0, 60).map(m =>
      `<code>${esc(m.exchange)}</code> <span>${esc(m.name || '')}</span>`).join('')}</div>
    ${state._swapped ? `<div class="mr-swap">Switched to <b>${esc(state._swapped)}</b>: Lactobacillaceae do not grow on minimal media.</div>` : ''}`;
  const btn = box.querySelector('.mr-missing');
  if (btn) btn.addEventListener('click', () => {
    const l = box.querySelector('.mr-list'); l.hidden = !l.hidden;
  });
}
