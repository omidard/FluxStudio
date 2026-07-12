/* ── Media DB client ──────────────────────────────────────────────────────────
   Reads the Media dataset (github.com/omidard/Media), a read-only static API on
   GitHub Pages with permissive CORS. Media are keyed by BiGG exchange
   reactions, with the same sign convention the solver uses: lower_bound < 0 is
   uptake, in mmol gDW-1 h-1.

     GET /data/index.json        catalog (444 KB gzipped, fetched once, cached)
     GET /data/media/{id}.json   one medium, with components[]
     GET /data/stats.json        counts only, a few hundred bytes

   Sources include DSMZ MediaDive (3,148 laboratory media), USDA FoodData
   Central, FooDB, and media curated from GEM papers and GrowthDB.
   ────────────────────────────────────────────────────────────────────────── */

export const MEDIA_BASE = 'https://omidard.github.io/Media';

let _catalog = null;
let _catalogPromise = null;
const _mediumCache = new Map();

/** Full catalog: [{id, name, category, organism_scope, aerobic, n_components, source_db, ...}] */
export async function catalog() {
  if (_catalog) return _catalog;
  if (!_catalogPromise) {
    _catalogPromise = fetch(`${MEDIA_BASE}/data/index.json`)
      .then(r => { if (!r.ok) throw new Error('Media catalog unavailable'); return r.json(); })
      .then(d => {
        _catalog = (d.media || []).map(m => ({
          ...m,
          _hay: `${m.id} ${m.name} ${m.source_db || ''} ${m.category || ''}`.toLowerCase(),
        }));
        _catalog.meta = { count: d.count, by_category: d.by_category, by_source_db: d.by_source_db };
        return _catalog;
      });
  }
  return _catalogPromise;
}

/* Counts, from their own endpoint, because the catalog is 6.5 MB and a label should
   not have to pull all of it. These were hardcoded before and went stale the moment
   a medium was added. Never blocks the UI: on failure the caller shows no number. */
let _statsPromise = null;
export async function stats() {
  if (!_statsPromise) {
    _statsPromise = fetch(`${MEDIA_BASE}/data/stats.json`)
      .then(r => r.ok ? r.json() : null)
      .catch(() => null);                      // offline, or the endpoint predates this
  }
  return _statsPromise;
}

/** Write the live media count into an element, once it is known. */
export function fillCount(el, fmt = (n) => n.toLocaleString()) {
  if (!el) return;
  stats().then(s => {
    const n = s ? s.count : (_catalog ? _catalog.length : null);
    if (n) el.textContent = fmt(n);
  });
}

/** Append "(4,198)" style counts to the category <option>s of a filter select. */
export function fillCategoryCounts(selectEl) {
  if (!selectEl) return;
  stats().then(s => {
    if (!s || !s.by_category) return;
    for (const o of selectEl.options) {
      const n = s.by_category[o.value];
      if (n) o.textContent = `${o.textContent.replace(/\s*\(.*\)$/, '')} (${n.toLocaleString()})`;
    }
  });
}

/** One medium, with components[]. Cached. */
export async function medium(id) {
  if (_mediumCache.has(id)) return _mediumCache.get(id);
  const r = await fetch(`${MEDIA_BASE}/data/media/${encodeURIComponent(id)}.json`);
  if (!r.ok) throw new Error(`Medium "${id}" not found`);
  const m = await r.json();
  // keep only components that actually name an exchange reaction
  m.components = (m.components || []).filter(c => c.exchange);
  _mediumCache.set(id, m);
  return m;
}

/** Rank-and-filter the catalog. Exact id, then name-prefix, then substring. */
export function search(list, q, filt = {}, limit = 60) {
  const query = (q || '').trim().toLowerCase();
  let out = list;

  if (filt.category) out = out.filter(m => m.category === filt.category);
  if (filt.oxygen === 'aerobic') out = out.filter(m => m.aerobic === true);
  if (filt.oxygen === 'anaerobic') out = out.filter(m => m.aerobic === false);
  if (filt.source) out = out.filter(m => (m.source_db || '') === filt.source);

  if (!query) return out.slice(0, limit);

  /* Match every word, not one contiguous string. Media names carry punctuation
     ("MRS (de Man, Rogosa, Sharpe) medium"), so a typed "MRS de Man" would never
     be a substring of it. */
  const terms = query.split(/\s+/).filter(Boolean);
  const scored = [];
  for (const m of out) {
    if (!terms.every(t => m._hay.includes(t))) continue;
    const name = (m.name || '').toLowerCase();
    let s = 4;
    if (m.id === query) s = 0;
    else if (name === query) s = 1;
    else if (name.startsWith(terms[0])) s = 2;
    else if (name.includes(query)) s = 3;
    scored.push([s, m]);
    if (scored.length > 4000) break;             // the catalog is large; do not scan forever
  }
  scored.sort((a, b) => a[0] - b[0] || (b[1].n_components || 0) - (a[1].n_components || 0));
  return scored.slice(0, limit).map(x => x[1]);
}

/** A handful of media worth surfacing before the user searches. */
export const FEATURED = [
  'bhi', 'biolog_if0_minimal', 'biospecimen_hmdb_blood', 'biospecimen_hmdb_urine',
  'biospecimen_hmdb_feces', 'biospecimen_hmdb_saliva', 'blood_agar_columbia',
];

/** Turn a legacy fba/media_presets.json entry into the Media-DB component shape,
    so the preset path and the database path go through the exact same binder. */
export function presetToComponents(preset) {
  return Object.entries(preset.bounds || {}).map(([exchange, lower_bound]) => ({
    exchange, lower_bound, upper_bound: 1000, name: exchange.replace(/^EX_/, '').replace(/_e$/, ''),
  }));
}
