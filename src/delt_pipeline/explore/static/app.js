// DELT Explorer: browse triaged DEL hits, building blocks and the resynthesis shortlist in the browser.
// All data comes from the inline JSON payload (see delt_pipeline/explore/payload.py); structures and
// substructure search use RDKit.js, loaded from a CDN. Without it the page still filters, plots and exports.
(() => {
  'use strict';

  const P = JSON.parse(document.getElementById('payload').textContent);
  const columnIndex = (columns) => Object.fromEntries(columns.map((name, i) => [name, i]));
  const H = columnIndex(P.hit_columns);
  const B = columnIndex(P.block_columns);
  const BS = columnIndex(P.blockset_columns);
  const hits = P.hits;
  const hitIndex = new Map(hits.map((row, i) => [row[H.hit_id], i]));
  const blockById = new Map(P.blocks.map((row) => [row[B.block_id], row]));
  const byScreen = Object.fromEntries(P.screens.map((s) => [s, []]));
  hits.forEach((row, i) => byScreen[row[H.screen]].push(i));

  const FP_BYTES = 256;
  const fingerprints = P.fingerprints.map((text) => Uint8Array.from(atob(text), (c) => c.charCodeAt(0)));
  const POPCOUNT = Uint8Array.from({ length: 256 }, (_, n) => n.toString(2).replaceAll('0', '').length);

  const TIERS = ['1', '2', '3', 'flagged'];
  const TIER_COLOURS = { 1: '#2a78d6', 2: '#1baf7a', 3: '#9a9890', flagged: '#eb6834' };
  const TIER_LABELS = { 1: 'Tier 1', 2: 'Tier 2', 3: 'Tier 3', flagged: 'Flagged (counter or bead)' };
  const FLAGS = ['counter_binder', 'matrix_binder', 'shared_target', 'replicate_support', 'disynthon_support'];
  const PAGE = 100;
  const SHOWN_BLOCKS = 24;
  const blockCount = { a: 1, b: 1, c: 1 };
  P.blocks.forEach((row) => {
    blockCount[row[B.position]] = Math.max(blockCount[row[B.position]], row[B.bb_id] + 1);
  });

  const state = {
    screen: P.screens.includes('mapk14') ? 'mapk14' : P.screens[0],
    tiers: new Set(TIERS),
    flag: '',
    minLB: 0,
    minCount: 0,
    picksOnly: false,
    marked: new Set(),
    current: null,
    sort: { key: 'rank', dir: 1 },
    page: 0,
    search: { mode: 'sub', text: '', blocks: null, sims: null, minSim: 0.4, error: '' },
  };

  const $ = (id) => document.getElementById(id);
  const esc = (text) => String(text ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
  const num = (x, digits = 3) => (x === null || x === undefined ? '–' : Number(x.toPrecision(digits)).toString());
  const blockId = (position, id) => position.toUpperCase() + String(id).padStart(4, '0');
  const blocksOf = (row) => [blockId('a', row[H.bb_a]), blockId('b', row[H.bb_b]), blockId('c', row[H.bb_c])];
  const textColour = () => getComputedStyle(document.body).color;

  // ---- RDKit.js (optional) ------------------------------------------------------------------------
  let RDKit = null;
  const molCache = new Map();

  function loadRDKit() {
    return new Promise((resolve) => {
      const script = document.createElement('script');
      script.src = window.DELT_CONFIG.rdkitJs;
      script.onload = () => {
        window.initRDKitModule({ locateFile: () => window.DELT_CONFIG.rdkitWasm }).then(resolve, () => resolve(null));
      };
      script.onerror = () => resolve(null);
      document.head.appendChild(script);
    });
  }

  function getMol(smiles) {
    if (!RDKit) return null;
    if (!molCache.has(smiles)) {
      let mol = null;
      try {
        mol = RDKit.get_mol(smiles);
        if (mol && !mol.is_valid()) {
          mol.delete();
          mol = null;
        }
      } catch (error) {
        mol = null;
      }
      molCache.set(smiles, mol);
    }
    return molCache.get(smiles);
  }

  function structure(smiles, width, height) {
    const mol = getMol(smiles);
    return mol ? `<span class="mol">${mol.get_svg(width, height)}</span>` : `<code>${esc(smiles)}</code>`;
  }

  // ---- search -------------------------------------------------------------------------------------
  function packBits(bitString) {
    const bytes = new Uint8Array(FP_BYTES);
    for (let i = 0; i < bitString.length; i++) if (bitString[i] === '1') bytes[i >> 3] |= 0x80 >> (i & 7);
    return bytes;
  }

  function tanimoto(a, b) {
    let common = 0;
    let total = 0;
    for (let i = 0; i < FP_BYTES; i++) {
      common += POPCOUNT[a[i] & b[i]];
      total += POPCOUNT[a[i]] + POPCOUNT[b[i]];
    }
    return total - common === 0 ? 0 : common / (total - common);
  }

  function substructureSearch(text) {
    const query = RDKit.get_qmol(text);
    if (!query || !query.is_valid()) throw new Error('not a valid SMARTS or SMILES');
    const matched = new Set();
    for (const row of P.blocks) {
      const mol = getMol(row[B.fragment_smiles]);
      if (mol && mol.get_substruct_match(query) !== '{}') matched.add(row[B.block_id]);
    }
    query.delete();
    return matched;
  }

  function similaritySearch(text) {
    const mol = RDKit.get_mol(text);
    if (!mol || !mol.is_valid()) throw new Error('not a valid SMILES');
    const query = packBits(mol.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 })));
    mol.delete();
    return new Map(fingerprints.map((print, i) => [i, tanimoto(query, print)]));
  }

  function runSearch() {
    const search = state.search;
    search.mode = $('mode').value;
    search.text = $('query').value.trim();
    search.minSim = Number($('min-sim').value) || 0.4;
    search.blocks = null;
    search.sims = null;
    search.error = '';
    if (search.text) {
      if (!RDKit) {
        search.error = 'Search needs RDKit.js, which could not be loaded (no network?).';
      } else {
        try {
          if (search.mode === 'sub') search.blocks = substructureSearch(search.text);
          else search.sims = similaritySearch(search.text);
        } catch (error) {
          search.error = `Cannot search for "${search.text}": ${error.message}.`;
        }
      }
    }
    state.page = 0;
    render();
  }

  // ---- filtering and sorting ----------------------------------------------------------------------
  function visibleIndices() {
    const { search } = state;
    const out = [];
    for (const i of byScreen[state.screen]) {
      const row = hits[i];
      if (!state.tiers.has(row[H.tier])) continue;
      if (row[H.poisson_lower] < state.minLB) continue;
      if (row[H.target_sum] < state.minCount) continue;
      if (state.picksOnly && !row[H.pick_rank]) continue;
      if (state.flag && !row[H[state.flag]]) continue;
      if (search.blocks && !blocksOf(row).some((id) => search.blocks.has(id))) continue;
      if (search.sims && !(search.sims.get(i) >= search.minSim)) continue;
      out.push(i);
    }
    return out;
  }

  const SORT_VALUE = {
    pick: (i) => hits[i][H.pick_rank] || Infinity,
    sim: (i) => state.search.sims?.get(i) ?? 0,
  };

  function sorted(indices) {
    const { key, dir } = state.sort;
    const value = SORT_VALUE[key] ?? ((i) => hits[i][H[key]]);
    return [...indices].sort((a, b) => {
      const x = value(a);
      const y = value(b);
      if (x === y) return hits[a][H.rank] - hits[b][H.rank];
      if (x === null) return 1;
      if (y === null) return -1;
      return (x < y ? -1 : 1) * dir;
    });
  }

  // ---- rendering ----------------------------------------------------------------------------------
  const hoverText = (i) => {
    const row = hits[i];
    return `${esc(row[H.hit_id].slice(0, 22))}<br>rank ${row[H.rank]} · ${blocksOf(row).join(' ')}<br>` +
      `target ${row[H.target_sum]} · control ${row[H.control_sum]} · LB ${num(row[H.poisson_lower])}<br>${esc(row[H.reasons])}`;
  };

  const chartLayout = (extra) => ({
    margin: { l: 56, r: 12, t: 28, b: 40 },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: textColour(), size: 12 },
    legend: { orientation: 'h', y: -0.18 },
    ...extra,
  });

  const plotConfig = { displaylogo: false, responsive: true };
  let scatterIndex = [];
  let cubeIndex = [];

  function tierTraces(vis, make) {
    const index = [];
    const traces = [];
    for (const tier of TIERS) {
      const members = vis.filter((i) => hits[i][H.tier] === tier);
      if (!members.length) continue;
      traces.push({ name: TIER_LABELS[tier], customdata: members, ...make(members, TIER_COLOURS[tier]) });
      index.push(members);
    }
    const marked = vis.filter((i) => state.marked.has(hits[i][H.hit_id]));
    traces.push({ name: 'Marked', hoverinfo: 'skip', ...make(marked, textColour(), true) });
    index.push(marked);
    return { traces, index };
  }

  function renderScatter(vis) {
    const { traces, index } = tierTraces(vis, (members, colour, ring) => ({
      type: 'scattergl',
      mode: 'markers',
      x: members.map((i) => hits[i][H.target_sum]),
      y: members.map((i) => hits[i][H.poisson_lower]),
      text: members.map(hoverText),
      hovertemplate: ring ? undefined : '%{text}<extra></extra>',
      marker: ring
        ? { symbol: 'circle-open', size: 14, color: colour, line: { color: colour, width: 2 } }
        : { color: colour, size: 6, opacity: 0.75 },
      selected: { marker: { opacity: ring ? 1 : 0.75 } },
      unselected: { marker: { opacity: ring ? 1 : 0.75 } },
    }));
    scatterIndex = index;
    Plotly.react(
      'scatter',
      traces,
      chartLayout({
        xaxis: { type: 'log', title: { text: 'Target count (replicates pooled)' } },
        yaxis: { type: 'log', title: { text: 'Poisson lower bound of the rate ratio' } },
        dragmode: 'lasso',
        uirevision: `scatter-${state.screen}`,
      }),
      plotConfig,
    );
  }

  function renderCube(vis) {
    const { traces, index } = tierTraces(vis, (members, colour, ring) => ({
      type: 'scatter3d',
      mode: 'markers',
      x: members.map((i) => hits[i][H.bb_a]),
      y: members.map((i) => hits[i][H.bb_b]),
      z: members.map((i) => hits[i][H.bb_c]),
      text: members.map(hoverText),
      hovertemplate: ring ? undefined : '%{text}<extra></extra>',
      marker: ring
        ? { symbol: 'circle-open', size: 7, color: colour, line: { color: colour, width: 3 } }
        : { color: colour, size: 3, opacity: 0.8 },
    }));
    cubeIndex = index;
    const axis = (title, position) => ({ title: { text: title }, range: [0, blockCount[position]] });
    Plotly.react(
      'cube',
      traces,
      chartLayout({
        margin: { l: 0, r: 0, t: 28, b: 0 },
        scene: { xaxis: axis('Building block A', 'a'), yaxis: axis('Building block B', 'b'), zaxis: axis('Building block C', 'c') },
        uirevision: `cube-${state.screen}`,
      }),
      plotConfig,
    );
  }

  const COLUMNS = [
    ['pick', 'Pick'],
    ['rank', 'Rank'],
    ['tier', 'Tier'],
    ['bb_a', 'A'],
    ['bb_b', 'B'],
    ['bb_c', 'C'],
    ['disynthon', 'Disynthon'],
    ['target_sum', 'Target'],
    ['control_sum', 'Control'],
    ['poisson_lower', 'Poisson LB'],
    ['eb_log2_ratio', 'EB log2'],
  ];

  function renderTable(vis) {
    const rows = sorted(vis);
    const pages = Math.max(1, Math.ceil(rows.length / PAGE));
    state.page = Math.min(state.page, pages - 1);
    const shown = rows.slice(state.page * PAGE, (state.page + 1) * PAGE);
    const showSim = Boolean(state.search.sims);
    const head = COLUMNS.concat(showSim ? [['sim', 'Similarity']] : [])
      .map(([key, label]) => `<th data-sort="${key}" aria-sort="${state.sort.key === key ? (state.sort.dir > 0 ? 'ascending' : 'descending') : 'none'}">${label}${state.sort.key === key ? (state.sort.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`)
      .join('');
    const body = shown.map((i) => {
      const row = hits[i];
      const id = row[H.hit_id];
      const marked = state.marked.has(id);
      return `<tr class="hit${id === state.current ? ' current' : ''}${marked ? ' marked' : ''}" data-id="${esc(id)}">` +
        `<td><input type="checkbox" class="mark" aria-label="Mark ${esc(id)}"${marked ? ' checked' : ''}></td>` +
        `<td>${row[H.pick_rank] || ''}</td><td>${row[H.rank]}</td><td>${esc(row[H.tier])}</td>` +
        blocksOf(row).map((b) => `<td>${b}</td>`).join('') +
        `<td>${esc(row[H.disynthon] ?? '')}</td><td>${row[H.target_sum]}</td><td>${row[H.control_sum]}</td>` +
        `<td>${num(row[H.poisson_lower])}</td><td>${num(row[H.eb_log2_ratio])}</td>` +
        (showSim ? `<td>${num(state.search.sims.get(i), 2)}</td>` : '') + '</tr>';
    }).join('');
    $('hits').innerHTML = `<thead><tr><th></th>${head}</tr></thead><tbody>${body}</tbody>`;
    $('page-info').textContent = `Page ${state.page + 1} of ${pages} · ${rows.length} hits`;
    $('prev').disabled = state.page === 0;
    $('next').disabled = state.page >= pages - 1;
  }

  function blockSummary(id) {
    const row = blockById.get(id);
    if (!row) return `<p>${id}</p>`;
    const stats = P.screens.concat(Object.keys(P.blocksets).filter((s) => !P.screens.includes(s)))
      .map((screen) => {
        const entry = P.blocksets[screen]?.[id];
        const label = screen.toUpperCase();
        return entry ? `${label} ${num(entry[BS.auc], 3)}${screen === state.screen ? '*' : ''}` : '';
      })
      .filter(Boolean)
      .join(' · ');
    return `<div class="card">${structure(row[B.fragment_smiles], 170, 120)}<br><b>${id}</b>` +
      `${row[B.n_same_fragment] > 1 ? ` <span class="tag">${row[B.n_same_fragment]} blocks, same fragment</span>` : ''}` +
      `${row[B.low_load] ? ' <span class="tag on">low load</span>' : ''}<br>` +
      `<span class="muted">MW ${num(row[B.mw], 4)} · logP ${num(row[B.clogp], 2)}</span><br>` +
      `<span class="muted" title="Set-enrichment AUC of the block (0.5 = none)">AUC ${stats || '–'}</span></div>`;
  }

  function renderBlocks() {
    const panel = $('blocks-panel');
    const matched = state.search.blocks;
    panel.hidden = !matched;
    if (!matched) return;
    const auc = (id) => P.blocksets[state.screen]?.[id]?.[BS.auc] ?? 0;
    const ids = [...matched].sort((a, b) => auc(b) - auc(a) || (a < b ? -1 : 1));
    $('blocks-summary').textContent = `${matched.size} of ${P.blocks.length} building blocks contain "${state.search.text}". ` +
      `Showing the ${Math.min(ids.length, SHOWN_BLOCKS)} with the highest set-enrichment AUC in ${state.screen.toUpperCase()}; ` +
      'the AUC of every screen is listed (* = current), including counter-screens.';
    $('blocks').innerHTML = ids.slice(0, SHOWN_BLOCKS).map(blockSummary).join('');
  }

  function renderDetail() {
    const panel = $('detail');
    const i = hitIndex.get(state.current);
    if (i === undefined) {
      panel.innerHTML = '<h2>Compound</h2><p class="muted">Click a row or a point.</p>';
      return;
    }
    const row = hits[i];
    const flags = FLAGS.map((flag) => `<span class="tag${row[H[flag]] && (flag === 'counter_binder' || flag === 'matrix_binder') ? ' on' : ''}">${row[H[flag]] ? '✓' : '✗'} ${flag.replaceAll('_', ' ')}</span>`).join('');
    const counts = (names) => names.map((name) => `<td>${row[H[name]]}</td>`).join('');
    panel.innerHTML = `<h2>Compound</h2>${structure(row[H.smiles], 280, 200)}` +
      `<p><code>${esc(row[H.molecule_hash].slice(0, 16))}</code> · rank ${row[H.rank]} · tier ${esc(row[H.tier])}` +
      `${row[H.pick_rank] ? ` · <b>shortlist pick ${row[H.pick_rank]}</b> (cluster ${row[H.cluster_id]})` : ''}</p>` +
      `<p>${esc(row[H.reasons])}</p><div>${flags}</div>` +
      '<h3>Counts</h3><table><tr><th></th><th>1</th><th>2</th><th>3</th></tr>' +
      `<tr><th>Target</th>${counts(['target_1', 'target_2', 'target_3'])}</tr>` +
      `<tr><th>Matrix</th>${counts(['matrix_1', 'matrix_2', 'matrix_3'])}</tr></table>` +
      `<p class="muted">Pre-selection load: ${row[H.load]}</p>` +
      '<h3>Scores</h3><table>' +
      `<tr><th>Poisson ratio</th><td>${num(row[H.poisson_ratio])} (${num(row[H.poisson_lower])} to ${num(row[H.poisson_upper])})</td></tr>` +
      `<tr><th>Enrichment ratio</th><td>${num(row[H.er])}</td></tr>` +
      `<tr><th>EB log2 ratio</th><td>${num(row[H.eb_log2_ratio])}</td></tr>` +
      `<tr><th>Faver z</th><td>${num(row[H.z_faver])}</td></tr>` +
      `<tr><th>Count difference</th><td>${num(row[H.count_diff])}</td></tr></table>` +
      `<h3>Building blocks</h3><div class="cards">${blocksOf(row).map(blockSummary).join('')}</div>`;
  }

  function renderStatus(vis) {
    const { search } = state;
    const parts = [`${vis.length} of ${byScreen[state.screen].length} hits shown`, `${state.marked.size} marked`];
    if (search.blocks) parts.push(`${search.blocks.size} matching building blocks`);
    if (!RDKit) parts.push(window.DELT.rdkitState === 'loading' ? 'loading structure tools…' : 'structures unavailable (RDKit.js not loaded)');
    $('status').className = search.error ? 'status error' : 'status';
    $('status').textContent = search.error || parts.join(' · ');
  }

  function render() {
    const vis = visibleIndices();
    renderStatus(vis);
    renderScatter(vis);
    renderCube(vis);
    renderTable(vis);
    renderBlocks();
    renderDetail();
  }

  // ---- interaction --------------------------------------------------------------------------------
  function mark(ids, replace) {
    if (replace) state.marked.clear();
    ids.forEach((id) => state.marked.add(id));
    render();
  }

  function toggleMark(id) {
    if (state.marked.has(id)) state.marked.delete(id);
    else state.marked.add(id);
    render();
  }

  const idsOf = (points, index) => points
    .filter((point) => point.curveNumber < index.length - 1) // the last trace is the ring of the marked hits
    .map((point) => hits[index[point.curveNumber][point.pointNumber]][H.hit_id]);

  function wirePlots() {
    $('scatter').on('plotly_selected', (event) => {
      if (event) mark(idsOf(event.points, scatterIndex), true);
    });
    $('scatter').on('plotly_deselect', () => mark([], true));
    $('scatter').on('plotly_click', (event) => {
      state.current = idsOf(event.points.slice(0, 1), scatterIndex)[0] ?? state.current;
      render();
    });
    $('cube').on('plotly_click', (event) => {
      const [id] = idsOf(event.points.slice(0, 1), cubeIndex);
      if (id) {
        state.current = id;
        toggleMark(id);
      }
    });
  }

  function csvCell(value) {
    const text = value === null || value === undefined ? '' : String(value);
    return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
  }

  function exportMarked() {
    const rows = [...state.marked].map((id) => hits[hitIndex.get(id)]);
    const lines = [P.hit_columns.join(',')].concat(rows.map((row) => row.map(csvCell).join(',')));
    const link = document.createElement('a');
    link.href = URL.createObjectURL(new Blob([`${lines.join('\n')}\n`], { type: 'text/csv' }));
    link.download = `delt_marked_${state.screen}.csv`;
    link.click();
    URL.revokeObjectURL(link.href);
  }

  function setup() {
    $('screen').innerHTML = P.screens.map((s) => `<option value="${esc(s)}">${esc(s.toUpperCase())}</option>`).join('');
    $('screen').value = state.screen;
    $('tiers').innerHTML = TIERS.map((tier) => `<label class="inline"><input type="checkbox" value="${tier}" checked> ${tier === 'flagged' ? 'Flagged' : `Tier ${tier}`}</label>`).join('');
    $('flag').innerHTML = '<option value="">any</option>' + FLAGS.map((f) => `<option value="${f}">${f.replaceAll('_', ' ')}</option>`).join('');
    const resetCount = () => {
      const top = Math.max(...byScreen[state.screen].map((i) => hits[i][H.target_sum]));
      $('min-count').max = top;
      $('min-count').value = 0;
      $('min-count-value').textContent = '0';
      state.minCount = 0;
    };
    resetCount();

    const refresh = () => {
      state.page = 0;
      render();
    };
    $('screen').addEventListener('change', (event) => {
      state.screen = event.target.value;
      state.marked.clear();
      state.current = null;
      state.search.blocks = null;
      state.search.sims = null;
      resetCount();
      refresh();
    });
    $('tiers').addEventListener('change', () => {
      state.tiers = new Set([...$('tiers').querySelectorAll('input:checked')].map((box) => box.value));
      refresh();
    });
    $('flag').addEventListener('change', (event) => {
      state.flag = event.target.value;
      refresh();
    });
    $('min-lb').addEventListener('input', (event) => {
      state.minLB = Number(event.target.value) || 0;
      refresh();
    });
    $('min-count').addEventListener('input', (event) => {
      state.minCount = Number(event.target.value);
      $('min-count-value').textContent = event.target.value;
      refresh();
    });
    $('picks-only').addEventListener('change', (event) => {
      state.picksOnly = event.target.checked;
      refresh();
    });
    $('search').addEventListener('click', runSearch);
    $('query').addEventListener('keydown', (event) => {
      if (event.key === 'Enter') runSearch();
    });
    $('clear-search').addEventListener('click', () => {
      $('query').value = '';
      runSearch();
    });
    $('clear-marks').addEventListener('click', () => mark([], true));
    $('export').addEventListener('click', exportMarked);
    $('prev').addEventListener('click', () => {
      state.page -= 1;
      render();
    });
    $('next').addEventListener('click', () => {
      state.page += 1;
      render();
    });
    $('hits').addEventListener('click', (event) => {
      const header = event.target.closest('th[data-sort]');
      if (header) {
        const key = header.dataset.sort;
        state.sort = { key, dir: state.sort.key === key ? -state.sort.dir : 1 };
        render();
        return;
      }
      const tr = event.target.closest('tr.hit');
      if (!tr) return;
      if (event.target.classList.contains('mark')) toggleMark(tr.dataset.id);
      else {
        state.current = tr.dataset.id;
        render();
      }
    });
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);
  }

  // The test hooks are read-only views of the state.
  window.DELT = {
    state,
    rdkitState: 'loading',
    visible: () => visibleIndices().map((i) => hits[i][H.hit_id]),
    fingerprintOf: (smiles) => {
      const mol = RDKit.get_mol(smiles);
      const bits = packBits(mol.get_morgan_fp(JSON.stringify({ radius: 2, nBits: 2048 })));
      mol.delete();
      return Array.from(bits);
    },
  };

  setup();
  render();
  wirePlots();
  window.DELT.ready = loadRDKit().then((module) => {
    RDKit = module;
    window.DELT.rdkitState = module ? 'ready' : 'failed';
    render();
  });
})();
