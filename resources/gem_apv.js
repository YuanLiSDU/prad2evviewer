// gem_apv.js — GEM APV waveform viewer tab
//
// Stacked per-GEM sections (one per detector, ordered by det_id), each
// with a faint background tint and a thin separator above it.  Inside
// each section is a responsive grid of small canvas panels — one per
// APV — drawing the 6 time-sample traces (blue → red), zero line,
// optional ±N-sigma threshold band, optional firmware-CM overlay, and
// two stacked tick rows: dim "fw_hits" (firmware survivors) above
// bright "hits" (software-cut survivors).
//
// Per-APV pedestal noise comes from /api/gem/calib (one-shot, cached;
// refetched only when the per-event calib_rev diverges from the cached
// value).  The σ multiplier is editable from the toolbar — POST goes to
// /api/gem/threshold and applies to the live reconstruction shared by
// every consumer (other browsers, histograms, clustering downstream).
//
// Mirrors the simplified controls from gem_event_viewer.py's RawApvTab:
// Process / Signal Only / Shared Y / FW Hits / Threshold / σ / CM /
// per-sample (t0…t5).  Clustering knobs are intentionally left out —
// the monitor is for at-a-glance live inspection.

'use strict';

// Per-GEM tints — same palette as gem_event_viewer.py's "All" tab so
// the desktop and web monitor agree on "GEM 1 is green".
const GEM_APV_TINTS = [
    'rgba(0, 180, 216, 0.13)',
    'rgba(81, 207, 102, 0.13)',
    'rgba(255, 146, 43, 0.13)',
    'rgba(204, 93, 232, 0.13)',
];
function gemApvTint(detId) {
    if (detId < 0) return 'transparent';
    return GEM_APV_TINTS[detId % GEM_APV_TINTS.length];
}

// Time-sample trace colours: HSV blue→red, matches ApvPanel._paint.
const GEM_APV_TS_COLORS = (() => {
    const out = [];
    for (let t = 0; t < 6; t++) {
        const frac = t / 5;
        // hue 0.66 (blue) → 0 (red)
        const h = 0.66 * (1 - frac);
        out.push(hsv2rgb(h, 0.85, 0.95));
    }
    return out;
})();
// CM overlay colours — same hue as the trace but desaturated so the
// firmware CM line reads as related to the matching time sample without
// fighting the data trace for visibility.
const GEM_APV_CM_COLORS = (() => {
    const out = [];
    for (let t = 0; t < 6; t++) {
        const frac = t / 5;
        const h = 0.66 * (1 - frac);
        out.push(hsv2rgb(h, 0.6, 1.0));
    }
    return out;
})();
function hsv2rgb(h, s, v) {
    const i = Math.floor(h * 6);
    const f = h * 6 - i;
    const p = v * (1 - s);
    const q = v * (1 - f * s);
    const t = v * (1 - (1 - f) * s);
    let r, g, b;
    switch (i % 6) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        case 5: r = v; g = p; b = q; break;
    }
    return `rgb(${Math.round(r*255)},${Math.round(g*255)},${Math.round(b*255)})`;
}

// localStorage helpers for tab preferences — mirror theme.js's pattern.
// Wrapped in try/catch so private-mode browsers don't TDZ the tab.
function gemApvPrefGet(key, dflt) {
    try {
        const v = localStorage.getItem('prad2.gem_apv.' + key);
        return v === null ? dflt : v;
    } catch (e) { return dflt; }
}
function gemApvPrefSet(key, v) {
    try { localStorage.setItem('prad2.gem_apv.' + key, String(v)); } catch (e) {}
}

// Tab state.
let gemApvData = null;          // last fetched per-event payload
let gemApvCurrentEvent = -1;
// Cached calibration: { rev, zs_sigma, noise: Map<id, Float32Array(128)> }.
// Populated lazily by ensureGemApvCalib(); refetched when a per-event
// payload arrives with calib_rev != gemApvCalib.rev.
let gemApvCalib = null;
let gemApvCalibInflight = false;
let gemApvShowProcessed = true;
let gemApvShowSignalOnly = false;
let gemApvSharedY = true;
let gemApvShowThreshold = false;
let gemApvShowFwHits = true;
let gemApvShowCm = false;
let gemApvSampleMask = [true, true, true, true, true, true];
// Per-detector visibility — index = det_id (0..3 cover all current PRad-II
// GEMs).  Out-of-range det_ids fall back to "show" so unexpected
// configurations don't disappear silently.
let gemApvDetMask = [true, true, true, true];

// Layout — 'sequential' (default; prad1 style, 128×6 = 768 points laid out
// in time order) or 'overlay' (6 TS stacked on the same channel axis).
// Defaults to the RC monitoring view so Ashot/Kondo see the expected plot
// without flipping a toggle; operators can still switch via the Layout
// select and the choice persists in localStorage.
let gemApvLayout = gemApvPrefGet('layout', 'sequential');
// Source — 'full' (default; locks onto the most recent full-readout
// monitoring event, served from a separate server-side snapshot via
// /api/gem/apv/latest_full) or 'current' (follows the live stream).
// 'full' is the requested RC monitoring view — bypasses online ZS so the
// entire pedestal/noise spectrum is visible across all 128 channels.
let gemApvSource = gemApvPrefGet('source', 'full');
// Normalize — older builds or a corrupted localStorage entry could leave
// these as anything; an unknown gemApvSource value would silently freeze
// the WS gate (both 'event' and 'full_event' branches would return).
if (gemApvLayout !== 'overlay' && gemApvLayout !== 'sequential') gemApvLayout = 'sequential';
if (gemApvSource !== 'current' && gemApvSource !== 'full')       gemApvSource = 'full';
// Pause — freezes auto-refresh on this tab.  WS new_event /
// gem_apv_full_event notifications skip the refetch while true; explicit
// user actions (σ change, navigation, source/layout change) still refresh.
let gemApvPaused = gemApvPrefGet('paused', '0') === '1';
// Staleness counter — bumped on every WS new_event we ignore while
// paused, displayed in the status line so the operator knows how far
// behind the live stream they are.
let gemApvSkippedSincePause = 0;
function gemApvDetVisible(detId) {
    if (detId < 0 || detId >= gemApvDetMask.length) return true;
    return gemApvDetMask[detId];
}
let gemApvBuiltKey = '';        // signature of section layout currently in DOM
const gemApvCanvases = new Map(); // apv_id → canvas element

// Panel size — driven by CSS minmax.  Canvas pixel size matches its
// CSS box at render time so traces stay crisp on HiDPI screens.
const GEM_APV_TITLE_H   = 16;
const GEM_APV_HIT_ROW_H = 6;
// Two stacked hit rows (fw above, sw below) with a 1px gap.  Reserved
// at all times so toggling FW Hits doesn't reflow the plot region.
const GEM_APV_HIT_BLOCK_H = 2 * GEM_APV_HIT_ROW_H + 1;

// =====================================================================
// Fetch + section build
// =====================================================================

// Apply a successfully-fetched APV payload: cache it, sync σ input, rebuild
// section skeleton if needed, refresh calib on rev mismatch, then redraw.
// Shared by per-event fetch and latest-full snapshot fetch so both paths
// stay consistent on calib_rev handling and σ display.
function applyGemApvData(data) {
    gemApvData = data;
    // Pull the event from the payload — the server stamps it and it's the
    // single source of truth for "what is being displayed" (in 'full' mode
    // this is the snapshot's seq, not the live currentEvent).
    gemApvCurrentEvent = (typeof data.event === 'number') ? data.event : -1;
    // Reflect the encode-time σ in the toolbar input so all viewers stay
    // in sync without firing a fresh POST.
    syncGemApvZsSigmaInput(data.zs_sigma);
    buildGemApvSections();
    // Load (or refresh) calibration on rev mismatch — the threshold band
    // needs noise, which lives on /api/gem/calib.  Explicit type check so
    // rev=0 isn't mis-read as "missing".
    const haveRev = (gemApvCalib && typeof gemApvCalib.rev === 'number')
        ? gemApvCalib.rev : null;
    const dataRev = (typeof data.calib_rev === 'number')
        ? data.calib_rev : null;
    if (dataRev !== null && dataRev !== haveRev) {
        ensureGemApvCalib(true /*force*/).then(renderGemApvPanels);
    } else {
        renderGemApvPanels();
    }
}

function fetchGemApvData(evnum) {
    if (typeof evnum !== 'number' || evnum <= 0) return Promise.resolve();
    return fetch(`/api/gem/apv/${evnum}`)
        .then(r => {
            if (!r.ok) throw new Error('http ' + r.status);
            return r.json();
        })
        .then(data => {
            if (data.error) { gemApvSetStatus(data.error); return; }
            applyGemApvData(data);
        })
        .catch(err => gemApvSetStatus('Fetch error: ' + err));
}

// Fetch the server's "latest full-readout" snapshot — the most recent
// monitoring event where firmware ZS was bypassed (so the entire pedestal
// spectrum is visible across all 128 channels per APV).  The server stamps
// the payload with the snapshot's event seq, which applyGemApvData picks up.
function fetchGemApvLatestFull() {
    return fetch('/api/gem/apv/latest_full')
        .then(r => {
            if (!r.ok) throw new Error('http ' + r.status);
            return r.json();
        })
        .then(data => {
            if (data.error) { gemApvSetStatus(data.error); return; }
            applyGemApvData(data);
        })
        .catch(err => gemApvSetStatus('Fetch error: ' + err));
}

// Entry point for the source-aware refresh used by tab-activation,
// source-change, layout-change, and explicit user re-fetch requests.  In
// 'current' mode it follows the live event the rest of the viewer is
// looking at; in 'full' mode it loads the latest server snapshot.  This is
// the bypass path for the pause gate — explicit user actions always run
// through here regardless of gemApvPaused.
function refreshGemApv(currentEventNum) {
    if (gemApvSource === 'full') {
        return fetchGemApvLatestFull();
    }
    return fetchGemApvData(currentEventNum);
}

// Called from online.js's WS handler when the server reports a new event
// (current mode) or a new full-readout event (full mode), and from
// loadEventData on every event swap.  Honours the pause gate when the
// caller signals an auto-refresh; user-initiated navigation (file-mode
// arrow keys / prev-next, online ring-buffer nav) passes manual=true so
// the operator's deliberate "show me this event" still goes through.
// 'full_event' always comes from the WS push and is never manual.
// Pre-update viewers without this function defined fall back to the
// existing direct fetchGemApvData call in viewer.js.
function gemApvOnLiveEvent(evnum, kind /* 'event' | 'full_event' */, manual) {
    if (kind === 'full_event' && gemApvSource !== 'full') return;
    if (kind === 'event'      && gemApvSource !== 'current') return;
    if (gemApvPaused && !manual) {
        gemApvSkippedSincePause++;
        // Cheap status refresh so the operator sees the staleness counter
        // tick even while traces are frozen.
        renderGemApvPanels();
        return;
    }
    if (kind === 'full_event') fetchGemApvLatestFull();
    else                       fetchGemApvData(evnum);
}

// Fetch /api/gem/calib once and cache.  Pass force=true to bypass the
// cache (e.g. on calib_rev mismatch).  Returns a Promise so callers can
// chain a redraw.  Multiple concurrent calls dedupe on the inflight flag.
function ensureGemApvCalib(force) {
    if (!force && gemApvCalib) return Promise.resolve(gemApvCalib);
    if (gemApvCalibInflight) return gemApvCalibInflight;
    const p = fetch('/api/gem/calib')
        .then(r => r.ok ? r.json() : Promise.reject('http ' + r.status))
        .then(j => {
            const noise = new Map();
            for (const a of (j.apvs || [])) {
                if (typeof a.id === 'number' && Array.isArray(a.noise))
                    noise.set(a.id, Float32Array.from(a.noise));
            }
            gemApvCalib = { rev: j.rev || 0, zs_sigma: j.zs_sigma || 0, noise };
            // First-time fetch happens before any event arrives — seed
            // the toolbar input from the calib response.
            syncGemApvZsSigmaInput(j.zs_sigma);
            return gemApvCalib;
        })
        .catch(err => {
            console.warn('gem calib fetch failed:', err);
            return null;
        })
        .finally(() => { gemApvCalibInflight = false; });
    gemApvCalibInflight = p;
    return p;
}

// Update the toolbar σ input without re-firing onchange (which would
// POST back to the server).  Skips if the user is currently editing.
function syncGemApvZsSigmaInput(sigma) {
    const el = document.getElementById('gem-apv-zs-sigma');
    if (!el || sigma == null) return;
    if (document.activeElement === el) return;   // don't clobber typing
    const v = (+sigma).toFixed(1);
    if (el.value !== v) el.value = v;
}

function gemApvSetStatus(text) {
    const el = document.getElementById('gem-apv-stats');
    if (el) el.textContent = text;
}

// Rebuild the section/grid skeleton if the detector layout changed
// (different file, different config).  Cheap when called repeatedly
// with the same data — the skeleton is keyed by det list signature.
function buildGemApvSections() {
    const body = document.getElementById('gem-apv-body');
    if (!body || !gemApvData) return;
    if (!gemApvData.enabled) {
        body.innerHTML = '<div style="padding:40px;text-align:center;color:var(--dim)">GEM not enabled</div>';
        gemApvBuiltKey = '_disabled_';
        return;
    }

    const dets = (gemApvData.detectors || []).slice()
        .sort((a, b) => (a.id - b.id));
    const apvs = gemApvData.apvs || [];

    // Group APVs by det_id, sorted by (plane, crate, mpd, adc).
    const byDet = new Map();
    for (const det of dets) byDet.set(det.id, []);
    for (const apv of apvs) {
        const arr = byDet.get(apv.det_id);
        if (arr) arr.push(apv);
    }
    for (const [, arr] of byDet) {
        arr.sort((a, b) =>
            (a.plane || '').localeCompare(b.plane || '') ||
            (a.crate - b.crate) ||
            (a.mpd   - b.mpd) ||
            (a.adc   - b.adc));
    }

    // Skeleton signature — only rebuild when the per-det APV list changes.
    const key = dets.map(d => `${d.id}:${(byDet.get(d.id) || []).map(a => a.id).join('-')}`).join('|');
    if (key === gemApvBuiltKey) return;
    gemApvBuiltKey = key;

    body.innerHTML = '';
    gemApvCanvases.clear();

    dets.forEach((det) => {
        // Section separators are drawn as a top border on the section
        // itself (see .gem-apv-section in viewer.css) so hiding a GEM
        // via the toolbar checkboxes also hides its separator naturally.
        // The topmost visible section gets .first-visible to suppress its
        // border — applied in renderGemApvPanels after visibility is set.
        const section = document.createElement('div');
        section.className = 'gem-apv-section';
        section.dataset.det = det.id;
        section.style.background = gemApvTint(det.id);

        const header = document.createElement('div');
        header.className = 'gem-apv-section-header';
        header.textContent = `GEM ${det.id} — ${det.name}   (${det.n_apvs} APVs)`;
        section.appendChild(header);

        const grid = document.createElement('div');
        grid.className = 'gem-apv-grid';
        section.appendChild(grid);

        const apvsHere = byDet.get(det.id) || [];
        for (const apv of apvsHere) {
            const panel = document.createElement('div');
            panel.className = 'gem-apv-panel';
            panel.dataset.apvId = apv.id;
            const canvas = document.createElement('canvas');
            canvas.className = 'gem-apv-canvas';
            panel.appendChild(canvas);
            grid.appendChild(panel);
            gemApvCanvases.set(apv.id, canvas);
        }
        body.appendChild(section);
    });
}

// =====================================================================
// Render — called whenever data or controls change
// =====================================================================

function renderGemApvPanels() {
    if (!gemApvData || !gemApvData.enabled) return;
    const apvs  = gemApvData.apvs || [];
    const field = gemApvShowProcessed ? 'processed' : 'raw';

    // Section visibility — toggle whole sections (header + grid + tint
    // background) for GEMs the user has unchecked.  The topmost visible
    // section gets .first-visible so its top border (which acts as the
    // separator above it) is hidden.
    const body = document.getElementById('gem-apv-body');
    let firstVisibleSection = null;
    if (body) {
        body.querySelectorAll('.gem-apv-section').forEach(sec => {
            const detId = parseInt(sec.dataset.det, 10);
            const visible = gemApvDetVisible(detId);
            sec.style.display = visible ? '' : 'none';
            sec.classList.remove('first-visible');
            if (visible && !firstVisibleSection) firstVisibleSection = sec;
        });
        if (firstVisibleSection) firstVisibleSection.classList.add('first-visible');
    }

    // Compute global Y range across visible (non-filtered) APVs.  Skip
    // APVs that didn't show up in this event — their frame is all zeros
    // and would otherwise pull the shared scale toward the origin.
    let yLo = Infinity, yHi = -Infinity;
    if (gemApvSharedY) {
        for (const apv of apvs) {
            if (!gemApvDetVisible(apv.det_id)) continue;
            if (!apv.present) continue;
            if (gemApvShowSignalOnly && !apvHasSignal(apv)) continue;
            const f = apv[field];
            if (!f) continue;
            for (let s = 0; s < f.length; s++) {
                for (let t = 0; t < 6; t++) {
                    if (!gemApvSampleMask[t]) continue;
                    const v = f[s][t];
                    if (v < yLo) yLo = v;
                    if (v > yHi) yHi = v;
                }
            }
        }
    }
    if (!isFinite(yLo) || !isFinite(yHi)) { yLo = 0; yHi = 1; }
    // In Process mode, force the Y axis symmetric around zero so the dashed
    // zero line is always centered and the threshold band reads as a
    // mirrored pair.  Raw mode keeps the auto-fit range — pedestal-level
    // traces (~1500-2500 ADC) all positive, so symmetric would waste half
    // the panel.
    if (gemApvShowProcessed) {
        const m = Math.max(Math.abs(yLo), Math.abs(yHi));
        yLo = -m; yHi = m;
    }
    if (yHi - yLo < 8) {
        const m = 0.5 * (yLo + yHi);
        yLo = m - 4; yHi = m + 4;
    }
    const pad = 0.08 * (yHi - yLo);
    yLo -= pad; yHi += pad;
    const sharedRange = gemApvSharedY ? [yLo, yHi] : null;

    let total = 0, shown = 0;
    for (const apv of apvs) {
        total++;
        const panel = panelOf(apv.id);
        if (!panel) continue;
        // Skip rendering work for APVs in a hidden GEM — the section is
        // already display:none, so any draw would be invisible anyway.
        if (!gemApvDetVisible(apv.det_id)) continue;
        const hasSig = apvHasSignal(apv);
        if (gemApvShowSignalOnly && !hasSig) {
            panel.style.display = 'none';
            continue;
        }
        panel.style.display = '';
        shown++;
        const canvas = gemApvCanvases.get(apv.id);
        if (canvas) drawApvCanvas(canvas, apv, field, sharedRange);
    }

    const mode = gemApvShowProcessed ? 'processed' : 'raw';
    const layoutLbl = (gemApvLayout === 'sequential') ? ' seq' : '';
    const srcLbl = (gemApvSource === 'full') ? ' mon' : '';
    const evlbl = gemApvCurrentEvent > 0 ? `evt ${gemApvCurrentEvent}` : '';
    const pauseLbl = gemApvPaused
        ? `  PAUSED${gemApvSkippedSincePause ? ` (skipped ${gemApvSkippedSincePause})` : ''}`
        : '';
    gemApvSetStatus(`${shown}/${total} APVs  [${mode}${layoutLbl}${srcLbl}]  ${evlbl}${pauseLbl}`);
}

function apvHasSignal(apv) {
    const h = apv.hits;
    if (!h) return false;
    for (let i = 0; i < h.length; i++) if (h[i]) return true;
    return false;
}

function panelOf(apvId) {
    const c = gemApvCanvases.get(apvId);
    return c ? c.parentElement : null;
}

// Draw a single APV: title, zero line, optional threshold band, traces,
// optional CM overlay, hit-tick rows, Y-range labels.  Non-present APVs
// (firmware sent nothing this event) draw as an empty placeholder so a
// missing chip is immediately spottable in the fixed grid.
function drawApvCanvas(canvas, apv, field, sharedRange) {
    // Match the canvas pixel buffer to its CSS box for sharp lines.
    const dpr = window.devicePixelRatio || 1;
    const cssW = canvas.clientWidth || 240;
    const cssH = canvas.clientHeight || 160;
    if (canvas.width  !== Math.round(cssW * dpr) ||
        canvas.height !== Math.round(cssH * dpr)) {
        canvas.width  = Math.round(cssW * dpr);
        canvas.height = Math.round(cssH * dpr);
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const W = cssW, H = cssH;
    // Canvas background: keep slightly darker than the section tint so
    // the panel reads as a "tile" sitting on the tinted GEM row.  Non-
    // present APVs use a desaturated tile so a missing chip stands out
    // against the bright tiles around it.
    const isMissing = !apv.present;
    ctx.fillStyle = isMissing
        ? (THEME && THEME.bgDim ? THEME.bgDim : '#0a0a14')
        : (THEME && THEME.canvas ? THEME.canvas : '#11112a');
    ctx.fillRect(0, 0, W, H);

    // Border — red if firmware reported full readout but no ZS hits,
    // accent if any ZS survivors, dim/dashed if APV didn't show up this
    // event, otherwise neutral.
    let borderCol = THEME && THEME.border ? THEME.border : '#333';
    let borderW = 1;
    let borderDash = null;
    if (isMissing) {
        borderCol = THEME && THEME.textDim ? THEME.textDim : '#666';
        borderDash = [3, 3];
    } else if (apv.no_hit_fr) {
        borderCol = THEME && THEME.danger ? THEME.danger : '#ff6b6b';
    } else if (apvHasSignal(apv) && !gemApvShowSignalOnly) {
        borderCol = THEME && THEME.accent ? THEME.accent : '#ffd166';
        borderW = 2;
    }
    ctx.strokeStyle = borderCol;
    ctx.lineWidth = borderW;
    if (borderDash) ctx.setLineDash(borderDash);
    ctx.strokeRect(0.5, 0.5, W - 1, H - 1);
    if (borderDash) ctx.setLineDash([]);

    // Title row.  Dim the title for missing APVs so the eye groups them
    // with the dashed border / dim background.
    const titleH = GEM_APV_TITLE_H;
    ctx.fillStyle = isMissing
        ? (THEME && THEME.textDim ? THEME.textDim : '#888')
        : (THEME && THEME.text ? THEME.text : '#e0e0e0');
    ctx.font = 'bold 10px ui-monospace, monospace';
    ctx.textBaseline = 'middle';
    const title = `c${apv.crate} m${apv.mpd} a${apv.adc}  ${apv.plane} p${apv.det_pos}`;
    ctx.fillText(title, 4, titleH / 2);
    if (isMissing) {
        ctx.fillStyle = THEME && THEME.textDim ? THEME.textDim : '#888';
        ctx.textAlign = 'right';
        ctx.fillText('no data', W - 4, titleH / 2);
        ctx.textAlign = 'start';
        ctx.textBaseline = 'alphabetic';
        return;   // skip plot/hit-rows — there is nothing to draw
    }
    if (apv.no_hit_fr) {
        ctx.fillStyle = THEME && THEME.danger ? THEME.danger : '#ff6b6b';
        ctx.textAlign = 'right';
        ctx.fillText('no hits', W - 4, titleH / 2);
        ctx.textAlign = 'start';
    }

    // Plot region — reserve two stacked hit rows (fw + sw) at the bottom
    // so toggling FW Hits doesn't reflow the trace area.
    const hitH = GEM_APV_HIT_ROW_H;
    const hitBlockH = GEM_APV_HIT_BLOCK_H;
    const plotX = 4, plotY = titleH + 2;
    const plotW = W - 8, plotH = H - titleH - hitBlockH - 6;
    if (plotW <= 0 || plotH <= 0) return;

    // Y range.
    const frame = apv[field];
    let yLo, yHi;
    if (sharedRange) {
        yLo = sharedRange[0]; yHi = sharedRange[1];
    } else {
        yLo = Infinity; yHi = -Infinity;
        if (frame) {
            for (let s = 0; s < frame.length; s++) {
                for (let t = 0; t < 6; t++) {
                    if (!gemApvSampleMask[t]) continue;
                    const v = frame[s][t];
                    if (v < yLo) yLo = v;
                    if (v > yHi) yHi = v;
                }
            }
        }
        if (!isFinite(yLo) || !isFinite(yHi)) { yLo = 0; yHi = 1; }
        // Process mode: clamp to a symmetric band around zero (see the
        // matching comment in renderGemApvPanels for rationale).
        if (gemApvShowProcessed) {
            const m = Math.max(Math.abs(yLo), Math.abs(yHi));
            yLo = -m; yHi = m;
        }
        if (yHi - yLo < 8) { const m = 0.5*(yLo+yHi); yLo = m-4; yHi = m+4; }
        const pad = 0.08 * (yHi - yLo);
        yLo -= pad; yHi += pad;
    }
    const ySpan = (yHi - yLo) || 1;
    const toY = v => plotY + plotH - (v - yLo) / ySpan * plotH;

    // Zero line if it's in range.
    if (yLo < 0 && yHi > 0) {
        ctx.strokeStyle = THEME && THEME.textDim ? THEME.textDim : '#888';
        ctx.lineWidth = 0.5;
        ctx.setLineDash([2, 2]);
        ctx.beginPath();
        const zy = toY(0);
        ctx.moveTo(plotX, zy); ctx.lineTo(plotX + plotW, zy);
        ctx.stroke();
        ctx.setLineDash([]);
    }

    // X-axis mapper — abstracts overlay vs sequential so threshold band,
    // traces, CM overlay and hit ticks all stay aligned.
    //   overlay:    x depends only on channel (all 6 TS share an axis)
    //   sequential: x = (ts * nStrips + ch) — channel-major within each TS
    //               block, TS-major across blocks (prad1 style)
    const N_TS = 6;
    const nStrips = 128;
    const seq = (gemApvLayout === 'sequential');
    const xSlotsOverlay    = nStrips;
    const xSlotsSequential = nStrips * N_TS;
    const xSlots = seq ? xSlotsSequential : xSlotsOverlay;
    const stepX  = plotW / Math.max(xSlots - 1, 1);
    const xAt = (ch, ts) => plotX + (seq ? (ts * nStrips + ch) : ch) * stepX;
    // List of TS blocks the per-TS overlays (threshold band, hit ticks)
    // should iterate over — only 6 in sequential layout, just [0] in overlay
    // (where one band/hit-row covers the shared axis).
    const tsBlocks = seq ? [0,1,2,3,4,5] : [0];

    // Threshold band: ±noise[ch]·zs_sigma, dashed grey, drawn before the
    // data traces so traces sit on top.  noise[] comes from the cached
    // /api/gem/calib payload; zs_sigma is per-event so the band always
    // tracks the σ used to encode this event's hits[].  Only meaningful
    // for processed view (raw view shows pre-pedestal-subtraction).
    const zsSigma = (gemApvData && gemApvData.zs_sigma) || 0;
    const noise = (gemApvCalib && gemApvCalib.noise) ? gemApvCalib.noise.get(apv.id) : null;
    if (gemApvShowThreshold && gemApvShowProcessed && noise && zsSigma > 0) {
        const noiseN = Math.min(nStrips, noise.length);
        ctx.strokeStyle = THEME && THEME.textDim ? THEME.textDim : '#888';
        ctx.lineWidth = 0.8;
        ctx.setLineDash([4, 3]);
        const drawBand = (sign, t) => {
            ctx.beginPath();
            for (let s = 0; s < noiseN; s++) {
                const x = xAt(s, t);
                const y = toY(sign * noise[s] * zsSigma);
                if (s === 0) ctx.moveTo(x, y);
                else         ctx.lineTo(x, y);
            }
            ctx.stroke();
        };
        for (const t of tsBlocks) { drawBand(+1, t); drawBand(-1, t); }
        ctx.setLineDash([]);
    }

    // Time-sample traces.
    if (frame && frame.length > 0) {
        const nS = Math.min(nStrips, frame.length);
        ctx.lineWidth = 0.9;
        for (let t = 0; t < N_TS; t++) {
            if (!gemApvSampleMask[t]) continue;
            ctx.strokeStyle = GEM_APV_TS_COLORS[t];
            ctx.beginPath();
            for (let s = 0; s < nS; s++) {
                const x = xAt(s, t);
                const y = toY(frame[s][t]);
                if (s === 0) ctx.moveTo(x, y);
                else         ctx.lineTo(x, y);
            }
            ctx.stroke();
        }
    }

    // Sequential dividers — thin dashed vertical lines between TS blocks.
    // Drawn after traces so they're visible on top, but kept dim so they
    // read as separators rather than data.
    if (seq) {
        ctx.strokeStyle = THEME && THEME.textDim ? THEME.textDim : '#888';
        ctx.lineWidth = 0.4;
        ctx.setLineDash([1, 3]);
        for (let t = 1; t < N_TS; t++) {
            const x = xAt(0, t);
            ctx.beginPath();
            ctx.moveTo(x, plotY);
            ctx.lineTo(x, plotY + plotH);
            ctx.stroke();
        }
        ctx.setLineDash([]);
    }

    // CM overlay: one horizontal dashed line per enabled time sample,
    // colour-matched (desaturated) with the trace so reader can pair
    // firmware CM with the same-colour strip waveform.  Drawn AFTER the
    // traces so it sits on top.  Skipped in Process mode (raw ADC values
    // would land off the pedestal-subtracted axis) and when the firmware
    // didn't emit type-0xD debug-header words (apv.cm == null).  In
    // sequential layout each line is constrained to its TS block so the
    // CM value lines up under the matching colour-coded trace block.
    if (gemApvShowCm && !gemApvShowProcessed && Array.isArray(apv.cm)) {
        ctx.lineWidth = 1.4;
        ctx.setLineDash([5, 3]);
        for (let t = 0; t < apv.cm.length && t < N_TS; t++) {
            if (!gemApvSampleMask[t]) continue;
            ctx.strokeStyle = GEM_APV_CM_COLORS[t];
            const y = toY(apv.cm[t]);
            ctx.beginPath();
            if (seq) {
                ctx.moveTo(xAt(0, t), y);
                ctx.lineTo(xAt(nStrips - 1, t), y);
            } else {
                ctx.moveTo(plotX, y);
                ctx.lineTo(plotX + plotW, y);
            }
            ctx.stroke();
        }
        ctx.setLineDash([]);
    }

    // Hit tick rows — bottom row = software-cut survivors (bright accent),
    // top row = firmware survivors (dim, gated by FW Hits checkbox).
    // Both rows reserved at all times so toggling FW Hits doesn't reflow.
    // hits[] / fw_hits[] are per-channel only (no TS dimension); in
    // sequential layout we repeat the same per-channel tick pattern under
    // each TS block so traces and ticks stay column-aligned.
    const swRowY = H - hitH - 2;
    const fwRowY = swRowY - hitH - 1;
    const accent = (THEME && THEME.accent) ? THEME.accent : '#ffd166';
    if (gemApvShowFwHits && Array.isArray(apv.fw_hits) && apv.fw_hits.length > 0) {
        const nS = Math.min(nStrips, apv.fw_hits.length);
        ctx.globalAlpha = 0.45;
        ctx.fillStyle = accent;
        for (const t of tsBlocks) {
            for (let s = 0; s < nS; s++) {
                if (apv.fw_hits[s]) {
                    const x = xAt(s, t);
                    ctx.fillRect(x - 0.8, fwRowY, 1.6, hitH);
                }
            }
        }
        ctx.globalAlpha = 1.0;
    }
    if (Array.isArray(apv.hits) && apv.hits.length > 0) {
        const nS = Math.min(nStrips, apv.hits.length);
        ctx.fillStyle = accent;
        for (const t of tsBlocks) {
            for (let s = 0; s < nS; s++) {
                if (apv.hits[s]) {
                    const x = xAt(s, t);
                    ctx.fillRect(x - 0.8, swRowY, 1.6, hitH);
                }
            }
        }
    }

    // Tiny Y-range labels in the plot corners — useful when shared Y is
    // off so each panel's auto-scale is visible at a glance.
    ctx.fillStyle = THEME && THEME.textDim ? THEME.textDim : '#888';
    ctx.font = '8px ui-monospace, monospace';
    ctx.textBaseline = 'top';
    ctx.fillText(fmtCompact(yHi), plotX + 2, plotY + 1);
    ctx.textBaseline = 'bottom';
    ctx.fillText(fmtCompact(yLo), plotX + 2, plotY + plotH - 1);
    ctx.textBaseline = 'alphabetic';
}

function fmtCompact(v) {
    if (Math.abs(v) < 1000) return v.toFixed(0);
    return (v / 1000).toFixed(1) + 'k';
}

// =====================================================================
// Controls
// =====================================================================

function setupGemApvControls() {
    const cb = (id, on) => {
        const el = document.getElementById(id);
        if (!el) return;
        el.checked = on;
        el.onchange = () => {
            switch (id) {
                case 'gem-apv-process':
                    gemApvShowProcessed = el.checked;
                    syncGemApvControlEnables();
                    break;
                case 'gem-apv-signal-only':
                    gemApvShowSignalOnly = el.checked; break;
                case 'gem-apv-shared-y':
                    gemApvSharedY = el.checked; break;
                case 'gem-apv-threshold':
                    gemApvShowThreshold = el.checked; break;
                case 'gem-apv-fw-hits':
                    gemApvShowFwHits = el.checked; break;
                case 'gem-apv-cm':
                    gemApvShowCm = el.checked; break;
            }
            renderGemApvPanels();
        };
    };
    cb('gem-apv-process',     gemApvShowProcessed);
    cb('gem-apv-signal-only', gemApvShowSignalOnly);
    cb('gem-apv-shared-y',    gemApvSharedY);
    cb('gem-apv-fw-hits',     gemApvShowFwHits);
    cb('gem-apv-threshold',   gemApvShowThreshold);
    cb('gem-apv-cm',          gemApvShowCm);
    syncGemApvControlEnables();

    // Layout select — overlay (default) vs sequential (prad1 style).
    // Pure render-side toggle, no refetch needed.
    const layoutEl = document.getElementById('gem-apv-layout');
    if (layoutEl) {
        layoutEl.value = gemApvLayout;
        layoutEl.onchange = () => {
            gemApvLayout = layoutEl.value;
            gemApvPrefSet('layout', gemApvLayout);
            renderGemApvPanels();
        };
    }

    // Source select — 'current' follows the live stream; 'full' locks the
    // panel onto the latest server-side full-readout snapshot.  Changing
    // source is an explicit user action so it bypasses the pause gate and
    // also resets the staleness counter (the "skipped N" only applies to
    // the prior source's live feed).
    const srcEl = document.getElementById('gem-apv-source');
    if (srcEl) {
        srcEl.value = gemApvSource;
        srcEl.onchange = () => {
            gemApvSource = srcEl.value;
            gemApvPrefSet('source', gemApvSource);
            gemApvSkippedSincePause = 0;
            // Pull a fresh frame from the new source so the panel updates
            // immediately instead of waiting for the next WS notification.
            // currentEvent is exposed by viewer.js as a global.
            const live = (typeof currentEvent === 'number') ? currentEvent : -1;
            refreshGemApv(live);
        };
    }

    // Pause button — toggles auto-refresh on this tab.  All other tabs
    // continue updating normally; explicit user actions on this tab also
    // bypass the pause via refreshGemApv.
    const pauseEl = document.getElementById('gem-apv-pause');
    if (pauseEl) {
        const syncPauseUi = () => {
            pauseEl.textContent = gemApvPaused ? '⏸ Paused' : '▶ Live';
            pauseEl.classList.toggle('gem-apv-paused', gemApvPaused);
        };
        syncPauseUi();
        pauseEl.onclick = () => {
            gemApvPaused = !gemApvPaused;
            gemApvPrefSet('paused', gemApvPaused ? '1' : '0');
            if (!gemApvPaused) {
                // Resume — pull the freshest frame for the active source
                // so the panel jumps straight to "now" instead of holding
                // the stale frame until the next WS notification.
                gemApvSkippedSincePause = 0;
                const live = (typeof currentEvent === 'number') ? currentEvent : -1;
                refreshGemApv(live);
            } else {
                // Pausing — just refresh the status line so the "PAUSED"
                // label appears immediately.
                renderGemApvPanels();
            }
            syncPauseUi();
        };
    }

    // σ input — POST on commit (change event fires on blur / Enter / arrows).
    // The server applies to the live reconstruction; new events arriving
    // through the WS will reflect the change in their per-event zs_sigma,
    // which then rolls back into syncGemApvZsSigmaInput for everyone.
    const zsEl = document.getElementById('gem-apv-zs-sigma');
    if (zsEl) {
        zsEl.onchange = () => {
            const v = parseFloat(zsEl.value);
            if (!isFinite(v) || v < 0) {
                // Restore the last known good value rather than POST garbage.
                if (gemApvCalib && gemApvCalib.zs_sigma != null)
                    zsEl.value = (+gemApvCalib.zs_sigma).toFixed(1);
                return;
            }
            postGemApvZsSigma(v);
        };
    }
    // Pre-fetch calibration so the threshold band is ready by the time
    // the first event renders (the band wouldn't draw without noise[]).
    ensureGemApvCalib(false);
    // Per-GEM filter (gem0…gem3) — hides whole sections, including the
    // separator above (which is a top border on the section itself).
    for (let d = 0; d < gemApvDetMask.length; d++) {
        const el = document.getElementById('gem-apv-d' + d);
        if (!el) continue;
        el.checked = gemApvDetMask[d];
        el.onchange = () => {
            gemApvDetMask[d] = el.checked;
            renderGemApvPanels();
        };
    }
    for (let t = 0; t < 6; t++) {
        const el = document.getElementById('gem-apv-t' + t);
        if (!el) continue;
        el.checked = gemApvSampleMask[t];
        el.onchange = () => {
            gemApvSampleMask[t] = el.checked;
            renderGemApvPanels();
        };
    }
}

// The two overlays have opposite Y-axis requirements:
//   Threshold band — needs a pedestal-subtracted axis (Process on)
//   CM overlay     — values are raw ADC counts, only fit the raw axis
// so we mirror that in the UI by disabling whichever checkbox is
// inapplicable in the current mode.
function syncGemApvControlEnables() {
    const thr = document.getElementById('gem-apv-threshold');
    if (thr) thr.disabled = !gemApvShowProcessed;
    const cm  = document.getElementById('gem-apv-cm');
    if (cm)  cm.disabled  =  gemApvShowProcessed;
}

// POST a new σ to the server.  Cached calib is updated optimistically
// so the next render uses the new value; a fresh fetch isn't needed
// (noise[] hasn't changed).  Other viewers pick the change up through
// the per-event zs_sigma echoed by the server in subsequent events.
function postGemApvZsSigma(v) {
    fetch('/api/gem/threshold', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ zs_sigma: v }),
    })
        .then(r => r.ok ? r.json() : Promise.reject('http ' + r.status))
        .then(j => {
            if (j.error) throw new Error(j.error);
            if (gemApvCalib) gemApvCalib.zs_sigma = j.zs_sigma;
            // Re-fetch from the active source so the hits[] (and the
            // threshold band's per-event zs_sigma) reflect the new σ.  In
            // 'current' mode this hits /api/gem/apv/<displayed-event> so
            // a paused panel keeps its frame and just rolls in the new σ
            // instead of jumping to the live stream; in 'full' mode
            // refreshGemApv ignores the arg and pulls the snapshot.
            // σ changes are explicit user actions so they bypass pause.
            refreshGemApv(gemApvCurrentEvent);
        })
        .catch(err => {
            console.warn('gem threshold POST failed:', err);
            // Rollback the input to the server's last known good value.
            if (gemApvCalib && gemApvCalib.zs_sigma != null)
                syncGemApvZsSigmaInput(gemApvCalib.zs_sigma);
        });
}

// Resize: just redraw — CSS auto-grid handles re-flow, canvas redraw
// adapts to new clientWidth/clientHeight on each render call.
function resizeGemApv() {
    if (!gemApvData || !gemApvData.enabled) return;
    renderGemApvPanels();
}

// Theme flip — every per-APV canvas reads THEME at draw time for the
// background tile, border, title, zero/threshold lines, and Y-range labels.
// Replay the full render so all of those pick up the new palette.
if (typeof onThemeChange === 'function') {
    onThemeChange(() => {
        if (gemApvData && gemApvData.enabled) renderGemApvPanels();
    });
}
