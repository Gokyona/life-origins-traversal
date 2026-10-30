/* Gokyona Environment - browser UI. Plain JS + canvas, no libraries (works offline).
   Every render function must tolerate missing / null fields: show "—", never throw. */
"use strict";

const POLL_MS = 2000, CHECK_MS = 6000, HW_MS = 30000;
const FALLBACK_BANNER = "Model, not life. Everything shown here comes from a simulation model, not from living systems. Length beyond 64 is unvalidated.";
const FALLBACK_TIME_NOTE = "Time scale follows the model's rate constants; it is not a calibrated Earth clock.";
const LEVELS = ["E1", "E2", "L1", "L2", "E3", "L3_loose", "L3_strict", "L4"];
const COL = { c1: "#5ab0ff", c2: "#f2b63c", c3: "#b58cff", c4: "#3ecf6e", red: "#ef5b5b", grid: "#252b33", axis: "#8b95a3", ctl: "#b58cff" };

const S = {
  params: [], presets: {}, run: "", runs: [],
  status: null, hist: [], histN: 0, perf: [],
  envVals: {}, touched: {}, check: null, hardware: {}, busy: false, token: 0,
};

/* ---------------------------------------------------------------- utilities */
const $ = (id) => document.getElementById(id);
function get(o, path) {
  let v = o;
  for (const k of path.split(".")) { if (v == null || typeof v !== "object") return undefined; v = v[k]; }
  return v;
}
const isNum = (v) => typeof v === "number" && isFinite(v);
function fmt(v, d = 2) {
  if (v === null || v === undefined || (typeof v === "number" && !isFinite(v))) return "—";
  if (typeof v === "boolean") return v ? "yes" : "no";
  if (typeof v !== "number") return String(v);
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(2);
  if (Number.isInteger(v)) return v.toLocaleString("en-US");
  return Number(v.toFixed(d)).toLocaleString("en-US", { maximumFractionDigits: d });
}
function fmtPhys(v) {
  if (!isNum(v)) return "—";
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-2)) return v.toExponential(2);
  return Number(v.toPrecision(3)).toString();
}
function fmtPct(v) { return isNum(v) ? (v * 100).toFixed(1) + " %" : "—"; }
function fmtDur(s) {
  if (!isNum(s)) return "—";
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60);
  return h > 0 ? `${h} h ${m} min` : `${m} min ${Math.floor(s % 60)} s`;
}
function kv(el, pairs) {
  el.textContent = "";
  for (const [k, v] of pairs) {
    const row = document.createElement("div");
    const a = document.createElement("span"); a.className = "k"; a.textContent = k;
    const b = document.createElement("span"); b.className = "v"; b.textContent = (v === undefined || v === null || v === "") ? "—" : v;
    row.append(a, b); el.append(row);
  }
}
function safe(name, fn) {
  try { fn(); } catch (e) { console.error("render " + name, e); }
}
async function api(path, body) {
  const opt = body === undefined ? {} : { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  const r = await fetch(path, opt);
  let j = null;
  try { j = await r.json(); } catch (e) { j = null; }
  if (!r.ok) throw new Error((j && j.error) || ("HTTP " + r.status));
  return j;
}
function msg(text, kind) { const el = $("msg"); el.textContent = text || ""; el.className = "msg " + (kind || ""); }
function xOf(h) { return isNum(h.sim_h) ? h.sim_h : (isNum(h.wall_s) ? h.wall_s : null); }
function series(key) {
  const out = [];
  for (const h of S.hist) { if (!h) continue; const x = xOf(h); out.push([x, isNum(h[key]) ? h[key] : null]); }
  return out;
}
function xLabel() {
  const h = S.hist.find((r) => r);
  return h && !isNum(h.sim_h) && isNum(h.wall_s) ? "wall s" : "sim h";
}

/* ---------------------------------------------------------------- canvas charts */
function prepCanvas(cv) {
  const dpr = window.devicePixelRatio || 1;
  if (!cv.dataset.h) cv.dataset.h = cv.getAttribute("height") || "150";
  const w = Math.max(50, cv.clientWidth), h = +cv.dataset.h;
  cv.style.height = h + "px";
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  }
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  ctx.font = "11px system-ui, sans-serif";
  return { ctx, w, h };
}
function niceTicks(lo, hi, n = 4) {
  if (!(hi > lo)) return [lo];
  const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => s >= raw) || raw;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
  return out;
}
function tickLabel(v) {
  const a = Math.abs(v);
  if (a >= 1e6 || (a > 0 && a < 1e-3)) return v.toExponential(0);
  if (a >= 1e4) return (v / 1000).toFixed(0) + "k";
  if (a >= 1000) return (v / 1000).toFixed(1) + "k";
  return String(+v.toPrecision(3));
}
function noData(ctx, w, h, text) {
  ctx.fillStyle = COL.axis; ctx.textAlign = "center"; ctx.fillText(text || "no data yet", w / 2, h / 2);
}
/* opts: {series:[{pts,color,dash,width}], logY, markers:[{x,color,filled}], hlines:[{y,color,dash}], xlabel, yMin, yMax} */
function lineChart(cv, opts) {
  const { ctx, w, h } = prepCanvas(cv);
  const L = 46, R = 10, T = 8, B = 22, pw = w - L - R, ph = h - T - B;
  const logY = !!opts.logY;
  const ty = (v) => (logY ? (v > 0 ? Math.log10(v) : null) : v);
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  const ser = (opts.series || []).map((s) => {
    let pts = (s.pts || []).filter((p) => p && isNum(p[0]));
    if (pts.length > 1500) { const st = Math.ceil(pts.length / 1500); pts = pts.filter((_, i) => i % st === 0 || i === pts.length - 1); }
    return Object.assign({}, s, { pts });
  });
  for (const s of ser) for (const [x, y] of s.pts) {
    xmin = Math.min(xmin, x); xmax = Math.max(xmax, x);
    const t = isNum(y) ? ty(y) : null;
    if (t !== null && isNum(t)) { ymin = Math.min(ymin, t); ymax = Math.max(ymax, t); }
  }
  for (const hl of opts.hlines || []) { const t = isNum(hl.y) ? ty(hl.y) : null; if (t !== null && isNum(t)) { ymin = Math.min(ymin, t); ymax = Math.max(ymax, t); } }
  for (const m of opts.markers || []) if (isNum(m.x)) { xmin = Math.min(xmin, m.x); xmax = Math.max(xmax, m.x); }
  if (!isFinite(xmin) || !isFinite(ymin)) { noData(ctx, w, h); return; }
  if (isNum(opts.yMin)) ymin = Math.min(ymin, ty(opts.yMin));
  if (isNum(opts.yMax)) ymax = Math.max(ymax, ty(opts.yMax));
  if (xmax === xmin) { xmax = xmin + 1; }
  if (ymax === ymin) { ymax += logY ? 0.5 : (Math.abs(ymax) * 0.1 || 1); ymin -= logY ? 0.5 : (Math.abs(ymin) * 0.1 || 1); }
  else if (!logY) { const pad = (ymax - ymin) * 0.06; ymax += pad; if (ymin !== 0) ymin -= pad; }
  const X = (x) => L + ((x - xmin) / (xmax - xmin)) * pw;
  const Y = (t) => T + ph - ((t - ymin) / (ymax - ymin)) * ph;

  // grid + axes
  ctx.strokeStyle = COL.grid; ctx.fillStyle = COL.axis; ctx.lineWidth = 1;
  const yt = logY ? (() => { const a = []; for (let e = Math.ceil(ymin); e <= Math.floor(ymax); e++) a.push(e); return a.length ? a : niceTicks(ymin, ymax); })() : niceTicks(ymin, ymax);
  ctx.textAlign = "right"; ctx.textBaseline = "middle";
  for (const t of yt) {
    const y = Y(t); ctx.beginPath(); ctx.moveTo(L, y); ctx.lineTo(L + pw, y); ctx.stroke();
    ctx.fillText(logY ? tickLabel(Math.pow(10, t)) : tickLabel(t), L - 4, y);
  }
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of niceTicks(xmin, xmax, Math.max(2, Math.floor(pw / 90)))) {
    const x = X(t); ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.fillText(tickLabel(t), x, T + ph + 4);
  }
  if (opts.xlabel) { ctx.textAlign = "right"; ctx.textBaseline = "top"; ctx.fillStyle = COL.axis; ctx.fillText(opts.xlabel, L + pw - 4, T + 3); }
  ctx.strokeStyle = "#3a414b"; ctx.strokeRect(L, T, pw, ph);

  ctx.save(); ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip();
  for (const hl of opts.hlines || []) {
    const t = isNum(hl.y) ? ty(hl.y) : null; if (t === null || !isNum(t)) continue;
    ctx.strokeStyle = hl.color || COL.axis; ctx.setLineDash(hl.dash ? [5, 4] : []); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(L, Y(t)); ctx.lineTo(L + pw, Y(t)); ctx.stroke();
  }
  for (const m of opts.markers || []) {
    if (!isNum(m.x)) continue;
    const x = X(m.x); ctx.strokeStyle = m.color || COL.c2; ctx.globalAlpha = 0.5; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(x, T); ctx.lineTo(x, T + ph); ctx.stroke(); ctx.globalAlpha = 1; ctx.setLineDash([]);
  }
  for (const s of ser) {
    ctx.strokeStyle = s.color || COL.c1; ctx.lineWidth = s.width || 1.6; ctx.setLineDash(s.dash ? [6, 4] : []);
    ctx.beginPath(); let pen = false;
    for (const [x, y] of s.pts) {
      const t = isNum(y) ? ty(y) : null;
      if (t === null || !isNum(t)) { pen = false; continue; }
      if (!pen) { ctx.moveTo(X(x), Y(t)); pen = true; } else ctx.lineTo(X(x), Y(t));
    }
    ctx.stroke();
    if (s.pts.length === 1 && isNum(s.pts[0][1])) { ctx.fillStyle = s.color; ctx.beginPath(); ctx.arc(X(s.pts[0][0]), Y(ty(s.pts[0][1])), 2.5, 0, 7); ctx.fill(); }
  }
  ctx.restore(); ctx.setLineDash([]);
  // sweep triangles on the time axis
  for (const m of opts.markers || []) {
    if (!isNum(m.x)) continue;
    const x = X(m.x), y = T + ph;
    ctx.beginPath(); ctx.moveTo(x, y - 9); ctx.lineTo(x - 5, y); ctx.lineTo(x + 5, y); ctx.closePath();
    ctx.strokeStyle = m.color || COL.c2; ctx.fillStyle = m.color || COL.c2; ctx.lineWidth = 1.5;
    if (m.filled) ctx.fill(); else ctx.stroke();
  }
}
function legend(el, items) {
  el.textContent = "";
  for (const it of items) {
    const s = document.createElement("span"); s.style.color = it.color;
    const i = document.createElement("i"); i.style.background = it.dash ? "" : it.color; if (it.dash) i.className = "dash";
    const t = document.createElement("span"); t.style.color = "var(--muted)"; t.textContent = it.label;
    s.append(i, t); el.append(s);
  }
}

/* ---------------------------------------------------------------- run control */
function sliderToVal(p, pos) {
  if (p.log) { const a = Math.log10(p.min), b = Math.log10(p.max); return Math.pow(10, a + (pos / 1000) * (b - a)); }
  return +pos;
}
function valToSlider(p, v) {
  if (!isNum(v)) v = p.log ? Math.sqrt(p.min * p.max) : (p.min + p.max) / 2;
  v = Math.min(p.max, Math.max(p.min, v));
  if (p.log) { const a = Math.log10(p.min), b = Math.log10(p.max); return Math.round(((Math.log10(v) - a) / (b - a)) * 1000); }
  return v;
}
function buildSliders() {
  const box = $("sliders"); box.textContent = "";
  for (const p of S.params) {
    if (!p || !p.name || !isNum(p.min) || !isNum(p.max)) continue;
    const row = document.createElement("div"); row.className = "slider";
    const lbl = document.createElement("span"); lbl.className = "lbl"; lbl.textContent = p.label || p.name; lbl.title = `${p.label || p.name} (${p.name})${p.log ? " — log scale" : ""}`;
    const inp = document.createElement("input"); inp.type = "range"; inp.id = "env-" + p.name;
    if (p.log) { inp.min = 0; inp.max = 1000; inp.step = 1; }
    else { inp.min = p.min; inp.max = p.max; inp.step = p.step || (p.max - p.min) / 1000; }
    const sv = document.createElement("span"); sv.className = "sv";
    inp.addEventListener("input", () => {
      S.envVals[p.name] = sliderToVal(p, +inp.value); S.touched[p.name] = true; updateSliderLabel(p);
    });
    row.append(lbl, inp, sv); box.append(row);
  }
}
function updateSliderLabel(p) {
  const inp = $("env-" + p.name); if (!inp) return;
  const row = inp.parentElement, v = S.envVals[p.name];
  row.querySelector(".sv").textContent = fmtPhys(v) + (p.unit ? " " + p.unit : "") + (p.log ? " ㏒" : "");
  row.querySelector(".lbl").classList.toggle("changed", !!S.touched[p.name]);
  $("env-changed").textContent = Object.keys(S.touched).length;
}
function applyPreset() {
  const pr = S.presets[$("preset").value] || {};
  const env = pr.env || {};
  S.envVals = {}; S.touched = {};
  for (const p of S.params) {
    if (!p || !p.name) continue;
    S.envVals[p.name] = isNum(env[p.name]) ? env[p.name] : null;
    const inp = $("env-" + p.name);
    if (!inp) continue;
    inp.value = valToSlider(p, env[p.name]);
    if (!isNum(S.envVals[p.name])) S.envVals[p.name] = sliderToVal(p, +inp.value);
    updateSliderLabel(p);
  }
  const t = pr.trait || {};
  $("t-enabled").checked = !!t.enabled;
  $("t-mode").value = ["self", "shared", "membrane", "stability", "division", "retention"].includes(t.mode) ? t.mode : "retention";
  $("t-pair").checked = t.pair_score !== false;
  $("t-s").value = isNum(t.s) ? t.s : 0.5; $("t-s-val").textContent = (+$("t-s").value).toFixed(2);
  $("t-lrep").value = isNum(t.L_rep) ? t.L_rep : 24;
  $("t-land").value = t.landscape === "additive-mean" ? "additive-mean" : "domain";
  $("t-epi").checked = !!t.epistatic;
  $("t-k").value = 6;
  $("t-seed").value = isNum(t.seed) ? t.seed : 2026;
  $("t-fid").checked = !!t.fidelity_link;
  $("t-sf").value = isNum(t.s_f) ? t.s_f : 0;
  syncTraitEnabled();
}
function syncTraitEnabled() {
  const on = $("t-enabled").checked;
  for (const id of ["t-mode", "t-s", "t-lrep", "t-seed", "t-fid", "t-land", "t-epi", "t-pair"]) $(id).disabled = !on;
  $("t-sf").disabled = !on || !$("t-fid").checked;
}
function collectCreate() {
  const env = {};
  for (const k of Object.keys(S.touched)) if (isNum(S.envVals[k])) env[k] = S.envVals[k];
  return {
    name: $("run-name").value.trim(), preset: $("preset").value, env,
    trait: {
      enabled: $("t-enabled").checked, mode: $("t-mode").value, s: +$("t-s").value,
      L_rep: parseInt($("t-lrep").value, 10), seed: parseInt($("t-seed").value, 10),
      fidelity_link: $("t-fid").checked, s_f: +$("t-sf").value || 0,
      landscape: $("t-land").value, epistatic: $("t-epi").checked, pair_score: $("t-pair").checked,
    },
    multi_start: $("multi-start").checked, ration: $("ration").value, census_years: +$("census-y").value || 100,
    lmax: parseInt($("lmax").value, 10), sites: parseInt($("sites").value, 10), slots: parseInt($("slots").value, 10),
    vesicle_slots: parseInt($("vslots").value, 10),
    migration_m: +$("mig-m").value || 0, migration_mode: $("mig-mode").value, half_pool_mix: $("half-mix").checked,
  };
}
async function doCreate() {
  const body = collectCreate();
  if (!/^[A-Za-z0-9_-]{1,40}$/.test(body.name)) { msg("Run name: use only letters, digits, _ and -, at most 40 characters.", "err"); return; }
  for (const k of ["lmax", "sites", "slots"]) if (!Number.isInteger(body[k])) { msg(k + " must be a whole number.", "err"); return; }
  try {
    await api("/api/create", body);
    msg(`Created run "${body.name}" (config.json written; not started).`, "ok");
    await refreshRuns(); selectRun(body.name);
  } catch (e) { msg("Create failed: " + e.message, "err"); }
}
async function doStart() {
  const run = S.run;
  if (!run) { msg("Select or create a run first.", "err"); return; }
  try {
    const c = await api("/api/check?run=" + encodeURIComponent(run));
    S.check = c; renderCheck();
    if (!c.ok) { msg("Start refused: " + c.reason, "err"); return; }
    if (!confirm(`Start the simulation engine for run "${run}" on this machine's GPU now?`)) return;
    msg("Starting…");
    const r = await api("/api/start", { run });
    msg(`Started run "${run}" (pid ${r.pid}). Log: ${r.log}`, "ok");
  } catch (e) { msg("Start failed: " + e.message, "err"); }
  checkNow();
}
async function doControl(cmd) {
  const run = S.run;
  if (!run) { msg("Select a run first.", "err"); return; }
  if (cmd === "stop" && !confirm(`Stop run "${run}"? The engine is asked to save and exit.`)) return;
  try { await api("/api/control", { run, cmd }); msg(`Sent "${cmd}" to run "${run}". The engine picks it up at its next report.`, "ok"); }
  catch (e) { msg(cmd + " failed: " + e.message, "err"); }
}
function renderCheck() {
  const el = $("check-line"), c = S.check;
  if (!c) { el.textContent = "Start check: —"; el.className = "check"; }
  else { el.textContent = "Start check: " + (c.reason || (c.ok ? "allowed" : "refused")); el.className = "check " + (c.ok ? "ok" : "no"); }
  const st = get(S.status, "state");
  const alive = st === "running" || st === "paused" || st === "starting";
  $("btn-start").disabled = !S.run || !c || !c.ok || alive;
  $("btn-pause").disabled = !S.run || st !== "running";
  $("btn-resume").disabled = !S.run || st !== "paused";
  $("btn-stop").disabled = !S.run || !alive;
}
async function checkNow() {
  try { S.check = await api("/api/check" + (S.run ? "?run=" + encodeURIComponent(S.run) : "")); }
  catch (e) { S.check = { ok: false, reason: "check failed: " + e.message }; }
  renderCheck();
}

/* ---------------------------------------------------------------- runs + polling */
async function refreshRuns() {
  const runs = await api("/api/runs");
  S.runs = Array.isArray(runs) ? runs : [];
  const sel = $("run-select"), cur = S.run;
  const sig = S.runs.map((r) => r.name + ":" + r.state).join("|");
  if (sel.dataset.sig === sig) return;
  sel.dataset.sig = sig; sel.textContent = "";
  if (!S.runs.length) { const o = document.createElement("option"); o.value = ""; o.textContent = "(no runs)"; sel.append(o); }
  for (const r of S.runs) {
    const o = document.createElement("option"); o.value = r.name; o.textContent = `${r.name} — ${r.state || "?"}`; sel.append(o);
  }
  if (cur && S.runs.some((r) => r.name === cur)) sel.value = cur;
  else if (S.runs.length) {
    const best = S.runs.slice().sort((a, b) => (b.mtime || 0) - (a.mtime || 0))[0];
    selectRun(best.name);
  } else selectRun("");
}
function selectRun(name) {
  if (name === S.run) { $("run-select").value = name; return; }
  S.run = name; S.token++; S.status = null; S.hist = []; S.histN = 0; S.perf = [];
  $("run-select").value = name;
  try { localStorage.setItem("gk-run", name); } catch (e) { /* ignore */ }
  renderAll(); checkNow();
}
async function poll() {
  if (S.busy) return;
  S.busy = true;
  try {
    await refreshRuns();
    if (S.run) {
      const tok = S.token, run = encodeURIComponent(S.run);
      let st = null;
      try { st = await api("/api/status?run=" + run); } catch (e) { st = null; }
      const hist = await api(`/api/history?run=${run}&since=${S.histN}`);
      if (tok === S.token) {
        S.status = st && typeof st === "object" ? st : null;
        if (Array.isArray(hist) && hist.length) {
          S.histN += hist.length; S.hist.push(...hist);
          if (S.hist.length > 20000) S.hist = S.hist.slice(-20000);
        }
        const p = get(S.status, "perf"), wall = get(S.status, "time.wall_s");
        if (p && isNum(p.step_ms) && isNum(wall) && (!S.perf.length || S.perf[S.perf.length - 1][0] !== wall)) {
          S.perf.push([wall, p.step_ms]); if (S.perf.length > 3000) S.perf.shift();
        }
      }
    }
    $("conn").textContent = "updated " + new Date().toLocaleTimeString();
    $("conn").className = "conn";
  } catch (e) {
    $("conn").textContent = "server not reachable (" + e.message + ")"; $("conn").className = "conn bad";
  } finally {
    S.busy = false;
  }
  renderAll();
}
async function pollHardware() {
  try { const h = await api("/api/hardware"); S.hardware = h && typeof h === "object" ? h : {}; } catch (e) { S.hardware = {}; }
  safe("hw", renderHeader);
}

/* ---------------------------------------------------------------- panel renders */
function renderAll() {
  safe("header", renderHeader); safe("time", renderTime); safe("chains", renderChains); safe("lineages", renderLineages);
  safe("vesicles", renderVesicles); safe("levels", renderLevels); safe("trait", renderTrait); safe("machine", renderMachine);
  safe("darwin", renderDarwin); safe("check", renderCheck);
}
function renderHeader() {
  const st = S.status || {};
  const b = typeof st.banner === "string" && st.banner.trim() ? st.banner : (S.run ? "" : "No run yet. ") + FALLBACK_BANNER;
  const t = typeof st.time_note === "string" && st.time_note.trim() ? st.time_note : FALLBACK_TIME_NOTE;
  if ($("banner").textContent !== b) $("banner").textContent = b;
  if ($("time-note").textContent !== t) $("time-note").textContent = t;
  const state = st.state || (S.run ? "—" : "no run");
  const pill = $("run-state"); pill.textContent = state + (st.error ? ": " + st.error : ""); pill.className = "pill " + (st.state || "");
  const hw = S.hardware || {};
  const dev = get(st, "perf.device") || hw.device;
  const note = get(st, "perf.hw_note") || hw.hw_note;
  let line = [dev, note].filter((x) => x !== undefined && x !== null && x !== "").map(String).join(" · ");
  if (!line) {
    const extra = Object.entries(hw).filter(([, v]) => v !== null && typeof v !== "object").slice(0, 6).map(([k, v]) => `${k}: ${v}`);
    line = extra.join(" · ");
  }
  $("hw-line").textContent = "Hardware: " + (line || "—");
}
function renderTime() {
  const t = get(S.status, "time") || {}, v = get(S.status, "vesicles") || {};
  kv($("time-kv"), [
    ["State", S.status ? S.status.state : null],
    ["Sim hours", fmt(t.sim_h, 1)], ["Sim days", fmt(t.sim_days, 2)], ["Sim years", fmt(t.sim_years, 4)],
    ["Progress to 1 y", fmtPct(isNum(t.sim_years) ? Math.min(1, t.sim_years) : null)],
    ["Cycles", fmt(t.cycles, 2)], ["Cycle period", isNum(t.period_h) ? fmt(t.period_h, 1) + " h" : "—"],
    ["Vesicle gen. (max)", fmt(v.gen_max)], ["Molecular gen.", fmt(v.molecular_generations, 1)],
    ["Speed", isNum(t.speed_sim_h_per_wall_s) ? fmt(t.speed_sim_h_per_wall_s, 3) + " sim h / s" : "—"],
    ["Wall time", fmtDur(t.wall_s)], ["ETA to 1 sim year", isNum(t.eta_wall_h_to_1y) ? fmt(t.eta_wall_h_to_1y, 1) + " wall h" : "—"],
  ]);
}
function renderChains() {
  const st = S.status || {}, c = st.chains || {};
  kv($("chain-kv"), [
    ["Chains", fmt(c.n)], ["Mean length", fmt(c.mean_len, 1)], ["Max length", fmt(c.max_len)],
    ["L_MAX", fmt(st.lmax)], ["L_MAX auto max", fmt(st.lmax_auto_max)], ["W", fmt(st.W)],
    ["Nucleotide share", fmtPct(c.nuc_share)], ["Enclosed share", fmtPct(c.enclosed_share)],
  ]);
  lineChart($("c-len"), { logY: true, xlabel: xLabel(), series: [
    { pts: series("mean_len"), color: COL.c1 }, { pts: series("max_len"), color: COL.c2 }] });
  legend($("leg-len"), [{ color: COL.c1, label: "mean length" }, { color: COL.c2, label: "max length" }]);
  // histogram
  const { ctx, w, h } = prepCanvas($("c-hist"));
  const edges = get(c, "len_hist.edges"), counts = get(c, "len_hist.counts");
  if (!Array.isArray(counts) || !counts.length) { noData(ctx, w, h); return; }
  const L = 34, B = 30, T = 8, pw = w - L - 8, ph = h - T - B, n = counts.length;
  const mx = Math.max(1, ...counts.map((x) => (isNum(x) ? x : 0)));
  const sc = (x) => Math.log10(1 + (isNum(x) ? Math.max(0, x) : 0)) / Math.log10(1 + mx);
  const bw = pw / n;
  ctx.fillStyle = COL.axis; ctx.textAlign = "right"; ctx.textBaseline = "middle";
  ctx.fillText(tickLabel(mx), L - 4, T + 4); ctx.fillText("0", L - 4, T + ph);
  ctx.save(); ctx.translate(10, T + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textAlign = "center"; ctx.fillText("log count", 0, 0); ctx.restore();
  for (let i = 0; i < n; i++) {
    const bh = sc(counts[i]) * ph;
    ctx.fillStyle = COL.c1; ctx.fillRect(L + i * bw + 1, T + ph - bh, Math.max(1, bw - 2), bh);
    ctx.fillStyle = COL.axis; ctx.textAlign = "center"; ctx.textBaseline = "top";
    let lab = "";
    if (Array.isArray(edges)) lab = edges.length === n + 1 ? (isNum(edges[i]) ? tickLabel(edges[i]) : "") : (isNum(edges[i]) ? tickLabel(edges[i]) : "");
    if (bw > 18 || i % Math.ceil(18 / bw) === 0) ctx.fillText(lab, L + i * bw + bw / 2, T + ph + 4);
  }
  ctx.textAlign = "right"; ctx.fillText("length ≥ bin edge", L + pw, T + ph + 16);
}
function stripCanvas(seq) {
  const cv = document.createElement("canvas");
  const n = Math.min(128, seq.length), cell = 4;
  cv.width = 128 * cell; cv.height = 12; cv.style.height = "12px";
  const ctx = cv.getContext("2d");
  for (let i = 0; i < n; i++) {
    const ch = seq[i];
    ctx.fillStyle = ch === "0" ? "#2f6fd0" : ch === "1" ? "#f2b63c" : "#555c66";
    ctx.fillRect(i * cell, 0, cell - (cell > 3 ? 0.6 : 0), 12);
  }
  return cv;
}
function renderLineages() {
  const lin = get(S.status, "lineages") || {};
  const top = Array.isArray(lin.top) ? lin.top : [];
  kv($("lin-kv"), [["Distinct (≥ L_H)", fmt(lin.distinct_ge_LH)], ["Top shown", top.length ? String(Math.min(top.length, 12)) : "—"]]);
  lineChart($("c-distinct"), { xlabel: xLabel(), series: [{ pts: series("distinct"), color: COL.c3 }] });
  const box = $("top-seqs");
  const sig = JSON.stringify(top.slice(0, 12).map((t) => t && [t.key, t.count, t.score, t.len]));
  if (box.dataset.sig === sig) return;
  box.dataset.sig = sig; box.textContent = "";
  if (!top.length) { const d = document.createElement("div"); d.className = "note"; d.textContent = "—"; box.append(d); return; }
  top.slice(0, 12).forEach((t, i) => {
    if (!t || typeof t !== "object") return;
    const row = document.createElement("div"); row.className = "seq";
    const meta = document.createElement("div"); meta.className = "meta";
    const seq = typeof t.seq === "string" ? t.seq : "";
    meta.textContent = `#${i + 1}  len ${fmt(t.len)} · count ${fmt(t.count)} · score ${fmt(t.score, 3)} · key ${t.key || "—"}` + (seq.length > 128 ? `  (+${seq.length - 128} more letters)` : "");
    row.append(meta);
    if (seq) row.append(stripCanvas(seq));
    box.append(row);
  });
}
function renderVesicles() {
  const v = get(S.status, "vesicles") || {};
  kv($("ves-kv"), [
    ["Vesicles per site", fmt(v.per_site, 2)], ["Generation max", fmt(v.gen_max)], ["Generation mean", fmt(v.gen_mean, 1)],
    ["Divisions", fmt(v.divisions)], ["Molecular generations", fmt(v.molecular_generations, 1)],
  ]);
  lineChart($("c-ves"), { logY: true, xlabel: xLabel(), series: [
    { pts: series("vesicles"), color: COL.c4 }, { pts: series("gen_max"), color: COL.c2 }] });
  legend($("leg-ves"), [{ color: COL.c4, label: "vesicles" }, { color: COL.c2, label: "generation max" }]);
}
function renderLevels() {
  const lv = get(S.status, "levels") || {};
  const box = $("lamps"); box.textContent = "";
  for (const k of LEVELS) {
    const v = lv[k];
    const d = document.createElement("div"); d.className = "lamp " + (v === true ? "on" : v === false ? "off" : "unk");
    const dot = document.createElement("div"); dot.className = "dot";
    const l = document.createElement("div"); l.textContent = k.replace("_", " ") + (v === undefined || v === null ? " —" : "");
    d.append(dot, l); box.append(d);
  }
}
function renderTrait() {
  const st = S.status || {}, tr = st.trait || {}, ts = st.trait_stats || {};
  const beyond = isNum(ts.rho_wz) && isNum(ts.null99) ? (ts.rho_wz > ts.null99 ? "yes" : "no") : "—";
  kv($("trait-kv"), [
    ["Trait", tr.enabled === undefined ? "—" : (tr.enabled ? "on" : "off")], ["Mode", tr.mode], ["s", fmt(tr.s, 2)],
    ["L_rep", fmt(tr.L_rep)], ["k", fmt(tr.k)], ["Seed", fmt(tr.seed)],
    ["Fidelity link", tr.fidelity_link === undefined ? "—" : fmt(tr.fidelity_link)], ["s_f", fmt(tr.s_f, 2)],
    ["Score mean", fmt(ts.score_mean, 3)], ["Score max", fmt(ts.score_max, 3)], ["Score dominant", fmt(ts.score_dominant, 3)],
    ["ρ(w,z)", fmt(ts.rho_wz, 3)], ["null99", fmt(ts.null99, 3)], ["ρ beyond null99", beyond],
  ]);
  lineChart($("c-score"), { xlabel: xLabel(), series: [
    { pts: series("score_mean"), color: COL.c1 }, { pts: series("score_max"), color: COL.c2 }, { pts: series("score_dominant"), color: COL.c4 }] });
  legend($("leg-score"), [{ color: COL.c1, label: "mean" }, { color: COL.c2, label: "max" }, { color: COL.c4, label: "dominant" }]);
  lineChart($("c-rho"), { xlabel: xLabel(), series: [{ pts: series("rho_wz"), color: COL.c3 }],
    hlines: [{ y: ts.null99, color: COL.red, dash: true }, { y: 0, color: "#3a414b" }] });
  legend($("leg-rho"), [{ color: COL.c3, label: "ρ(w,z)" }, { color: COL.red, label: "null99 (current)", dash: true }]);
}
function renderMachine() {
  const p = get(S.status, "perf") || {};
  const temp = p.gpu_temp_c;
  kv($("mach-kv"), [
    ["Device", p.device], ["Steps", fmt(p.steps)], ["Step time", isNum(p.step_ms) ? fmt(p.step_ms, 1) + " ms" : "—"],
    ["GPU temp", isNum(temp) ? fmt(temp, 0) + " °C" + (temp >= 85 ? " (hot)" : "") : "—"],
    ["VRAM", isNum(p.vram_gib) ? fmt(p.vram_gib, 2) + " / " + fmt(p.vram_budget_gib, 1) + " GiB" : "—"],
    ["Hardware note", p.hw_note],
  ]);
  const frac = isNum(p.vram_gib) && isNum(p.vram_budget_gib) && p.vram_budget_gib > 0 ? Math.min(1, p.vram_gib / p.vram_budget_gib) : 0;
  const bar = $("vram-bar"); bar.style.width = (frac * 100).toFixed(1) + "%"; bar.style.background = frac > 0.9 ? COL.red : frac > 0.75 ? COL.c2 : COL.c1;
  const hs = S.hist.filter((r) => r && isNum(r.wall_s) && isNum(r.step_ms)).map((r) => [r.wall_s, r.step_ms]);
  lineChart($("c-step"), { xlabel: "wall s", series: [{ pts: hs.length ? hs : S.perf, color: COL.c2 }] });
}
function scoreColor(s, smax, smin = 0) {
  const t = isNum(s) ? Math.max(0, Math.min(1, smax > smin ? (s - smin) / (smax - smin) : 0)) : null;
  if (t === null) return "#666";
  // dark blue -> teal -> yellow
  const a = [40, 60, 140], b = [40, 170, 160], c = [250, 220, 70];
  const mix = (u, v, k) => u.map((x, i) => Math.round(x + (v[i] - x) * k));
  const rgb = t < 0.5 ? mix(a, b, t * 2) : mix(b, c, (t - 0.5) * 2);
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}
function renderDarwin() {
  const d = get(S.status, "darwin") || {};
  const status = typeof d.status === "string" ? d.status : null;
  const light = $("d-light");
  const cls = !status ? "off" : /sustained/i.test(status) ? "green" : /one sweep|selection on the trait/i.test(status) ? "amber" : /no selection/i.test(status) ? "red" : "off";
  light.className = "light " + cls;
  $("d-status").textContent = status || "—";
  const tr = d.trait_response || {};
  const trTxt = tr.diff === undefined ? "" : ` | trait response: real − control z = ${fmt(tr.diff, 3)}, control SD ${fmt(tr.sd_ctl, 3)}, ` +
    `consecutive passing reports ${fmt(tr.run)}, shifts ${fmt(tr.shifts)}${(tr.plateaus || []).length ? " (plateaus " + tr.plateaus.join(", ") + ")" : ""}`;
  $("d-note").textContent = (d.status_note || "sustained = ≥ 2 successive sweeps by a fitter lineage, each beyond the neutral null") + trTxt;
  // halves
  const box = $("d-halves"); box.textContent = "";
  for (const [key, title, color] of [["real", "Real (trait linked to sequence)", COL.c4], ["control", "Control (shuffled landscape)", COL.ctl]]) {
    const h = d[key] || {};
    const card = document.createElement("div"); card.className = "half";
    const h4 = document.createElement("h4"); h4.textContent = title; h4.style.color = color;
    const k = document.createElement("div"); k.className = "kv";
    kv(k, [["Score mean", fmt(h.score_mean, 3)], ["Score max", fmt(h.score_max, 3)], ["Score dominant", fmt(h.score_dominant, 3)],
      ["L4", h.L4 === undefined || h.L4 === null ? "—" : (h.L4 ? "reached" : "no")]]);
    const r = document.createElement("div");
    const known = isNum(h.rho_wz) && isNum(h.null99);
    r.className = "rho " + (known && h.rho_wz > h.null99 ? "beyond" : "within");
    r.textContent = `ρ(w,z) ${fmt(h.rho_wz, 3)} vs null99 ${fmt(h.null99, 3)} → ` + (known ? (h.rho_wz > h.null99 ? "beyond null" : "within null") : "—");
    card.append(h4, k, r); box.append(card);
  }
  // charts
  const sweeps = Array.isArray(d.sweeps) ? d.sweeps.filter((s) => s && typeof s === "object") : [];
  const markers = sweeps.map((s) => ({ x: s.sim_h, color: s.half === "control" ? COL.ctl : COL.c4, filled: !!s.beyond_null }));
  const maxR = get(d, "real.score_max"), maxC = get(d, "control.score_max");
  lineChart($("c-darwin-score"), { xlabel: xLabel(), markers, series: [
    { pts: series("score_mean_real"), color: COL.c1 }, { pts: series("score_mean_ctl"), color: COL.c1, dash: true },
    { pts: series("score_dom_real"), color: COL.c2 }, { pts: series("score_dom_ctl"), color: COL.c2, dash: true }],
    hlines: [{ y: maxR, color: COL.c4, dash: true }, { y: maxC, color: COL.ctl, dash: true }] });
  legend($("leg-darwin-score"), [
    { color: COL.c1, label: "mean (real)" }, { color: COL.c1, dash: true, label: "mean (control)" },
    { color: COL.c2, label: "dominant (real)" }, { color: COL.c2, dash: true, label: "dominant (control)" },
    { color: COL.c4, dash: true, label: "max now (real)" }, { color: COL.ctl, dash: true, label: "max now (control)" },
    { color: COL.c4, label: "▲ sweep real (filled = beyond null)" }, { color: COL.ctl, label: "▲ sweep control" }]);
  lineChart($("c-darwin-rho"), { xlabel: xLabel(), markers, series: [
    { pts: series("rho_wz_real"), color: COL.c4 }, { pts: series("null99_real"), color: COL.c4, dash: true },
    { pts: series("rho_wz_ctl"), color: COL.ctl }, { pts: series("null99_ctl"), color: COL.ctl, dash: true }],
    hlines: [{ y: 0, color: "#3a414b" }] });
  legend($("leg-darwin-rho"), [
    { color: COL.c4, label: "ρ real" }, { color: COL.c4, dash: true, label: "null99 real" },
    { color: COL.ctl, label: "ρ control" }, { color: COL.ctl, dash: true, label: "null99 control" }]);
  renderTree(get(d, "tree.nodes"));
  const last = sweeps.slice(-5).reverse();
  $("d-sweeps").textContent = sweeps.length
    ? `Sweeps: ${sweeps.length} total (${sweeps.filter((s) => s.half === "real").length} real, ${sweeps.filter((s) => s.half === "control").length} control). Latest: ` +
      last.map((s) => `${fmt(s.sim_h, 0)} h ${s.half || "?"} ${s.key || "?"} score ${fmt(s.score, 2)} share ${fmtPct(s.share)}${s.beyond_null ? " ✓null" : ""}`).join("; ")
    : "Sweeps: —";
  renderClades(d);
}
function renderClades(d) {
  const cl = d.clades || {};
  const ev = Array.isArray(d.clade_sweeps) ? d.clade_sweeps.filter((e) => e && typeof e === "object") : [];
  const markers = ev.map((e) => ({ x: e.sim_h, color: e.half === "control" ? COL.ctl : COL.c4, filled: e.verdict === "counts" }));
  lineChart($("c-clades"), { xlabel: xLabel(), markers, series: [
    { pts: series("clade_top1_share_real"), color: COL.c4 }, { pts: series("clade_top1_share_control"), color: COL.ctl, dash: true }],
    hlines: [{ y: 0.3, color: "#3a414b", dash: true }].concat(isNum(d.clade_ctl_p95) ? [{ y: d.clade_ctl_p95, color: COL.ctl }] : []) });
  legend($("leg-clades"), [
    { color: COL.c4, label: "top clade share (real)" }, { color: COL.ctl, dash: true, label: "top clade share (control)" },
    { color: "#3a414b", dash: true, label: "30 %" }, { color: COL.ctl, label: "control 95th percentile" },
    { color: COL.c4, label: "▲ clade event real (filled = counts)" }, { color: COL.ctl, label: "▲ clade event control" }]);
  const half = (h, name) => h ? `${name}: top clade ${fmtPct(h.top1_share)} over ${fmt(h.top1_sites)} sites, z ${fmt(h.z_mean, 2)}, ` +
    `${fmt(h.n_clades)} clades (census ${fmt(h.census)})` : `${name}: —`;
  const real = ev.filter((e) => e.half === "real");
  const last = real.slice(-5).reverse();
  $("d-clades").textContent = (cl.real || cl.control)
    ? `${half(cl.real, "Real")} | ${half(cl.control, "Control")} | control 95th pct ${isNum(d.clade_ctl_p95) ? fmtPct(d.clade_ctl_p95) : "pending"} ` +
      `(${fmt(d.clade_ctl_reports)} reports) | events: ${real.length} real (${fmt(d.clade_counted || 0)} counted), ` +
      `${ev.length - real.length} control` + (last.length ? ". Latest real: " + last.map((e) =>
        `${fmt(e.sim_h, 0)} h ${e.tag} ${fmtPct(e.share)} z ${fmt(e.z_mean, 2)} vs ${fmt(e.census_z_real, 2)} → ${e.verdict || "?"}`).join("; ") : "")
    : "Clades: — (no census yet)";
}
function renderTree(nodesIn) {
  const { ctx, w, h } = prepCanvas($("c-tree"));
  let nodes = Array.isArray(nodesIn) ? nodesIn.filter((n) => n && typeof n === "object" && n.id != null) : [];
  if (!nodes.length) { noData(ctx, w, h, "no lineage tree yet"); return; }
  if (nodes.length > 600) nodes = nodes.slice().sort((a, b) => (b.share || 0) - (a.share || 0)).slice(0, 600);
  const byId = new Map(nodes.map((n) => [String(n.id), n]));
  const halves = [["real", nodes.filter((n) => n.half !== "control")], ["control", nodes.filter((n) => n.half === "control")]].filter(([, a]) => a.length);
  let xmin = Infinity, xmax = -Infinity, smax = -Infinity, smin = Infinity;
  for (const n of nodes) { if (isNum(n.born_h)) { xmin = Math.min(xmin, n.born_h); xmax = Math.max(xmax, n.born_h); } if (isNum(n.score)) { smax = Math.max(smax, n.score); smin = Math.min(smin, n.score); } }
  if (!isFinite(smax)) { smin = 0; smax = 1; } if (smin > 0) smin = 0; if (smax <= smin) smax = smin + 1;
  const simh = get(S.status, "time.sim_h"); if (isNum(simh)) xmax = Math.max(xmax, simh);
  if (!isFinite(xmin)) { xmin = 0; xmax = 1; } if (xmax <= xmin) xmax = xmin + 1;
  const L = 70, R = 14, T = 8, B = 20, pw = w - L - R, ph = h - T - B;
  const X = (x) => L + ((isNum(x) ? x : xmin) - xmin) / (xmax - xmin) * pw;
  // order nodes per half by DFS so children sit next to their parent
  const pos = new Map(); let row = 0; const total = nodes.length + (halves.length - 1) * 2;
  const bands = [];
  for (const [name, arr] of halves) {
    const ids = new Set(arr.map((n) => String(n.id)));
    const kids = new Map();
    for (const n of arr) { const p = n.parent != null ? String(n.parent) : null; if (p && ids.has(p) && p !== String(n.id)) { if (!kids.has(p)) kids.set(p, []); kids.get(p).push(n); } }
    const roots = arr.filter((n) => !(n.parent != null && ids.has(String(n.parent)) && String(n.parent) !== String(n.id)));
    const seen = new Set(), start = row;
    const visit = (n, depth) => {
      const id = String(n.id); if (seen.has(id) || depth > 2000) return; seen.add(id);
      pos.set(id, row++);
      (kids.get(id) || []).slice().sort((a, b) => (a.born_h || 0) - (b.born_h || 0)).forEach((c) => visit(c, depth + 1));
    };
    roots.sort((a, b) => (a.born_h || 0) - (b.born_h || 0)).forEach((r) => visit(r, 0));
    for (const n of arr) if (!seen.has(String(n.id))) visit(n, 0);  // cycles
    bands.push([name, start, row - 1]); row += 2;
  }
  const Y = (r) => T + (total <= 1 ? ph / 2 : (r / (total - 1)) * ph);
  // band labels + axis
  ctx.textBaseline = "middle"; ctx.textAlign = "left";
  for (const [name, a, b] of bands) {
    ctx.fillStyle = name === "control" ? COL.ctl : COL.c4; ctx.fillText(name, 6, (Y(a) + Y(b)) / 2);
    ctx.strokeStyle = COL.grid; ctx.beginPath(); ctx.moveTo(L - 6, Y(a) - 3); ctx.lineTo(L - 6, Y(b) + 3); ctx.stroke();
  }
  ctx.fillStyle = COL.axis; ctx.textAlign = "center"; ctx.textBaseline = "top";
  for (const t of niceTicks(xmin, xmax, Math.max(2, Math.floor(pw / 90)))) ctx.fillText(tickLabel(t), X(t), T + ph + 4);
  ctx.textAlign = "left"; ctx.fillText("born (sim h)", L, T + 2);
  const rowH = total > 1 ? ph / (total - 1) : ph;
  // edges
  ctx.strokeStyle = "#4a535f"; ctx.lineWidth = 1;
  for (const n of nodes) {
    const p = n.parent != null ? byId.get(String(n.parent)) : null;
    if (!p || p === n || (p.half === "control") !== (n.half === "control")) continue;
    const x0 = X(p.born_h), y0 = Y(pos.get(String(p.id))), x1 = X(n.born_h), y1 = Y(pos.get(String(n.id)));
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0, y1); ctx.lineTo(x1, y1); ctx.stroke();
  }
  // sweep markers on the time axis
  const sw = get(S.status, "darwin.sweeps");
  if (Array.isArray(sw)) for (const s of sw) {
    if (!s || !isNum(s.sim_h)) continue;
    ctx.strokeStyle = s.half === "control" ? COL.ctl : COL.c4; ctx.globalAlpha = 0.35; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(X(s.sim_h), T); ctx.lineTo(X(s.sim_h), T + ph); ctx.stroke(); ctx.globalAlpha = 1; ctx.setLineDash([]);
  }
  // nodes
  const rmax = Math.max(3, Math.min(11, rowH * 0.9));
  for (const n of nodes) {
    const share = isNum(n.share) ? Math.max(0, Math.min(1, n.share)) : 0;
    const r = 2 + (rmax - 2) * Math.sqrt(share);
    ctx.fillStyle = scoreColor(n.score, smax, smin);
    ctx.beginPath(); ctx.arc(X(n.born_h), Y(pos.get(String(n.id))), r, 0, 7); ctx.fill();
    if (share > 0.3) { ctx.strokeStyle = "#fff"; ctx.lineWidth = 1; ctx.stroke(); }
  }
  // colour key
  const kx = L + pw - 130, ky = T + 2;
  for (let i = 0; i < 100; i++) { ctx.fillStyle = scoreColor(i / 99, 1); ctx.fillRect(kx + i, ky, 1, 7); }
  ctx.fillStyle = COL.axis; ctx.textAlign = "right"; ctx.textBaseline = "top"; ctx.fillText("score " + fmt(smin, 2), kx - 3, ky - 1);
  ctx.textAlign = "left"; ctx.fillText(fmt(smax, 2), kx + 103, ky - 1);
}

/* ---------------------------------------------------------------- init */
async function init() {
  try {
    const [params, presets] = await Promise.all([api("/api/params"), api("/api/presets")]);
    S.params = Array.isArray(params) ? params : [];
    S.presets = presets && typeof presets === "object" ? presets : {};
  } catch (e) { msg("Could not load params/presets: " + e.message, "err"); }
  const sel = $("preset"); sel.textContent = "";
  for (const [k, p] of Object.entries(S.presets)) {
    const o = document.createElement("option"); o.value = k; o.textContent = (p && p.label) || k; sel.append(o);
  }
  buildSliders(); safe("preset", applyPreset);
  sel.addEventListener("change", () => safe("preset", applyPreset));
  $("reset-env").addEventListener("click", () => safe("preset", applyPreset));
  $("t-s").addEventListener("input", () => { $("t-s-val").textContent = (+$("t-s").value).toFixed(2); });
  $("t-enabled").addEventListener("change", syncTraitEnabled);
  $("t-fid").addEventListener("change", syncTraitEnabled);
  $("btn-create").addEventListener("click", doCreate);
  $("btn-start").addEventListener("click", doStart);
  $("btn-pause").addEventListener("click", () => doControl("pause"));
  $("btn-resume").addEventListener("click", () => doControl("resume"));
  $("btn-stop").addEventListener("click", () => doControl("stop"));
  $("run-select").addEventListener("change", (e) => selectRun(e.target.value));
  let saved = ""; try { saved = localStorage.getItem("gk-run") || ""; } catch (e) { /* ignore */ }
  try { await refreshRuns(); } catch (e) { /* poll() reports it */ }
  if (saved && S.runs.some((r) => r.name === saved)) selectRun(saved);
  let rt = null;
  window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(renderAll, 150); });
  poll(); checkNow(); pollHardware();
  setInterval(poll, POLL_MS); setInterval(checkNow, CHECK_MS); setInterval(pollHardware, HW_MS);
}
init();
