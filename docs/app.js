"use strict";

const C = 299792458.0;
const OMEGA_EARTH = 7.2921150e-5;   // rad/s
const GPS_UNIX_EPOCH = 315964800;   // unix seconds at GPS epoch 1980-01-06
const BG = "#d6dae0", SCENE_BG = "#dde1e7", FG = "#1a1a1a";
const N_T_MAX = 500, N_F_MAX = 400;

let EPH = null, DETS = null, COAST = null;

// ---------------------------------------------------------------------------
// Time and geometry
// ---------------------------------------------------------------------------

function leapOffset(t) {            // GPS - UTC, valid from 2009 onward
  if (t >= 1167264018) return 18;   // 2017-01-01
  if (t >= 1119744017) return 17;   // 2015-07-01
  if (t >= 1025136016) return 16;   // 2012-07-01
  return 15;
}

function gmst(t) {                  // Greenwich mean sidereal time (rad), IAU 1982
  const jd = 2444244.5 + (t - leapOffset(t)) / 86400.0;
  const T = (jd - 2451545.0) / 36525.0;
  const s = 67310.54841 + (876600.0 * 3600 + 8640184.812866) * T + 0.093104 * T * T - 6.2e-6 * T * T * T;
  const r = (s * 2 * Math.PI / 86400.0) % (2 * Math.PI);
  return r < 0 ? r + 2 * Math.PI : r;
}

function gpsToUtcString(t) {
  const d = new Date((GPS_UNIX_EPOCH + t - leapOffset(t)) * 1000);
  return d.toISOString().slice(0, 19).replace("T", " ") + " UTC";
}

function catmullRom(arr, idx, comp, u) {
  const n = EPH.n;
  const i1 = Math.min(Math.max(idx, 0), n - 1);
  const i0 = Math.max(i1 - 1, 0), i2 = Math.min(i1 + 1, n - 1), i3 = Math.min(i1 + 2, n - 1);
  const p0 = arr[3 * i0 + comp], p1 = arr[3 * i1 + comp], p2 = arr[3 * i2 + comp], p3 = arr[3 * i3 + comp];
  return 0.5 * ((2 * p1) + (-p0 + p2) * u + (2 * p0 - 5 * p1 + 4 * p2 - p3) * u * u + (-p0 + 3 * p1 - 3 * p2 + p3) * u * u * u);
}

function earthPosVel(t) {           // [x,y,z] AU and [vx,vy,vz] m/s, barycentric, ICRS axes
  let x = (t - EPH.t0) / EPH.dt;
  x = Math.min(Math.max(x, 0), EPH.n - 1.000001);
  const i = Math.floor(x), u = x - i;
  const pos = [0, 1, 2].map(c => catmullRom(EPH.pos_au, i, c, u));
  const vel = [0, 1, 2].map(c => catmullRom(EPH.vel_ms, i, c, u));
  return [pos, vel];
}

function sitePosVel(t, det) {       // detector position (m) and velocity (m/s) about Earth's centre
  const g = gmst(t), cg = Math.cos(g), sg = Math.sin(g);
  const [x, y, z] = det.location_m;
  const r = [x * cg - y * sg, x * sg + y * cg, z];
  const v = [-OMEGA_EARTH * r[1], OMEGA_EARTH * r[0], 0];
  return [r, v];
}

function unitVec(alpha, delta) {
  return [Math.cos(delta) * Math.cos(alpha), Math.cos(delta) * Math.sin(alpha), Math.sin(delta)];
}

function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function norm(a) { return Math.sqrt(dot(a, a)); }
function scale(a, s) { return [a[0] * s, a[1] * s, a[2] * s]; }

function antenna(t, det, ra, dec, psi) {   // F+, Fx  (same algorithm as LAL ComputeDetAMResponse)
  const gha = gmst(t) - ra;
  const cg = Math.cos(gha), sg = Math.sin(gha), cd = Math.cos(dec), sd = Math.sin(dec);
  const cp = Math.cos(psi), sp = Math.sin(psi);
  const X = [-cp * sg - sp * cg * sd, -cp * cg + sp * sg * sd, sp * cd];
  const Y = [sp * sg - cp * cg * sd, sp * cg + cp * sg * sd, cp * cd];
  const D = det.response;
  let fp = 0, fc = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    fp += (X[i] * X[j] - Y[i] * Y[j]) * D[i][j];
    fc += (X[i] * Y[j] + Y[i] * X[j]) * D[i][j];
  }
  return [fp, fc];
}

function dopplerTrack(ts, p) {
  const det = DETS[p.detector], n = unitVec(p.Alpha, p.Delta);
  const out = new Float64Array(ts.length);
  for (let i = 0; i < ts.length; i++) {
    const t = ts[i];
    const [, vo] = earthPosVel(t);
    let v = vo;
    if (p.spin) { const [, vs] = sitePosVel(t, det); v = [vo[0] + vs[0], vo[1] + vs[1], vo[2] + vs[2]]; }
    const beta = dot(v, n) / C;
    const dt = t - p.tref;
    out[i] = (p.F0 + p.F1 * dt + 0.5 * p.F2 * dt * dt) * (1 + beta);
  }
  return out;
}

function mulberry32(seed) {
  return function () {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function emulate(p, tRange, fRange) {
  const Tsft = p.Tsft, det = DETS[p.detector];
  const k0 = Math.max(Math.floor((tRange[0] - p.tstart) / Tsft), 0);
  const k1 = Math.max(Math.ceil((tRange[1] - p.tstart) / Tsft), k0 + 1);
  let ks = [];
  const nk = k1 - k0 + 1;
  if (nk > N_T_MAX) {
    const seen = new Set();
    for (let i = 0; i < N_T_MAX; i++) { const k = Math.round(k0 + (k1 - k0) * i / (N_T_MAX - 1)); if (!seen.has(k)) { seen.add(k); ks.push(k); } }
  } else for (let k = k0; k <= k1; k++) ks.push(k);
  const t = ks.map(k => p.tstart + (k + 0.5) * Tsft);

  const dF = 1.0 / Tsft;
  const j0 = Math.floor(fRange[0] / dF), j1 = Math.max(Math.ceil(fRange[1] / dF), j0 + 1);
  let f = [];
  const nj = j1 - j0 + 1;
  if (nj > N_F_MAX) for (let i = 0; i < N_F_MAX; i++) f.push((j0 + (j1 - j0) * i / (N_F_MAX - 1)) * dF);
  else for (let j = j0; j <= j1; j++) f.push(j * dF);
  const dfPix = f.length > 1 ? f[1] - f[0] : dF;
  const width = Math.max(dF, dfPix);

  const track = dopplerTrack(t, p);
  const aPlus = p.h0 * (1 + p.cosi * p.cosi) / 2, aCross = p.h0 * p.cosi;
  const rng = mulberry32(12345);
  const power = new Array(t.length);
  for (let i = 0; i < t.length; i++) {
    const [fp, fc] = antenna(t[i], det, p.Alpha, p.Delta, p.psi);
    const a2 = (aPlus * fp) ** 2 + (aCross * fc) ** 2;
    const amp2 = a2 * Tsft / (2 * p.sqrtSX * p.sqrtSX);
    const row = new Float32Array(f.length);
    for (let j = 0; j < f.length; j++) {
      const x = Math.PI * (f[j] - track[i]) / width;
      const s = Math.abs(x) < 1e-9 ? 1 : Math.sin(x) / x;
      row[j] = amp2 * s * s + (p.noise ? -Math.log(1 - rng()) : 0);
    }
    power[i] = row;
  }
  return { t, f, power, track };
}

// ---------------------------------------------------------------------------
// State and UI
// ---------------------------------------------------------------------------

const p = {
  F0: 100, F1: -1e-9, F2: 0, Alpha: 1.5, Delta: 0.5, h0: 1e-23, sqrtSX: 1e-23, cosi: 0, psi: 0,
  Tsft: 1800, duration: 500 * 86400, tstart: 1400000000, tref: 1400000000,
  detector: "H1", spin: true, track: true, noise: true,
};
let view = null;     // {t: [d0, d1] days, f: [f0, f1] Hz}
let tmark = 0;       // days

function fullView() {
  const n = 400, ts = new Float64Array(n);
  for (let i = 0; i < n; i++) ts[i] = p.tstart + p.duration * i / (n - 1);
  const tr = dopplerTrack(ts, p);
  let lo = Infinity, hi = -Infinity;
  for (const v of tr) { if (v < lo) lo = v; if (v > hi) hi = v; }
  const pad = Math.max(0.15 * (hi - lo), 20 / p.Tsft);
  return { t: [0, p.duration / 86400], f: [lo - pad, hi + pad] };
}

function percentile(rows, q) {
  const sample = [];
  const step = Math.max(1, Math.floor(rows.length * rows[0].length / 20000));
  let c = 0;
  for (const r of rows) for (let j = 0; j < r.length; j++) { if (c++ % step === 0) sample.push(r[j]); }
  sample.sort((a, b) => a - b);
  return sample[Math.min(sample.length - 1, Math.floor(q * sample.length))];
}

function renderSpec() {
  const t0 = p.tstart + view.t[0] * 86400, t1 = p.tstart + view.t[1] * 86400;
  const r = emulate(p, [t0, t1], view.f);
  const tDays = r.t.map(v => (v - p.tstart) / 86400);
  const zmax = Math.max(percentile(r.power, 0.995), 1);
  const data = [{
    type: "heatmap", x: tDays, y: r.f, z: r.power, transpose: true,
    colorscale: "Viridis", zmin: 0, zmax, colorbar: { title: { text: "power / noise" } },
    hovertemplate: "t = %{x:.3f} d<br>f = %{y:.6f} Hz<br>P = %{z:.1f}<extra></extra>",
  }];
  if (p.track) data.push({ type: "scatter", mode: "lines", x: tDays, y: Array.from(r.track),
    line: { color: "#ff4444", width: 1 }, hoverinfo: "skip" });
  const layout = {
    shapes: [{ type: "line", x0: tmark, x1: tmark, y0: 0, y1: 1, yref: "paper",
               line: { color: "white", dash: "dash", width: 1.5 } }],
    xaxis: { title: { text: "Time from start (days)" }, range: view.t.slice() },
    yaxis: { title: { text: "Frequency (Hz)" }, range: view.f.slice(), tickformat: ".6f", exponentformat: "none" },
    margin: { l: 90, r: 10, t: 30, b: 40 }, paper_bgcolor: BG, plot_bgcolor: "#440154",
    font: { color: FG }, showlegend: false, dragmode: "zoom",
    title: { text: `${p.detector}   Tsft = ${p.Tsft} s   bin = ${(1e3 / p.Tsft).toFixed(3)} mHz   pixels ${r.t.length} x ${r.f.length}`,
             font: { size: 13 } },
  };
  Plotly.react("spec", data, layout, { scrollZoom: true, displaylogo: false, responsive: true });
}

// ---- 3D views -------------------------------------------------------------

function sphere(radius, n = 40) {
  const x = [], y = [], z = [];
  for (let i = 0; i < n; i++) {
    const u = 2 * Math.PI * i / (n - 1); const rx = [], ry = [], rz = [];
    for (let j = 0; j < n / 2; j++) {
      const v = Math.PI * j / (n / 2 - 1);
      rx.push(radius * Math.cos(u) * Math.sin(v)); ry.push(radius * Math.sin(u) * Math.sin(v)); rz.push(radius * Math.cos(v));
    }
    x.push(rx); y.push(ry); z.push(rz);
  }
  return { x, y, z };
}
const SPHERE = sphere(1.0);

function arrow(tail, head, color, name, width = 5) {
  const d = [head[0] - tail[0], head[1] - tail[1], head[2] - tail[2]];
  return [
    { type: "scatter3d", mode: "lines", x: [tail[0], head[0]], y: [tail[1], head[1]], z: [tail[2], head[2]],
      line: { color, width }, name, hoverinfo: "name" },
    { type: "cone", x: [head[0]], y: [head[1]], z: [head[2]], u: [d[0]], v: [d[1]], w: [d[2]],
      sizemode: "absolute", sizeref: 0.25 * norm(d), anchor: "tip",
      colorscale: [[0, color], [1, color]], showscale: false, hoverinfo: "skip" },
  ];
}

function sceneLayout(uirev, rng, title) {
  const ax = { visible: false, range: [-rng, rng] };
  return {
    uirevision: uirev, title: { text: title, font: { size: 13 } },
    scene: { xaxis: ax, yaxis: ax, zaxis: ax, aspectmode: "cube", bgcolor: SCENE_BG, uirevision: uirev },
    paper_bgcolor: BG, font: { color: FG }, margin: { l: 0, r: 0, t: 30, b: 0 },
    legend: { x: 0, y: 1, font: { size: 10 }, bgcolor: "rgba(255,255,255,0.6)" },
  };
}

function lonlatToXyz(lonDeg, latDeg, g, r = 1.003) {
  const lam = lonDeg * Math.PI / 180 + g, phi = latDeg * Math.PI / 180;
  return [r * Math.cos(phi) * Math.cos(lam), r * Math.cos(phi) * Math.sin(lam), r * Math.sin(phi)];
}

function coastTrace(g) {
  const x = [], y = [], z = [];
  for (const seg of COAST) {
    for (const [lon, lat] of seg) { const q = lonlatToXyz(lon, lat, g); x.push(q[0]); y.push(q[1]); z.push(q[2]); }
    x.push(null); y.push(null); z.push(null);
  }
  return { type: "scatter3d", mode: "lines", x, y, z, line: { color: "#e8f0d8", width: 2 }, hoverinfo: "skip", showlegend: false };
}

function renderEarth(tgps) {
  const n = unitVec(p.Alpha, p.Delta);
  const [posE] = earthPosVel(tgps);
  const g = gmst(tgps);
  const sunDir = scale(posE, -1 / norm(posE));
  const sites = {};
  for (const name in DETS) { const [r] = sitePosVel(tgps, DETS[name]); sites[name] = scale(r, 1.02 / norm(r)); }

  const traces = [{
    type: "surface", x: SPHERE.x, y: SPHERE.y, z: SPHERE.z, surfacecolor: SPHERE.z,
    colorscale: [[0, "#0b3d91"], [1, "#2d7dd2"]], showscale: false, hoverinfo: "skip", lighting: { ambient: 0.7 },
  }, coastTrace(g)];
  const gx = [], gy = [], gz = [];
  for (let i = 0; i < 60; i++) { const q = lonlatToXyz(0, -90 + 180 * i / 59, g); gx.push(q[0]); gy.push(q[1]); gz.push(q[2]); }
  traces.push({ type: "scatter3d", mode: "lines", x: gx, y: gy, z: gz, line: { color: "#ffcc00", width: 3 }, hoverinfo: "skip", name: "Greenwich meridian" });
  const ex = [], ey = [], ez = [];
  for (let i = 0; i < 90; i++) { const q = lonlatToXyz(360 * i / 89, 0, g); ex.push(q[0]); ey.push(q[1]); ez.push(q[2]); }
  traces.push({ type: "scatter3d", mode: "lines", x: ex, y: ey, z: ez, line: { color: "#ffffff", width: 2 }, hoverinfo: "skip", name: "Equator" });

  const others = Object.keys(sites).filter(k => k !== p.detector);
  traces.push({ type: "scatter3d", mode: "markers+text", x: others.map(k => sites[k][0]), y: others.map(k => sites[k][1]), z: others.map(k => sites[k][2]),
    marker: { size: 5, color: "#000000", line: { color: "#fff", width: 1 } }, text: others, textposition: "top center",
    textfont: { color: "#000000", size: 11 }, name: "Other detectors", hoverinfo: "text" });
  const s = sites[p.detector];
  traces.push({ type: "scatter3d", mode: "markers+text", x: [s[0]], y: [s[1]], z: [s[2]],
    marker: { size: 8, color: "#ff2d2d", line: { color: "#fff", width: 1 } }, text: [p.detector], textposition: "top center",
    textfont: { color: "#ff2d2d", size: 13 }, name: `${p.detector} (selected)` });
  traces.push(...arrow(scale(n, 1.9), scale(n, 1.15), "#ff4444", "CW propagation (from source)"));
  traces.push(...arrow([0, 0, 0], scale(sunDir, 1.5), "#ff9900", "To Sun", 4));
  traces.push({ type: "scatter3d", mode: "lines", x: [0, 1.4 * n[0]], y: [0, 1.4 * n[1]], z: [0, 1.4 * n[2]],
    line: { color: "#ff4444", width: 2, dash: "dot" }, showlegend: false, hoverinfo: "skip" });

  Plotly.react("earth", traces, sceneLayout("earth", 1.5, `Earth view (GCRS axes)   GPS ${tgps.toFixed(0)}   ${gpsToUtcString(tgps)}`),
    { displaylogo: false, responsive: true });
}

let ringCache = { key: null, ring: null };
function orbitRing() {
  if (ringCache.key !== p.tstart) {
    const ring = [];
    for (let i = 0; i < 200; i++) ring.push(earthPosVel(p.tstart + 365.25 * 86400 * i / 199)[0]);
    ringCache = { key: p.tstart, ring };
  }
  return ringCache.ring;
}

function renderOrbit(tgps) {
  const n = unitVec(p.Alpha, p.Delta);
  const [posE] = earthPosVel(tgps);
  const ring = orbitRing();
  const path = [];
  for (let i = 0; i < 120; i++) path.push(earthPosVel(p.tstart + Math.max(tgps - p.tstart, 1) * i / 119)[0]);
  const traces = [
    { type: "scatter3d", mode: "lines", x: ring.map(q => q[0]), y: ring.map(q => q[1]), z: ring.map(q => q[2]),
      line: { color: "#99a", width: 2 }, name: "Orbit (1 yr)", hoverinfo: "skip" },
    { type: "scatter3d", mode: "lines", x: path.map(q => q[0]), y: path.map(q => q[1]), z: path.map(q => q[2]),
      line: { color: "#1aa34a", width: 4 }, name: "Observed so far", hoverinfo: "skip" },
    { type: "scatter3d", mode: "markers", x: [0], y: [0], z: [0], marker: { size: 12, color: "#ff9900" }, name: "Sun" },
    { type: "scatter3d", mode: "markers", x: [posE[0]], y: [posE[1]], z: [posE[2]], marker: { size: 6, color: "#2d7dd2" }, name: "Earth" },
    ...arrow(scale(n, 2.6), scale(n, 1.6), "#ff4444", "CW propagation (from source)"),
  ];
  Plotly.react("orbit", traces, sceneLayout("orbit", 2.6, "Orbit view (SSB, AU)"), { displaylogo: false, responsive: true });
}

// ---- wiring ---------------------------------------------------------------

let specPending = false, viewsPending = false;
function scheduleSpec(reset) {
  if (reset) view = fullView();
  if (specPending) return;
  specPending = true;
  requestAnimationFrame(() => { specPending = false; renderSpec(); });
}
function scheduleViews() {
  if (viewsPending) return;
  viewsPending = true;
  requestAnimationFrame(() => {
    viewsPending = false;
    const tgps = p.tstart + tmark * 86400;
    renderEarth(tgps); renderOrbit(tgps);
  });
}

function setParam(id, value) {
  if (id === "duration") p.duration = value * 86400; else p[id] = value;
}

function wirePanel() {
  for (const ctl of document.querySelectorAll(".ctl")) {
    const id = ctl.dataset.id, log = ctl.dataset.log === "1";
    const slider = ctl.querySelector("input[type=range]"), box = ctl.querySelector("input[type=number]");
    slider.addEventListener("input", () => {
      const v = log ? Math.pow(10, +slider.value) : +slider.value;
      box.value = log ? v.toExponential(2) : v;
      setParam(id, v); afterParamChange(id);
    });
    box.addEventListener("change", () => {
      const v = +box.value;
      if (!isFinite(v)) return;
      if (log) { if (v <= 0) return; slider.value = Math.log10(v); } else slider.value = v;
      setParam(id, v); afterParamChange(id);
    });
  }
  const det = document.getElementById("detector");
  det.addEventListener("change", () => { p.detector = det.value; afterParamChange("detector"); });
  const ts = document.getElementById("tstart");
  const showDate = () => { document.getElementById("tstart-date").textContent = gpsToUtcString(p.tstart); };
  ts.addEventListener("change", () => { const v = +ts.value; if (!isFinite(v)) return; p.tstart = v; p.tref = v; showDate(); afterParamChange("tstart"); });
  showDate();
  for (const [id, key] of [["opt-spin", "spin"], ["opt-track", "track"], ["opt-noise", "noise"]]) {
    const el = document.getElementById(id);
    el.addEventListener("change", () => { p[key] = el.checked; afterParamChange(key); });
  }
  const tsl = document.getElementById("tslider"), tbox = document.getElementById("tbox");
  const setMark = v => {
    tmark = Math.min(Math.max(v, 0), p.duration / 86400);
    tsl.value = tmark; tbox.value = tmark.toFixed(2);
    Plotly.relayout("spec", { "shapes[0].x0": tmark, "shapes[0].x1": tmark });
    scheduleViews();
  };
  tsl.addEventListener("input", () => setMark(+tsl.value));
  tbox.addEventListener("change", () => setMark(+tbox.value));
}

function afterParamChange(id) {
  if (p.Tsft <= 0 || p.duration <= 0 || p.sqrtSX <= 0 || p.h0 < 0 || p.F0 <= 0) return;
  const dur = p.duration / 86400;
  const tsl = document.getElementById("tslider");
  tsl.max = dur;
  if (tmark > dur) { tmark = dur; tsl.value = dur; document.getElementById("tbox").value = dur.toFixed(2); }
  scheduleSpec(true);
  scheduleViews();
}

function wireSpec() {
  const el = document.getElementById("spec");
  el.on("plotly_relayout", ev => {
    if (!ev || view === null) return;
    let changed = false;
    const nv = { t: view.t.slice(), f: view.f.slice() };
    if ("xaxis.range[0]" in ev) { nv.t = [ev["xaxis.range[0]"], ev["xaxis.range[1]"]]; changed = true; }
    else if ("xaxis.range" in ev) { nv.t = ev["xaxis.range"].slice(); changed = true; }
    if ("yaxis.range[0]" in ev) { nv.f = [ev["yaxis.range[0]"], ev["yaxis.range[1]"]]; changed = true; }
    else if ("yaxis.range" in ev) { nv.f = ev["yaxis.range"].slice(); changed = true; }
    if (ev["xaxis.autorange"] || ev["yaxis.autorange"]) {
      const full = fullView();
      if (ev["xaxis.autorange"]) nv.t = full.t;
      if (ev["yaxis.autorange"]) nv.f = full.f;
      changed = true;
    }
    if (!changed) return;
    nv.t[0] = Math.max(nv.t[0], 0); nv.t[1] = Math.min(nv.t[1], p.duration / 86400);
    if (nv.t[1] <= nv.t[0] || nv.f[1] <= nv.f[0]) return;
    view = nv;
    scheduleSpec(false);
  });
}

async function main() {
  const [eph, dets, coast] = await Promise.all([
    fetch("data/ephemeris.json").then(r => r.json()),
    fetch("data/detectors.json").then(r => r.json()),
    fetch("data/coastlines_110m.json").then(r => r.json()),
  ]);
  EPH = eph; DETS = dets; COAST = coast;
  wirePanel();
  view = fullView();
  renderSpec();
  wireSpec();
  scheduleViews();
}

window.cwsky = { dopplerTrack, antenna, earthPosVel, sitePosVel, gmst, emulate, p: () => p };
main();
