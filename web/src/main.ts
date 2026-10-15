/**
 * Design: ONE WebGLRenderer on a full-window fixed canvas; the two views are transparent DOM
 * <div>s laid out by CSS grid (side by side, stacked below 800px). Each frame the renderer's
 * viewport+scissor are set to each div's bounding rect and the view's scene is rendered there.
 * OrbitControls attach to the divs (which sit above the canvas), so each view has its own camera.
 */
import * as THREE from "three";
import { EARTH, J2000 } from "./math/bodies";
import { Epicycle, geocentricChain } from "./math/expansion";
import { position } from "./math/kepler";
import { Clock, MAX_SPEED, MIN_SPEED } from "./sim/clock";
import { AppState } from "./sim/state";
import { chainTipInto, computeError, ErrorReadout, formatAU, formatLonError, jdToDate } from "./sim/geometry";
import {
  BIG_N_WARN, BODY_ELEMENTS, BODY_LABEL, BodyId, NAFF_HELIO_ELEMENTS, NAFF_INTERIOR, NAFF_WINDOW,
  PresetId, circleCountFor, ptolemaicChain, clampToWindow, fourierN, nMaxFor, parseParams, presetDefaults, selectTopK, studyRows,
} from "./sim/presets";
import { LoadedChain, loadBundledNaffChain } from "./data/loader";
import { buildInfo } from "./ui/info";
import { HelioView } from "./view/heliocentric";
import { GeoView } from "./view/geocentric";

const $ = <T extends HTMLElement>(id: string): T => document.getElementById(id) as T;

function webglAvailable(): boolean {
  try {
    const c = document.createElement("canvas");
    return !!(c.getContext("webgl2") || c.getContext("webgl"));
  } catch { return false; }
}

async function main(): Promise<void> {
  if (!webglAvailable()) { $("nogl").style.display = "flex"; return; }
  const init = parseParams(location.search);
  const state = new AppState(new Clock(init.t, init.speed, !init.paused), init.N);
  const clock = state.clock;

  const canvas = $<HTMLCanvasElement>("gl");
  let renderer: THREE.WebGLRenderer;
  try { renderer = new THREE.WebGLRenderer({ canvas, antialias: true }); }
  catch { $("nogl").style.display = "flex"; return; }
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0x0b0d14, 1);
  renderer.setScissorTest(true);
  const resize = (): void => renderer.setSize(window.innerWidth, window.innerHeight, false);
  window.addEventListener("resize", resize);
  resize();

  const vL = $("vL"), vR = $("vR"), infoEl = $("info");
  const info = buildInfo(infoEl);
  const helio = new HelioView(vL, vL, EARTH, BODY_ELEMENTS.mars);
  const geo = new GeoView(vR, vR, EARTH, BODY_ELEMENTS.mars);

  // controls
  const play = $<HTMLButtonElement>("play"), speed = $<HTMLInputElement>("speed"), scrub = $<HTMLInputElement>("scrub");
  const nsl = $<HTMLInputElement>("nsl"), bodySel = $<HTMLSelectElement>("body"), warn = $("warn");
  const chkCircles = $<HTMLInputElement>("chkCircles"), chkError = $<HTMLInputElement>("chkError");
  const lg = Math.log10(MAX_SPEED / MIN_SPEED);
  speed.value = String(Math.log10(clock.speed / MIN_SPEED) / lg);

  let body: BodyId = "mars";
  let naff: LoadedChain | null = null;
  let naffFull: Epicycle[] = [];
  let ptolemaic = false;
  let pendingFrame: "ref" | "opp" | null = "ref";
  // Re-run the last automatic framing when the geocentric viewport is resized (e.g. the info panel
  // finishes laying out) unless the user has since moved the camera.
  let lastFrame: "ref" | "opp" = "ref", userMoved = false;
  geo.controls.addEventListener("start", () => { userMoved = true; });
  new ResizeObserver(() => { if (!userMoved && pendingFrame === null) pendingFrame = lastFrame; }).observe(vR);
  let fourierStart: number | null = null;
  let warnTimer = 0;
  const flash = (msg: string): void => { warn.textContent = msg; clearTimeout(warnTimer); warnTimer = window.setTimeout(updateWarn, 4000); };

  const rebuild = (): void => {
    if (body === "naff") geo.setChain(selectTopK(naffFull, state.N));
    else {
      const full = geocentricChain(BODY_ELEMENTS[body], EARTH, state.N, state.N);
      geo.setChain(ptolemaic && state.N === 1 ? ptolemaicChain(full) : full);
    }
    updateWarn();
  };
  function updateWarn(): void {
    const c = circleCountFor(body, state.N);
    warn.textContent = body !== "naff" && state.N > BIG_N_WARN
      ? `Warning: N=${state.N} gives ${c} circle terms (4N rotating + 2 static, one combined offset); rendering may be slow.` : "";
  }

  async function setBody(id: BodyId): Promise<void> {
    if (id === "naff" && !naff) {
      try { naff = loadBundledNaffChain(); naffFull = naff.epicycles; }
      catch (e) { flash(`Could not load bundled perturbed_ceres chain: ${(e as Error).message}`); id = "mars"; }
    }
    body = id;
    bodySel.value = id;
    state.setNMax(nMaxFor(id));
    nsl.max = String(state.nMax);
    $("nlabel").textContent = id === "naff" ? "K" : "N";
    nsl.setAttribute("aria-label", id === "naff"
      ? "Keep the top K harmonics: 2K+1 largest NAFF circles (arrow keys)" : "Number of epicycle harmonics N (arrow keys)");
    if (id === "naff") {
      geo.setBody(NAFF_HELIO_ELEMENTS, "Ceres (NAFF chain)", (tau, out) => chainTipInto(naffFull, tau, out));
      helio.setBody(NAFF_HELIO_ELEMENTS, "Ceres (osculating Kepler)");
      scrub.min = String(NAFF_WINDOW[0] - J2000); scrub.max = String(NAFF_WINDOW[1] - J2000); scrub.step = "5";
      clock.scrub(clampToWindow(clock.t, id));
    } else {
      const el = BODY_ELEMENTS[id];
      geo.setBody(el, BODY_LABEL[id]);
      helio.setBody(el, BODY_LABEL[id]);
      scrub.min = "-7305"; scrub.max = "7305"; scrub.step = "0.5";
    }
    nsl.value = String(state.N);
    rebuild();
    pendingFrame = "ref";
    const rows = naff && id === "naff" ? studyRows(naff.study) : [];
    const note = id === "naff"
      ? "Perturbed Ceres: N-body integration (Sun, Earth, Jupiter, Saturn) fitted by NAFF into 401 circles; the slider keeps the top K by radius (2K+1 circle terms including the static offset). Reference path = full NAFF reconstruction (all 401 circles), not the N-body truth; see paper for N-body truth error and the table below. Left view shows the unperturbed osculating Kepler ellipse for context."
      : `${BODY_LABEL[id]}: analytic Bessel/Fourier expansion of the Kepler orbit, computed live. Reference path = exact geocentric Kepler path. ${circleCountFor(id, state.N)} circle terms.`;
    const win = id === "naff"
      ? `Study window ${jdToDate(NAFF_WINDOW[0])} to ${jdToDate(NAFF_WINDOW[1])} (clock limited to it); interior (middle 80%) ${jdToDate(NAFF_INTERIOR[0])} to ${jdToDate(NAFF_INTERIOR[1])}.` : "";
    info.setMode(note, rows, win);
  }

  function setErrorShown(on: boolean): void { chkError.checked = on; geo.setShowError(on); }
  const chkPtol = $<HTMLInputElement>("chkPtol");
  function setPtolemaic(on: boolean): void { ptolemaic = on; chkPtol.checked = on; rebuild(); }

  async function applyPreset(p: PresetId): Promise<void> {
    const d = presetDefaults(p);
    fourierStart = null;
    state.setNMax(nMaxFor(d.body ?? "mars"));
    state.setN(d.N ?? 20);
    await setBody(d.body ?? "mars");
    state.setN(d.N ?? 20); nsl.value = String(state.N);
    setPtolemaic(d.ptolemaic ?? false);
    clock.scrub(d.t ?? clock.t); clock.play();
    setErrorShown(d.showError ?? false);
    if (d.fourier) fourierStart = performance.now();
  }

  await setBody(init.body);
  state.setN(init.N); nsl.value = String(state.N); rebuild();
  geo.setShowCircles(init.showCircles); chkCircles.checked = init.showCircles;
  setErrorShown(init.showError);
  setPtolemaic(init.ptolemaic);
  if (init.fourier) { state.setN(1); nsl.value = "1"; rebuild(); fourierStart = performance.now(); }
  if (init.zoom === "opposition") pendingFrame = "opp";
  const showInfo = init.info && window.innerWidth > 800;
  infoEl.classList.toggle("hidden", !showInfo);
  $("infoBtn").setAttribute("aria-expanded", String(showInfo));

  const stepN = (d: number): void => { fourierStart = null; if (state.setN(state.N + d)) { nsl.value = String(state.N); rebuild(); } };
  play.onclick = (): void => { clock.toggle(); };
  speed.oninput = (): void => clock.setSpeed(MIN_SPEED * Math.pow(10, Number(speed.value) * lg));
  scrub.oninput = (): void => clock.scrub(clampToWindow(J2000 + Number(scrub.value), body));
  nsl.oninput = (): void => { fourierStart = null; if (state.setN(Number(nsl.value))) rebuild(); };
  bodySel.onchange = (): void => { fourierStart = null; setPtolemaic(false); void setBody(bodySel.value as BodyId); };
  chkCircles.onchange = (): void => geo.setShowCircles(chkCircles.checked);
  chkPtol.onchange = (): void => setPtolemaic(chkPtol.checked);
  chkError.onchange = (): void => geo.setShowError(chkError.checked);
  $("zoomOpp").onclick = (): void => { pendingFrame = "opp"; };
  $("frameAll").onclick = (): void => { pendingFrame = "ref"; };
  $("pPtolemy").onclick = (): void => { void applyPreset("ptolemy"); };
  $("pFourier").onclick = (): void => { void applyPreset("fourier"); };
  $("pNaff").onclick = (): void => { void applyPreset("ceres-naff"); };
  $("infoBtn").onclick = (): void => {
    const hide = !infoEl.classList.contains("hidden");
    infoEl.classList.toggle("hidden", hide);
    $("infoBtn").setAttribute("aria-expanded", String(!hide));
  };
  document.addEventListener("keydown", (ev) => {
    const el = ev.target as HTMLElement, tag = el.tagName;
    if (ev.ctrlKey || ev.metaKey || ev.altKey) return;
    if (ev.key === " ") {
      if (tag === "BUTTON" || tag === "SELECT" || (tag === "INPUT" && (el as HTMLInputElement).type === "checkbox")) return;
      ev.preventDefault(); clock.toggle();
    } else if (tag !== "INPUT" && tag !== "SELECT" && tag !== "TEXTAREA") {
      if (ev.key === "ArrowRight" || ev.key === "ArrowUp") { ev.preventDefault(); stepN(1); }
      else if (ev.key === "ArrowLeft" || ev.key === "ArrowDown") { ev.preventDefault(); stepN(-1); }
    }
  });

  const err: ErrorReadout = { au: 0, lonArcsec: 0 };
  const tip = new Float64Array(3);
  const truth = new Float64Array(3);
  const tmpB = new Float64Array(3);
  const vp = new THREE.Vector4();
  let last = performance.now(), frames = 0, fpsT = last, fps = 0;
  const readout = $("readout"), date = $("date"), speedv = $("speedv"), nv = $("nv");

  function truthAt(t: number): void {
    if (body === "naff") { chainTipInto(naffFull, t, truth); return; }
    const pm = position(BODY_ELEMENTS[body], t), pe = position(EARTH, t);
    tmpB[0] = pm[0] - pe[0]; tmpB[1] = pm[1] - pe[1]; tmpB[2] = pm[2] - pe[2];
    truth.set(tmpB);
  }

  function drawView(el: HTMLElement, scene: THREE.Scene, cam: THREE.PerspectiveCamera): { w: number; h: number } {
    const r = el.getBoundingClientRect();
    const w = Math.max(1, Math.floor(r.width)), h = Math.max(1, Math.floor(r.height));
    vp.set(Math.floor(r.left), Math.floor(window.innerHeight - r.bottom), w, h);
    renderer.setViewport(vp.x, vp.y, vp.z, vp.w);
    renderer.setScissor(vp.x, vp.y, vp.z, vp.w);
    cam.aspect = w / h;
    cam.updateProjectionMatrix();
    renderer.render(scene, cam);
    return { w, h };
  }

  function frame(now: number): void {
    const dt = Math.min((now - last) / 1000, 0.25);
    last = now;
    clock.step(dt);
    if (body === "naff" && (clock.t > NAFF_WINDOW[1] || clock.t < NAFF_WINDOW[0])) clock.scrub(NAFF_WINDOW[0]);
    if (fourierStart !== null) {
      const n = fourierN((now - fourierStart) / 1000);
      if (state.setN(n)) { nsl.value = String(state.N); rebuild(); }
      if (n >= 20) fourierStart = null;
    }
    const t = clock.t;
    truthAt(t);

    const rl = vL.getBoundingClientRect(), rr = vR.getBoundingClientRect();
    helio.update(t, rl.width, rl.height);
    geo.update(t, truth, tip, rr.width, rr.height);
    if (pendingFrame !== null) {
      let framed: "ref" | "opp" = pendingFrame;
      if (pendingFrame === "opp") {
        if (!geo.zoomToOpposition(t)) { flash("No retrograde loop within +-800 days; framing the whole path."); framed = "ref"; }
      } else geo.frameReference(t);
      lastFrame = framed; userMoved = false;
      pendingFrame = null;
    }
    computeError(tip, truth, err);

    drawView(vL, helio.scene, helio.camera);
    drawView(vR, geo.scene, geo.camera);

    frames++;
    if (now - fpsT > 1000) { fps = (frames * 1000) / (now - fpsT); frames = 0; fpsT = now; }
    play.textContent = clock.playing ? "Pause" : "Play";
    scrub.value = String(t - J2000);
    date.textContent = jdToDate(t);
    speedv.textContent = `${clock.speed.toFixed(0)} d/s`;
    nv.textContent = `${state.N}`;
    readout.textContent =
      `circle terms ${geo.circleCount} | err ${formatAU(err.au)} | lon err ${formatLonError(err.lonArcsec)} | ${fps.toFixed(0)} fps` +
      (body === "naff" ? " | ref = full NAFF" : "");
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}

main().catch((e) => { console.error(e); });
