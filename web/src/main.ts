/**
 * Design: ONE WebGLRenderer on a full-window fixed canvas; the two views are transparent DOM
 * <div>s laid out by CSS grid (side by side, stacked below 800px). Each frame the renderer's
 * viewport+scissor are set to each div's bounding rect and the view's scene is rendered there.
 * OrbitControls attach to the divs (which sit above the canvas), so each view has its own camera.
 */
import * as THREE from "three";
import { EARTH, MARS, J2000 } from "./math/bodies";
import { geocentricChain } from "./math/expansion";
import { position } from "./math/kepler";
import { Clock, MAX_SPEED, MIN_SPEED } from "./sim/clock";
import { AppState } from "./sim/state";
import { computeError, ErrorReadout, formatAU, formatLonError, jdToDate } from "./sim/geometry";
import { HelioView } from "./view/heliocentric";
import { GeoView } from "./view/geocentric";

const $ = <T extends HTMLElement>(id: string): T => document.getElementById(id) as T;

function webglAvailable(): boolean {
  try {
    const c = document.createElement("canvas");
    return !!(c.getContext("webgl2") || c.getContext("webgl"));
  } catch { return false; }
}

function main(): void {
  if (!webglAvailable()) { $("nogl").style.display = "flex"; return; }
  const q = new URLSearchParams(location.search);
  const num = (k: string, d: number): number => { const v = Number(q.get(k)); return q.has(k) && Number.isFinite(v) ? v : d; };
  const state = new AppState(new Clock(num("t", 2452880), num("speed", 30), q.get("paused") !== "1"), num("N", 20));
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

  const vL = $("vL"), vR = $("vR");
  const helio = new HelioView(vL, vL, EARTH, MARS);
  const geo = new GeoView(vR, vR, EARTH, MARS);
  const rebuild = (): void => geo.setChain(geocentricChain(MARS, EARTH, state.N, state.N));
  rebuild();

  // controls
  const play = $<HTMLButtonElement>("play"), speed = $<HTMLInputElement>("speed"), scrub = $<HTMLInputElement>("scrub");
  const nsl = $<HTMLInputElement>("nsl");
  const lg = Math.log10(MAX_SPEED / MIN_SPEED);
  speed.value = String(Math.log10(clock.speed / MIN_SPEED) / lg);
  nsl.value = String(state.N);
  play.onclick = (): void => { clock.toggle(); };
  speed.oninput = (): void => clock.setSpeed(MIN_SPEED * Math.pow(10, Number(speed.value) * lg));
  scrub.oninput = (): void => clock.scrub(J2000 + Number(scrub.value));
  nsl.oninput = (): void => { if (state.setN(Number(nsl.value))) rebuild(); };

  const err: ErrorReadout = { au: 0, lonArcsec: 0 };
  const tip = new Float64Array(3);
  const truth = [0, 0, 0];
  const vp = new THREE.Vector4();
  let last = performance.now(), frames = 0, fpsT = last, fps = 0;
  const readout = $("readout"), date = $("date"), speedv = $("speedv"), nv = $("nv");

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
    const t = clock.t;

    const pm = position(MARS, t), pe = position(EARTH, t);
    truth[0] = pm[0] - pe[0]; truth[1] = pm[1] - pe[1]; truth[2] = pm[2] - pe[2];

    const rl = vL.getBoundingClientRect(), rr = vR.getBoundingClientRect();
    helio.update(t, rl.width, rl.height);
    geo.update(t, truth, tip, rr.width, rr.height);
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
      `circles ${geo.circleCount} | err ${formatAU(err.au)} | lon err ${formatLonError(err.lonArcsec)} | ${fps.toFixed(0)} fps`;
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}

main();
