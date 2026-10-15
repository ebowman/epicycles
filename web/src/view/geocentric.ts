import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { Elements, position } from "../math/kepler";
import type { Epicycle } from "../math/expansion";
import { chainArmsInto, chainTipInto, circleMatrix, fitDistance } from "../sim/geometry";
import { RingTrail } from "../sim/trail";
import { Label, makeGrid } from "./labels";

const INITIAL_POOL = 242;     // 2 * (2*60 + 1); grows on demand (Phaethon N=300 -> 1202)
const TRAIL_CAP = 400;
const TRAIL_STEP = 2;         // days between trail samples
const REF_HALF = 500;         // reference path covers t +- REF_HALF days
const REF_PTS = 401;
/** Circle opacity tiers by radius ratio r/rmax: >0.1, >0.01, >0.001, smaller. */
const TIER_OPACITY = [0.5, 0.26, 0.12, 0.05];
function tierOf(ratio: number): number { return ratio > 0.1 ? 0 : ratio > 0.01 ? 1 : ratio > 0.001 ? 2 : 3; }

function unitCircleGeometry(): THREE.BufferGeometry {
  const n = 64, a = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) { a[3 * i] = Math.cos((2 * Math.PI * i) / n); a[3 * i + 1] = Math.sin((2 * Math.PI * i) / n); }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(a, 3));
  return g;
}

export class GeoView {
  readonly scene = new THREE.Scene();
  readonly camera = new THREE.PerspectiveCamera(40, 1, 0.01, 200);
  readonly controls: OrbitControls;
  private chain: readonly Epicycle[] = [];
  private pool: THREE.LineLoop[] = [];
  private circ: THREE.BufferGeometry;
  private mats: Record<string, THREE.LineBasicMaterial[]> = {};
  private arms = new Float64Array(3 * (INITIAL_POOL + 1));
  private armPos = new Float32Array(3 * (INITIAL_POOL + 1));
  private showCircles = true;
  private errLine: THREE.Line;
  private errPos = new Float32Array(6);
  private aspect = 1;
  private earthMark: THREE.Mesh;
  private refFn: (tau: number, out: Float64Array) => void;
  private scan = new Float64Array(3);
  private armLine: THREE.Line;
  private tip = new THREE.Mesh(new THREE.SphereGeometry(0.03, 14, 10), new THREE.MeshBasicMaterial({ color: 0xff5ad0 }));
  private truthMark = new THREE.Mesh(new THREE.SphereGeometry(0.018, 10, 8), new THREE.MeshBasicMaterial({ color: 0xffffff }));
  private trail = new RingTrail(TRAIL_CAP);
  private trailBuf = new Float32Array(3 * (TRAIL_CAP + 1));
  private trailGeo = new THREE.BufferGeometry();
  private trailT: number | null = null;
  private refGeo = new THREE.BufferGeometry();
  private refPos = new Float32Array(3 * REF_PTS);
  private refCenter = Number.NaN;
  private tmp = new Float64Array(3);
  private lTip: Label; private lEarth: Label;
  private labels: Label[];

  constructor(host: HTMLElement, dom: HTMLElement, private earthEl: Elements, private bodyEl: Elements) {
    this.camera.far = 400; this.camera.updateProjectionMatrix();
    const sc = this.scene;
    sc.add(makeGrid(8, 16));
    const earth = this.earthMark = new THREE.Mesh(new THREE.SphereGeometry(0.05, 16, 12), new THREE.MeshBasicMaterial({ color: 0x4aa3ff }));
    sc.add(earth, this.tip, this.truthMark);

    // Shared unit circle; per-circle matrices are written directly (matrixAutoUpdate = false).
    // Circles are bucketed into opacity tiers by radius relative to the largest, so the many tiny
    // circles at high N fade out instead of cluttering the retrograde loops.
    this.circ = unitCircleGeometry();
    const colours: Record<string, number> = { body: 0xff9a5a, naff: 0xff9a5a, earth: 0x5ab0ff };
    for (const src of Object.keys(colours)) {
      this.mats[src] = TIER_OPACITY.map((o) => new THREE.LineBasicMaterial({ color: colours[src]!, transparent: true, opacity: o, depthWrite: false }));
    }
    this.growPool(INITIAL_POOL);
    const ag = new THREE.BufferGeometry();
    ag.setAttribute("position", new THREE.BufferAttribute(this.armPos, 3));
    this.armLine = new THREE.Line(ag, new THREE.LineBasicMaterial({ color: 0xffe08a }));
    this.armLine.frustumCulled = false;
    sc.add(this.armLine);

    this.errLine = new THREE.Line(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: 0xff3b3b }));
    this.errLine.geometry.setAttribute("position", new THREE.BufferAttribute(this.errPos, 3));
    this.errLine.frustumCulled = false;
    this.errLine.visible = false;
    sc.add(this.errLine);

    // Tip trail with per-vertex fade (oldest = background colour).
    this.trailGeo.setAttribute("position", new THREE.BufferAttribute(this.trailBuf, 3));
    const col = new Float32Array(3 * (TRAIL_CAP + 1));
    for (let i = 0; i <= TRAIL_CAP; i++) {
      const f = 0.25 + 0.75 * Math.pow(i / TRAIL_CAP, 1.5);
      col[3 * i] = 1.0 * f; col[3 * i + 1] = 0.35 * f; col[3 * i + 2] = 0.82 * f;
    }
    this.trailGeo.setAttribute("color", new THREE.BufferAttribute(col, 3));
    const trailLine = new THREE.Line(this.trailGeo, new THREE.LineBasicMaterial({ vertexColors: true }));
    trailLine.frustumCulled = false;
    sc.add(trailLine);

    // True geocentric path (reference).
    this.refGeo.setAttribute("position", new THREE.BufferAttribute(this.refPos, 3));
    const refLine = new THREE.Line(this.refGeo, new THREE.LineBasicMaterial({ color: 0x9fffd0, transparent: true, opacity: 0.7 }));
    refLine.frustumCulled = false;
    sc.add(refLine);

    this.camera.up.set(0, 0, 1);
    this.camera.position.set(0, -3.4, 6.8);
    this.controls = new OrbitControls(this.camera, dom);
    this.controls.enableDamping = true;
    this.refFn = (tau, out) => {
      const b = position(this.bodyEl, tau), e = position(this.earthEl, tau);
      out[0] = b[0] - e[0]; out[1] = b[1] - e[1]; out[2] = b[2] - e[2];
    };
    this.lEarth = new Label(host, "Earth (origin)");
    this.lTip = new Label(host, bodyEl.name ?? "Body");
    this.frameReference(2452880);
    this.labels = [this.lEarth, this.lTip];
  }
  get circleCount(): number { return this.chain.length; }

  private growPool(n: number): void {
    while (this.pool.length < n) {
      const l = new THREE.LineLoop(this.circ, this.mats.body![0]!);
      l.matrixAutoUpdate = false;
      l.frustumCulled = false;
      l.visible = false;
      this.pool.push(l);
      this.scene.add(l);
    }
    if (this.arms.length < 3 * (n + 1)) {
      this.arms = new Float64Array(3 * (n + 1));
      this.armPos = new Float32Array(3 * (n + 1));
      const old = this.armLine.geometry;
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(this.armPos, 3));
      this.armLine.geometry = g;
      old.dispose();
    }
  }

  /** Set the epicycle chain (any length; the circle pool grows as needed). */
  setChain(chain: readonly Epicycle[]): void {
    this.chain = chain;
    this.growPool(chain.length);
    let rmax = 0;
    for (const c of chain) rmax = Math.max(rmax, c.radius);
    for (let i = 0; i < this.pool.length; i++) {
      const l = this.pool[i]!;
      l.visible = this.showCircles && i < chain.length;
      if (i < chain.length) {
        const c = chain[i]!;
        l.material = this.mats[c.source]![tierOf(rmax > 0 ? c.radius / rmax : 0)]!;
      }
    }
    this.armLine.geometry.setDrawRange(0, chain.length + 1);
    this.trailT = null; // trail is rebuilt (prefilled) on next update
  }

  setShowCircles(on: boolean): void {
    this.showCircles = on;
    for (let i = 0; i < this.pool.length; i++) this.pool[i]!.visible = on && i < this.chain.length;
  }
  setShowError(on: boolean): void { this.errLine.visible = on; }

  /** Switch the reference (true) body; resets the reference path. refFn overrides it (NAFF mode). */
  setBody(bodyEl: Elements, label: string, refFn?: (tau: number, out: Float64Array) => void): void {
    this.bodyEl = bodyEl;
    this.lTip.el.textContent = label;
    this.refFn = refFn ?? ((tau, out) => {
      const b = position(this.bodyEl, tau), e = position(this.earthEl, tau);
      out[0] = b[0] - e[0]; out[1] = b[1] - e[1]; out[2] = b[2] - e[2];
    });
    this.refCenter = Number.NaN;
    this.trailT = null;
  }

  private updateReference(t: number): void {
    if (!(Math.abs(t - this.refCenter) < REF_HALF / 4)) {
      this.refCenter = t;
      const o = this.scan;
      for (let i = 0; i < REF_PTS; i++) {
        this.refFn(t - REF_HALF + (2 * REF_HALF * i) / (REF_PTS - 1), o);
        this.refPos[3 * i] = o[0]!; this.refPos[3 * i + 1] = o[1]!; this.refPos[3 * i + 2] = o[2]!;
      }
      (this.refGeo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    }
  }

  /** Place the camera to view an axis-aligned box (AU) top-down with a slight tilt. */
  private frameBox(min: ArrayLike<number>, max: ArrayLike<number>): void {
    const cx = (min[0]! + max[0]!) / 2, cy = (min[1]! + max[1]!) / 2, cz = (min[2]! + max[2]!) / 2;
    const hw = Math.max(0.02, (max[0]! - min[0]!) / 2), hh = Math.max(0.02, (max[1]! - min[1]!) / 2);
    const dist = fitDistance(hw, hh, this.camera.fov, this.aspect);
    const d = new THREE.Vector3(0, -0.3, 0.954);
    this.controls.target.set(cx, cy, cz);
    this.camera.position.set(cx + d.x * dist, cy + d.y * dist, cz + d.z * dist);
    this.camera.lookAt(this.controls.target);
    this.controls.update();
  }

  /** Frame the camera to the bounding box of the true geocentric path around time t. */
  frameReference(t: number): void {
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    const o = this.scan;
    for (let i = 0; i < REF_PTS; i++) {
      this.refFn(t - REF_HALF + (2 * REF_HALF * i) / (REF_PTS - 1), o);
      for (let k = 0; k < 3; k++) { min[k] = Math.min(min[k]!, o[k]!); max[k] = Math.max(max[k]!, o[k]!); }
    }
    this.includeCircles(t, min, max);
    this.frameBox(min, max);
  }

  /** Grow a framing box so the epicycle circles drawn at time t (centre +- radius, in-plane) stay in view. */
  private includeCircles(t: number, min: number[], max: number[]): void {
    const arms = new Float64Array(3 * (this.chain.length + 1));
    chainArmsInto(this.chain, t, arms);
    for (let j = 0; j < this.chain.length; j++) {
      const r = this.chain[j]!.radius;
      for (let k = 0; k < 2; k++) {
        min[k] = Math.min(min[k]!, arms[3 * j + k]! - r); max[k] = Math.max(max[k]!, arms[3 * j + k]! + r);
      }
    }
  }

  /**
   * Centre/zoom on the retrograde loop (decreasing geocentric longitude) containing t, or the nearest one
   * within +-REF_HALF days. Returns false (and frames the whole reference path) if there is none.
   */
  zoomToOpposition(t: number): boolean {
    const o = this.scan, step = 2;
    const n = Math.floor((2 * REF_HALF) / step);
    const lon = new Float64Array(n + 1), pts = new Float64Array(3 * (n + 1));
    let prev = 0;
    for (let i = 0; i <= n; i++) {
      this.refFn(t - REF_HALF + i * step, o);
      pts[3 * i] = o[0]!; pts[3 * i + 1] = o[1]!; pts[3 * i + 2] = o[2]!;
      let a = Math.atan2(o[1]!, o[0]!);
      if (i > 0) { while (a - prev > Math.PI) a -= 2 * Math.PI; while (a - prev < -Math.PI) a += 2 * Math.PI; }
      lon[i] = a; prev = a;
    }
    let best = -1, bestDist = Infinity, start = -1;
    for (let i = 1; i <= n + 1; i++) {
      const retro = i <= n && lon[i]! < lon[i - 1]!;
      if (retro && start < 0) start = i - 1;
      if (!retro && start >= 0) {
        const end = i - 1, mid = (start + end) / 2, d = Math.abs(t - (t - REF_HALF + mid * step));
        if (end > start && d < bestDist) { bestDist = d; best = start * 100000 + end; }
        start = -1;
      }
    }
    if (best < 0) { this.frameReference(t); return false; }
    const pad = Math.round(40 / step);
    const s0 = Math.max(0, Math.floor(best / 100000) - pad), e0 = Math.min(n, (best % 100000) + pad);
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    for (let i = s0; i <= e0; i++) for (let k = 0; k < 3; k++) {
      min[k] = Math.min(min[k]!, pts[3 * i + k]!); max[k] = Math.max(max[k]!, pts[3 * i + k]!);
    }
    // keep the box at least half as tall as wide (and vice versa) so a thin arc is not over-zoomed
    const dx = max[0]! - min[0]!, dy = max[1]! - min[1]!, big = Math.max(dx, dy, 0.3);
    for (const k of [0, 1] as const) {
      const ext = max[k]! - min[k]!, want = Math.max(1.5 * ext, 0.8 * big), c = (max[k]! + min[k]!) / 2;
      min[k] = c - want / 2; max[k] = c + want / 2;
    }
    this.frameBox(min, max);
    return true;
  }

  private updateTrail(t: number): void {
    const tmp = this.tmp;
    if (this.trailT === null || t < this.trailT || t - this.trailT > 40 * TRAIL_STEP) {
      this.trail.clear();
      for (let i = TRAIL_CAP - 1; i >= 0; i--) {
        chainTipInto(this.chain, t - i * TRAIL_STEP, tmp);
        this.trail.push(tmp[0]!, tmp[1]!, tmp[2]!);
      }
      this.trailT = t;
    } else {
      while (this.trailT + TRAIL_STEP <= t) {
        this.trailT += TRAIL_STEP;
        chainTipInto(this.chain, this.trailT, tmp);
        this.trail.push(tmp[0]!, tmp[1]!, tmp[2]!);
      }
    }
    const n = this.trail.copyOrdered(this.trailBuf);
    const a = this.arms, m = 3 * this.chain.length;
    this.trailBuf[3 * n] = a[m]!; this.trailBuf[3 * n + 1] = a[m + 1]!; this.trailBuf[3 * n + 2] = a[m + 2]!;
    this.trailGeo.setDrawRange(0, n + 1);
    (this.trailGeo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
  }

  /** Updates everything for time t; returns the chain tip via `tipOut` (x,y,z). */
  update(t: number, truth: ArrayLike<number>, tipOut: Float64Array, w: number, h: number): void {
    const chain = this.chain, arms = this.arms, n = chain.length;
    this.aspect = w / Math.max(1, h);
    const ms = Math.min(1.5, Math.max(0.02, this.camera.position.distanceTo(this.controls.target) / 7));
    this.earthMark.scale.setScalar(ms); this.tip.scale.setScalar(ms); this.truthMark.scale.setScalar(ms);
    chainArmsInto(chain, t, arms);
    for (let j = 0; j < n; j++) {
      const l = this.pool[j]!;
      circleMatrix(chain[j]!, arms[3 * j]!, arms[3 * j + 1]!, arms[3 * j + 2]!, l.matrix.elements);
      l.matrixWorldNeedsUpdate = true;
    }
    const ap = this.armPos;
    for (let i = 0; i < 3 * (n + 1); i++) ap[i] = arms[i]!;
    (this.armLine.geometry.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    tipOut[0] = arms[3 * n]!; tipOut[1] = arms[3 * n + 1]!; tipOut[2] = arms[3 * n + 2]!;
    this.tip.position.set(tipOut[0], tipOut[1], tipOut[2]);
    this.truthMark.position.set(truth[0]!, truth[1]!, truth[2]!);
    if (this.errLine.visible) {
      const e = this.errPos;
      e[0] = tipOut[0]!; e[1] = tipOut[1]!; e[2] = tipOut[2]!; e[3] = truth[0]!; e[4] = truth[1]!; e[5] = truth[2]!;
      (this.errLine.geometry.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    }
    this.updateReference(t);
    this.updateTrail(t);
    this.lTip.pos.copy(this.tip.position);
    this.controls.update();
    for (const l of this.labels) l.update(this.camera, w, h);
  }
}
