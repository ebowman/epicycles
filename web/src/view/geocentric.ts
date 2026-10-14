import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { Elements, position } from "../math/kepler";
import type { Epicycle } from "../math/expansion";
import { chainArmsInto, chainTipInto, circleMatrix } from "../sim/geometry";
import { RingTrail } from "../sim/trail";
import { Label, makeGrid } from "./labels";

const MAX_CIRCLES = 242;      // 2 * (2*60 + 1)
const TRAIL_CAP = 400;
const TRAIL_STEP = 2;         // days between trail samples
const REF_HALF = 800;         // reference path covers t +- REF_HALF days
const REF_PTS = 401;

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
  private arms = new Float64Array(3 * (MAX_CIRCLES + 1));
  private armPos = new Float32Array(3 * (MAX_CIRCLES + 1));
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
    const sc = this.scene;
    sc.add(makeGrid(8, 16));
    const earth = new THREE.Mesh(new THREE.SphereGeometry(0.05, 16, 12), new THREE.MeshBasicMaterial({ color: 0x4aa3ff }));
    sc.add(earth, this.tip, this.truthMark);

    // Shared unit circle; per-circle matrices are written directly (matrixAutoUpdate = false).
    const circ = unitCircleGeometry();
    const matBody = new THREE.LineBasicMaterial({ color: 0xff9a5a, transparent: true, opacity: 0.28 });
    const matEarth = new THREE.LineBasicMaterial({ color: 0x5ab0ff, transparent: true, opacity: 0.28 });
    this.matBody = matBody; this.matEarth = matEarth;
    for (let i = 0; i < MAX_CIRCLES; i++) {
      const l = new THREE.LineLoop(circ, matBody);
      l.matrixAutoUpdate = false;
      l.frustumCulled = false;
      l.visible = false;
      this.pool.push(l);
      sc.add(l);
    }
    const ag = new THREE.BufferGeometry();
    ag.setAttribute("position", new THREE.BufferAttribute(this.armPos, 3));
    this.armLine = new THREE.Line(ag, new THREE.LineBasicMaterial({ color: 0xffe08a }));
    this.armLine.frustumCulled = false;
    sc.add(this.armLine);

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
    this.lEarth = new Label(host, "Earth (origin)");
    this.lTip = new Label(host, bodyEl.name ?? "Body");
    this.labels = [this.lEarth, this.lTip];
  }
  private matBody: THREE.Material; private matEarth: THREE.Material;

  get circleCount(): number { return this.chain.length; }

  setChain(chain: readonly Epicycle[]): void {
    this.chain = chain.length > MAX_CIRCLES ? chain.slice(0, MAX_CIRCLES) : chain;
    for (let i = 0; i < MAX_CIRCLES; i++) {
      const l = this.pool[i]!;
      l.visible = i < this.chain.length;
      if (l.visible) l.material = this.chain[i]!.source === "earth" ? this.matEarth : this.matBody;
    }
    this.armLine.geometry.setDrawRange(0, this.chain.length + 1);
    this.trailT = null; // trail is rebuilt (prefilled) on next update
  }

  private updateReference(t: number): void {
    if (!(Math.abs(t - this.refCenter) < REF_HALF / 4)) {
      this.refCenter = t;
      for (let i = 0; i < REF_PTS; i++) {
        const tau = t - REF_HALF + (2 * REF_HALF * i) / (REF_PTS - 1);
        const b = position(this.bodyEl, tau), e = position(this.earthEl, tau);
        this.refPos[3 * i] = b[0] - e[0]; this.refPos[3 * i + 1] = b[1] - e[1]; this.refPos[3 * i + 2] = b[2] - e[2];
      }
      (this.refGeo.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    }
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
    this.updateReference(t);
    this.updateTrail(t);
    this.lTip.pos.copy(this.tip.position);
    this.controls.update();
    for (const l of this.labels) l.update(this.camera, w, h);
  }
}
