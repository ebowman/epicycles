import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { Elements, position, meanMotion } from "../math/kepler";
import { Label, makeGrid } from "./labels";

function orbitLine(el: Elements, color: number): THREE.LineLoop {
  const n = 256, P = (2 * Math.PI) / meanMotion(el);
  const pts = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) {
    const p = position(el, el.epoch + (P * i) / n);
    pts[3 * i] = p[0]; pts[3 * i + 1] = p[1]; pts[3 * i + 2] = p[2];
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(pts, 3));
  return new THREE.LineLoop(g, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.6 }));
}

export class HelioView {
  readonly scene = new THREE.Scene();
  readonly camera = new THREE.PerspectiveCamera(40, 1, 0.01, 200);
  readonly controls: OrbitControls;
  private earth = new THREE.Mesh(new THREE.SphereGeometry(0.035, 16, 12), new THREE.MeshBasicMaterial({ color: 0x4aa3ff }));
  private body = new THREE.Mesh(new THREE.SphereGeometry(0.03, 16, 12), new THREE.MeshBasicMaterial({ color: 0xff8a4a }));
  private link: THREE.Line;
  private bodyOrbit: THREE.LineLoop;
  private linkPos = new Float32Array(6);
  private labels: Label[];
  private lEarth: Label; private lBody: Label;

  constructor(host: HTMLElement, dom: HTMLElement, earthEl: Elements, bodyEl: Elements) {
    const sc = this.scene;
    sc.add(makeGrid(6, 12));
    const sun = new THREE.Mesh(new THREE.SphereGeometry(0.07, 20, 14), new THREE.MeshBasicMaterial({ color: 0xffd84a }));
    this.bodyOrbit = orbitLine(bodyEl, 0xff8a4a);
    sc.add(sun, this.earth, this.body, orbitLine(earthEl, 0x4aa3ff), this.bodyOrbit);
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(this.linkPos, 3));
    this.link = new THREE.Line(g, new THREE.LineBasicMaterial({ color: 0xffffff }));
    this.link.frustumCulled = false;
    sc.add(this.link);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(0, -3.6, 4.2);
    this.controls = new OrbitControls(this.camera, dom);
    this.controls.enableDamping = true;
    const lSun = new Label(host, "Sun");
    this.lEarth = new Label(host, "Earth"); this.lBody = new Label(host, bodyEl.name ?? "Body");
    this.labels = [lSun, this.lEarth, this.lBody];
    this.earthEl = earthEl; this.bodyEl = bodyEl;
  }
  private earthEl: Elements; private bodyEl: Elements;

  /** Switch the body whose orbit is drawn; fits the camera distance to the aphelion. */
  setBody(bodyEl: Elements, label: string): void {
    this.bodyEl = bodyEl;
    this.scene.remove(this.bodyOrbit);
    this.bodyOrbit.geometry.dispose();
    this.bodyOrbit = orbitLine(bodyEl, 0xff8a4a);
    this.scene.add(this.bodyOrbit);
    this.lBody.el.textContent = label;
    const aph = bodyEl.a * (1 + bodyEl.e);
    const dist = Math.max(3, 4.0 * aph);
    this.controls.target.set(0, 0, 0);
    this.camera.position.set(0, -0.65 * dist, 0.76 * dist);
    this.controls.update();
  }

  update(t: number, w: number, h: number): void {
    const pe = position(this.earthEl, t), pb = position(this.bodyEl, t);
    this.earth.position.set(pe[0], pe[1], pe[2]);
    this.body.position.set(pb[0], pb[1], pb[2]);
    const a = this.linkPos;
    a[0] = pe[0]; a[1] = pe[1]; a[2] = pe[2]; a[3] = pb[0]; a[4] = pb[1]; a[5] = pb[2];
    (this.link.geometry.getAttribute("position") as THREE.BufferAttribute).needsUpdate = true;
    this.lEarth.pos.copy(this.earth.position); this.lBody.pos.copy(this.body.position);
    this.controls.update();
    for (const l of this.labels) l.update(this.camera, w, h);
  }
}
