import * as THREE from "three";

/** DOM text label pinned to a 3D point; positioned via transform each frame (no allocation). */
export class Label {
  readonly el: HTMLElement;
  readonly pos = new THREE.Vector3();
  private v = new THREE.Vector3();
  constructor(host: HTMLElement, text: string, cls = "lbl") {
    this.el = document.createElement("span");
    this.el.className = cls;
    this.el.textContent = text;
    host.appendChild(this.el);
  }
  update(camera: THREE.Camera, w: number, h: number): void {
    this.v.copy(this.pos).project(camera);
    if (this.v.z > 1 || this.v.z < -1) { this.el.style.display = "none"; return; }
    this.el.style.display = "";
    this.el.style.transform =
      `translate(${((this.v.x * 0.5 + 0.5) * w + 6).toFixed(1)}px,${((-this.v.y * 0.5 + 0.5) * h - 6).toFixed(1)}px)`;
  }
}

export function makeGrid(size: number, divisions: number): THREE.GridHelper {
  const g = new THREE.GridHelper(size, divisions, 0x3a4560, 0x232a3d);
  g.rotation.x = Math.PI / 2; // GridHelper lies in XZ; ecliptic is XY (z = ecliptic north)
  return g;
}
