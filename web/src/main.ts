import * as THREE from "three";

// Placeholder: proves three.js imports and a canvas renders. Real views come later.
const app = document.getElementById("app")!;
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
app.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.1, 100);
camera.position.z = 3;
const ring = new THREE.Mesh(
  new THREE.TorusGeometry(1, 0.01, 8, 128),
  new THREE.MeshBasicMaterial({ color: 0x66aaff }),
);
scene.add(ring);
function frame(): void {
  ring.rotation.x += 0.01;
  renderer.render(scene, camera);
  requestAnimationFrame(frame);
}
frame();
