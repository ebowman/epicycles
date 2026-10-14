/** Kepler propagation in the ecliptic J2000 frame (port of epicycles/kepler.py). */

export const GAUSS_K = 0.01720209895;
export const MU_SUN = GAUSS_K * GAUSS_K;

export type Vec3 = [number, number, number];

export interface Elements {
  a: number;
  e: number;
  i: number;
  Omega: number;
  omega: number;
  M0: number;
  epoch: number;
  mu: number;
  name: string;
}

/** Solve M = E - e sin E (Newton; same algorithm/tolerances as kepler.py). */
export function solveKepler(M: number, e: number, tol = 1e-15, maxiter = 100): number {
  if (!(e >= 0 && e < 1)) throw new RangeError("eccentricity must satisfy 0 <= e < 1");
  const twoPi = 2 * Math.PI;
  // numpy remainder semantics (result has the sign of the divisor)
  let x = (M + Math.PI) % twoPi;
  if (x < 0) x += twoPi;
  const Mr = x - Math.PI;
  let E: number;
  if (e > 0.8) E = Mr === 0 ? 0 : Math.PI * Math.sign(Mr);
  else E = Mr + e * Math.sin(Mr);
  for (let it = 0; it < maxiter; it++) {
    const f = E - e * Math.sin(E) - Mr;
    let dE = f / (1 - e * Math.cos(E));
    dE = Math.max(-1, Math.min(1, dE));
    E -= dE;
    if (Math.abs(dE) < tol || Math.abs(E - e * Math.sin(E) - Mr) < 1e-15) break;
  }
  return E + (M - Mr);
}

export function meanMotion(el: Elements): number {
  return Math.sqrt(el.mu / el.a ** 3);
}

/** Rows of R = R_z(Omega) R_x(i) R_z(omega) (perifocal -> ecliptic). */
export function rotation(el: Elements): [Vec3, Vec3, Vec3] {
  const cO = Math.cos(el.Omega), sO = Math.sin(el.Omega);
  const ci = Math.cos(el.i), si = Math.sin(el.i);
  const cw = Math.cos(el.omega), sw = Math.sin(el.omega);
  // R_x(i) R_z(omega)
  const A: [Vec3, Vec3, Vec3] = [
    [cw, -sw, 0],
    [ci * sw, ci * cw, -si],
    [si * sw, si * cw, ci],
  ];
  const R: [Vec3, Vec3, Vec3] = [
    [cO * A[0][0] - sO * A[1][0], cO * A[0][1] - sO * A[1][1], cO * A[0][2] - sO * A[1][2]],
    [sO * A[0][0] + cO * A[1][0], sO * A[0][1] + cO * A[1][1], sO * A[0][2] + cO * A[1][2]],
    [A[2][0], A[2][1], A[2][2]],
  ];
  return R;
}

/** Plane basis (P, Q, N) = columns of the rotation matrix. */
export function planeBasis(el: Elements): { P: Vec3; Q: Vec3; N: Vec3 } {
  const R = rotation(el);
  return {
    P: [R[0][0], R[1][0], R[2][0]],
    Q: [R[0][1], R[1][1], R[2][1]],
    N: [R[0][2], R[1][2], R[2][2]],
  };
}

export function checkElements(el: Elements): void {
  if (!(el.e >= 0 && el.e < 1)) throw new RangeError("eccentricity must satisfy 0 <= e < 1");
  if (!(el.a > 0)) throw new RangeError("semi-major axis must be positive");
}

/** Heliocentric ecliptic J2000 position [AU] at Julian date t. */
export function position(el: Elements, t: number): Vec3 {
  checkElements(el);
  const M = el.M0 + meanMotion(el) * (t - el.epoch);
  const E = solveKepler(M, el.e);
  const s = Math.sqrt(1 - el.e ** 2);
  const x = el.a * (Math.cos(E) - el.e);
  const y = el.a * s * Math.sin(E);
  const R = rotation(el);
  return [
    R[0][0] * x + R[0][1] * y,
    R[1][0] * x + R[1][1] * y,
    R[2][0] * x + R[2][1] * y,
  ];
}
