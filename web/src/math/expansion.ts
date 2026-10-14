/** Analytic epicycle expansion of Kepler orbits (port of epicycles/expansion.py). */
import { besselJArray } from "./bessel";
import { Elements, Vec3, checkElements, meanMotion, planeBasis } from "./kepler";

export type Source = "body" | "earth" | "naff";

export interface Epicycle {
  radius: number;
  omega: number;
  phase: number;
  normal: Vec3;
  u: Vec3;
  v: Vec3;
  tRef: number;
  source: Source;
  k: number | null;
}

export interface ComplexArray {
  re: Float64Array;
  im: Float64Array;
}

/**
 * Fourier coefficients c_k, k=-N..N of z/a = sum c_k exp(ikM). Index k+N.
 * (All coefficients are real, so im is identically zero.)
 */
export function keplerFourier(e: number, N: number): ComplexArray {
  if (!(e >= 0 && e < 1)) throw new RangeError("eccentricity must satisfy 0 <= e < 1");
  if (!Number.isInteger(N) || N < 0) throw new RangeError("N must be a non-negative integer");
  const re = new Float64Array(2 * N + 1);
  const im = new Float64Array(2 * N + 1);
  if (e === 0) {
    if (N >= 1) re[N + 1] = 1;
    return { re, im };
  }
  const s = Math.sqrt(1 - e * e);
  re[N] = -1.5 * e;
  for (let k = 1; k <= N; k++) {
    const J = besselJArray(k + 1, k * e);
    const jk = J[k]!, jk1 = J[k + 1]!;
    re[N + k] = (((1 + s) / e) * jk - jk1) / k;
    re[N - k] = ((e / (1 + s)) * jk - jk1) / k;
  }
  return { re, im };
}

/** Truncated (|k|<=N) epicycle expansion of the heliocentric orbit, sorted by radius desc. */
export function toEpicycles(el: Elements, N: number, tRef?: number,
                            source: Source = "body"): Epicycle[] {
  checkElements(el);
  const t0 = tRef ?? el.epoch;
  const c = keplerFourier(el.e, N);
  const { P, Q, N: Nn } = planeBasis(el);
  const n = meanMotion(el);
  const Mref = el.M0 + n * (t0 - el.epoch);
  const out: Epicycle[] = [];
  for (let idx = 0; idx < c.re.length; idx++) {
    const k = idx - N;
    const ck = c.re[idx]!;
    if (ck === 0) continue;
    out.push({
      radius: el.a * Math.abs(ck),
      omega: k * n,
      phase: k * Mref + (ck < 0 ? Math.PI : 0),
      normal: [...Nn] as Vec3,
      u: [...P] as Vec3,
      v: [...Q] as Vec3,
      tRef: t0,
      source,
      k,
    });
  }
  out.sort((a, b) => b.radius - a.radius);
  return out;
}

/** Epicycles of r_body - r_earth (Earth terms negated: phase + pi). */
export function geocentricChain(body: Elements, earth: Elements, Nbody: number, Nearth: number,
                                tRef?: number): Epicycle[] {
  const t0 = tRef ?? body.epoch;
  const eps = toEpicycles(earth, Nearth, t0, "earth").map((ep) => ({ ...ep, phase: ep.phase + Math.PI }));
  const chain = [...toEpicycles(body, Nbody, t0, "body"), ...eps];
  chain.sort((a, b) => b.radius - a.radius);
  return chain;
}

/** Sum of all circles at Julian date t [AU]. */
export function evaluateChain(chain: readonly Epicycle[], t: number): Vec3 {
  let x = 0, y = 0, z = 0;
  for (const c of chain) {
    const th = c.phase + c.omega * (t - c.tRef);
    const a = c.radius * Math.cos(th), b = c.radius * Math.sin(th);
    x += a * c.u[0] + b * c.v[0];
    y += a * c.u[1] + b * c.v[1];
    z += a * c.u[2] + b * c.v[2];
  }
  return [x, y, z];
}

/** Cumulative arm tips, flat Float64Array of 3*(n+1); point 0 = origin. */
export function chainArms(chain: readonly Epicycle[], t: number): Float64Array {
  const out = new Float64Array(3 * (chain.length + 1));
  let x = 0, y = 0, z = 0;
  chain.forEach((c, j) => {
    const th = c.phase + c.omega * (t - c.tRef);
    const a = c.radius * Math.cos(th), b = c.radius * Math.sin(th);
    x += a * c.u[0] + b * c.v[0];
    y += a * c.u[1] + b * c.v[1];
    z += a * c.u[2] + b * c.v[2];
    out[3 * (j + 1)] = x;
    out[3 * (j + 1) + 1] = y;
    out[3 * (j + 1) + 2] = z;
  });
  return out;
}
