/** Allocation-free geometry helpers for drawing the epicycle chain. */
import type { Epicycle } from "../math/expansion";

type Writable = { [i: number]: number };

/** Same values as chainArms(), written into `out` (length >= 3*(chain.length+1)). */
export function chainArmsInto(chain: readonly Epicycle[], t: number, out: Writable): void {
  let x = 0, y = 0, z = 0;
  out[0] = 0; out[1] = 0; out[2] = 0;
  for (let j = 0; j < chain.length; j++) {
    const c = chain[j]!;
    const th = c.phase + c.omega * (t - c.tRef);
    const a = c.radius * Math.cos(th), b = c.radius * Math.sin(th);
    x += a * c.u[0] + b * c.v[0];
    y += a * c.u[1] + b * c.v[1];
    z += a * c.u[2] + b * c.v[2];
    out[3 * (j + 1)] = x; out[3 * (j + 1) + 1] = y; out[3 * (j + 1) + 2] = z;
  }
}

/** Chain sum at t written to out[0..2] (same as evaluateChain, no allocation). */
export function chainTipInto(chain: readonly Epicycle[], t: number, out: Writable): void {
  let x = 0, y = 0, z = 0;
  for (let j = 0; j < chain.length; j++) {
    const c = chain[j]!;
    const th = c.phase + c.omega * (t - c.tRef);
    const a = c.radius * Math.cos(th), b = c.radius * Math.sin(th);
    x += a * c.u[0] + b * c.v[0];
    y += a * c.u[1] + b * c.v[1];
    z += a * c.u[2] + b * c.v[2];
  }
  out[0] = x; out[1] = y; out[2] = z;
}

/**
 * Column-major 4x4 transform mapping the unit circle (cos s, sin s, 0) to the epicycle circle
 * centred at (cx,cy,cz): columns radius*u, radius*v, normal, centre.
 */
export function circleMatrix(c: Epicycle, cx: number, cy: number, cz: number, m: Writable): void {
  const r = c.radius;
  m[0] = r * c.u[0]; m[1] = r * c.u[1]; m[2] = r * c.u[2]; m[3] = 0;
  m[4] = r * c.v[0]; m[5] = r * c.v[1]; m[6] = r * c.v[2]; m[7] = 0;
  m[8] = c.normal[0]; m[9] = c.normal[1]; m[10] = c.normal[2]; m[11] = 0;
  m[12] = cx; m[13] = cy; m[14] = cz; m[15] = 1;
}

export interface ErrorReadout { au: number; lonArcsec: number }

const ARCSEC_PER_RAD = (180 / Math.PI) * 3600;

/** |chain tip - truth| in AU and the wrapped ecliptic-longitude difference in arcsec (signed). */
export function computeError(tip: ArrayLike<number>, truth: ArrayLike<number>, out: ErrorReadout): ErrorReadout {
  const dx = tip[0]! - truth[0]!, dy = tip[1]! - truth[1]!, dz = tip[2]! - truth[2]!;
  out.au = Math.sqrt(dx * dx + dy * dy + dz * dz);
  let d = Math.atan2(tip[1]!, tip[0]!) - Math.atan2(truth[1]!, truth[0]!);
  d = Math.atan2(Math.sin(d), Math.cos(d));
  out.lonArcsec = d * ARCSEC_PER_RAD;
  return out;
}

export function formatLonError(arcsec: number): string {
  const a = Math.abs(arcsec);
  if (a >= 3600) return `${(a / 3600).toFixed(2)}°`;
  if (a >= 100) return `${a.toFixed(0)}″`;
  return `${a.toPrecision(3)}″`;
}

export function formatAU(x: number): string {
  return x === 0 ? "0 AU" : `${x.toExponential(2)} AU`;
}

export function jdToDate(jd: number): string {
  const ms = (jd - 2440587.5) * 86400000;
  const d = new Date(ms);
  return Number.isNaN(d.getTime()) ? "n/a" : d.toISOString().slice(0, 10);
}
