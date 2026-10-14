/**
 * Bessel functions of the first kind, integer order, via Miller backward
 * recurrence normalised with J_0 + 2 sum J_{2k} = 1. Valid for n >= 0, x >= 0
 * (tested to x ~ 700, n ~ 700).
 */

/** Returns [J_0(x), ..., J_nmax(x)]. */
export function besselJArray(nmax: number, x: number): Float64Array {
  if (!Number.isInteger(nmax) || nmax < 0) throw new RangeError("order must be a non-negative integer");
  if (!(x >= 0) || !Number.isFinite(x)) throw new RangeError("x must be finite and >= 0");
  const out = new Float64Array(nmax + 1);
  if (x === 0) {
    out[0] = 1;
    return out;
  }
  const big = Math.max(nmax, x);
  let m = Math.ceil(big + 40 + 10 * Math.cbrt(big));
  if (m % 2) m += 1;
  const tmp = new Float64Array(m + 2);
  const RESCALE = 1e250;
  tmp[m + 1] = 0;
  tmp[m] = 1e-300;
  let norm = 0; // sum for even orders j>=2 (to be doubled)
  const twoOverX = 2 / x;
  for (let j = m; j >= 1; j--) {
    const prev = j * twoOverX * tmp[j]! - tmp[j + 1]!; // J_{j-1}
    tmp[j - 1] = prev;
    if (Math.abs(prev) > RESCALE) {
      for (let q = j - 1; q <= m + 1; q++) tmp[q] *= 1 / RESCALE;
      norm *= 1 / RESCALE;
    }
    if ((j - 1) % 2 === 0 && j - 1 >= 2) norm += tmp[j - 1]!;
  }
  const total = tmp[0]! + 2 * norm;
  for (let j = 0; j <= nmax; j++) out[j] = tmp[j]! / total;
  return out;
}

export function besselJ(n: number, x: number): number {
  return besselJArray(n, x)[n]!;
}

/** J_n'(x) = (J_{n-1} - J_{n+1}) / 2 (n >= 1), J_0' = -J_1. */
export function besselJp(n: number, x: number): number {
  const a = besselJArray(n + 1, x);
  return n === 0 ? -a[1]! : 0.5 * (a[n - 1]! - a[n + 1]!);
}
