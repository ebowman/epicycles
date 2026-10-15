/** App state shared by both views: the clock plus the epicycle count N (per body). */
import { Clock } from "./clock";

export const N_MIN = 1;
export const N_MAX = 60;

export class AppState {
  readonly clock: Clock;
  N: number;
  constructor(clock = new Clock(), N = 20) {
    this.clock = clock;
    this.N = AppState.clampN(N);
  }
  /** Upper bound for N in the current mode (per-body, see presets.nMaxFor). Default N_MAX. */
  nMax = N_MAX;
  static clampN(n: number, max = N_MAX): number {
    if (!Number.isFinite(n)) return N_MIN;
    return Math.min(max, Math.max(N_MIN, Math.round(n)));
  }
  /** Change the upper bound (and re-clamp N). Returns true if N changed. */
  setNMax(max: number): boolean {
    this.nMax = Math.max(N_MIN, Math.round(max));
    return this.setN(this.N);
  }
  /** Change N. Never touches the clock. Returns true if N actually changed. */
  setN(n: number): boolean {
    const m = AppState.clampN(n, this.nMax);
    if (m === this.N) return false;
    this.N = m;
    return true;
  }
}
