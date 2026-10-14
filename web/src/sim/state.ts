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
  static clampN(n: number): number {
    if (!Number.isFinite(n)) return N_MIN;
    return Math.min(N_MAX, Math.max(N_MIN, Math.round(n)));
  }
  /** Change N. Never touches the clock. Returns true if N actually changed. */
  setN(n: number): boolean {
    const m = AppState.clampN(n);
    if (m === this.N) return false;
    this.N = m;
    return true;
  }
}
