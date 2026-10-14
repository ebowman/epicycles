/** Shared simulation clock. Pure TS, no DOM. t is a Julian date [days]. */
export const MIN_SPEED = 1;
export const MAX_SPEED = 200;

export class Clock {
  t: number;
  playing: boolean;
  speed: number; // simulated days per real second

  constructor(t = 2451545.0, speed = 30, playing = true) {
    this.t = t;
    this.speed = Clock.clampSpeed(speed);
    this.playing = playing;
  }

  static clampSpeed(s: number): number {
    if (!Number.isFinite(s)) return MIN_SPEED;
    return Math.min(MAX_SPEED, Math.max(MIN_SPEED, s));
  }

  /** Advance by dtReal real seconds (ignored when paused or dt is not a finite non-negative number). */
  step(dtReal: number): void {
    if (!this.playing || !Number.isFinite(dtReal) || dtReal < 0) return;
    this.t += this.speed * dtReal;
  }
  setSpeed(s: number): void { this.speed = Clock.clampSpeed(s); }
  play(): void { this.playing = true; }
  pause(): void { this.playing = false; }
  toggle(): boolean { this.playing = !this.playing; return this.playing; }
  scrub(t: number): void { if (Number.isFinite(t)) this.t = t; }
}
