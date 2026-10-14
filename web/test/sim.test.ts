import { describe, expect, it } from "vitest";
import { Clock } from "../src/sim/clock";
import { AppState } from "../src/sim/state";
import { RingTrail } from "../src/sim/trail";
import { chainArmsInto, chainTipInto, circleMatrix, computeError, fitDistance, jdToDate, ErrorReadout } from "../src/sim/geometry";
import { EARTH, MARS, J2000 } from "../src/math/bodies";
import { chainArms, evaluateChain, geocentricChain } from "../src/math/expansion";
import { position } from "../src/math/kepler";

describe("clock", () => {
  it("steps only while playing", () => {
    const c = new Clock(1000, 10, true);
    c.step(0.5);
    expect(c.t).toBeCloseTo(1005);
    c.pause();
    c.step(1);
    expect(c.t).toBeCloseTo(1005);
    c.play();
    c.step(1);
    expect(c.t).toBeCloseTo(1015);
    expect(c.toggle()).toBe(false);
  });
  it("clamps speed and ignores bad input", () => {
    const c = new Clock(0, 10, true);
    c.setSpeed(1e6); expect(c.speed).toBe(200);
    c.setSpeed(-3); expect(c.speed).toBe(1);
    c.setSpeed(NaN); expect(c.speed).toBe(1);
    c.step(NaN); c.step(-1); expect(c.t).toBe(0);
    c.scrub(NaN); expect(c.t).toBe(0);
    c.scrub(42); expect(c.t).toBe(42);
  });
  it("N change does not touch the clock", () => {
    const s = new AppState(new Clock(2452880, 25, true), 10);
    s.clock.step(2);
    const t = s.clock.t;
    expect(s.setN(60)).toBe(true);
    expect(s.clock.t).toBe(t);
    expect(s.clock.playing).toBe(true);
    expect(s.clock.speed).toBe(25);
    expect(s.setN(60)).toBe(false);
    s.setN(500); expect(s.N).toBe(60);
    s.setN(0); expect(s.N).toBe(1);
    s.setN(NaN); expect(s.N).toBe(1);
  });
});

describe("geometry", () => {
  const chain = geocentricChain(MARS, EARTH, 12, 12);
  const t = J2000 + 1234.5;
  it("chainArmsInto matches chainArms and tip matches evaluateChain", () => {
    const ref = chainArms(chain, t);
    const out = new Float64Array(ref.length);
    chainArmsInto(chain, t, out);
    for (let i = 0; i < ref.length; i++) expect(out[i]).toBe(ref[i]);
    const tip = new Float64Array(3);
    chainTipInto(chain, t, tip);
    const e = evaluateChain(chain, t);
    for (let i = 0; i < 3; i++) expect(tip[i]).toBeCloseTo(e[i]!, 14);
  });
  it("circle matrix maps unit circle onto the epicycle circle at the arm point", () => {
    const arms = chainArms(chain, t);
    const m = new Float64Array(16);
    for (const j of [0, 3, chain.length - 1]) {
      const c = chain[j]!;
      circleMatrix(c, arms[3 * j]!, arms[3 * j + 1]!, arms[3 * j + 2]!, m);
      // unit-circle point at the current angle lands on next arm point
      const th = c.phase + c.omega * (t - c.tRef);
      const x = Math.cos(th), y = Math.sin(th);
      for (let k = 0; k < 3; k++) {
        const p = m[k]! * x + m[4 + k]! * y + m[12 + k]!;
        expect(p).toBeCloseTo(arms[3 * (j + 1) + k]!, 12);
      }
    }
  });
  it("error readout is ~0 for N=60 Mars and sane for N=1", () => {
    const big = geocentricChain(MARS, EARTH, 60, 60);
    expect(big.length).toBeLessThanOrEqual(242);
    const r: ErrorReadout = { au: 0, lonArcsec: 0 };
    const tip = new Float64Array(3);
    for (const dt of [0, 500, -3000, 6000]) {
      const tt = J2000 + dt;
      chainTipInto(big, tt, tip);
      const pm = position(MARS, tt), pe = position(EARTH, tt);
      computeError(tip, [pm[0] - pe[0], pm[1] - pe[1], pm[2] - pe[2]], r);
      expect(r.au).toBeLessThan(1e-9);
      expect(Math.abs(r.lonArcsec)).toBeLessThan(1);
    }
    const small = geocentricChain(MARS, EARTH, 1, 1);
    chainTipInto(small, J2000, tip);
    const pm = position(MARS, J2000), pe = position(EARTH, J2000);
    computeError(tip, [pm[0] - pe[0], pm[1] - pe[1], pm[2] - pe[2]], r);
    expect(r.au).toBeGreaterThan(1e-3);
  });
  it("longitude error wraps across +-180 deg and is exact for a known offset", () => {
    const r: ErrorReadout = { au: 0, lonArcsec: 0 };
    const a = 1e-3; // rad
    computeError([Math.cos(Math.PI - a), Math.sin(Math.PI - a), 0], [-1, 0, 0], r);
    expect(Math.abs(r.lonArcsec)).toBeCloseTo(a * 206264.806247, 3);
  });
  it("jdToDate", () => { expect(jdToDate(2451545)).toBe("2000-01-01"); });
});

describe("RingTrail", () => {
  it("keeps newest points in chronological order after wrap", () => {
    const tr = new RingTrail(3);
    for (let i = 1; i <= 5; i++) tr.push(i, 0, 0);
    const out = new Float32Array(9);
    expect(tr.copyOrdered(out)).toBe(3);
    expect([out[0], out[3], out[6]]).toEqual([3, 4, 5]);
    tr.clear(); expect(tr.count).toBe(0);
  });
});

describe("fitDistance", () => {
  it("keeps a wide box inside a narrow viewport, with margin", () => {
    const fov = 40, hw = 2, hh = 2.3, aspect = 0.72;
    const d = fitDistance(hw, hh, fov, aspect);
    const tanH = Math.tan((fov * Math.PI) / 360);
    expect(d * tanH * aspect).toBeGreaterThanOrEqual(1.25 * hw);
    expect(d * tanH).toBeGreaterThanOrEqual(1.25 * hh);
  });
  it("is height-limited in a wide viewport and clamps tiny aspects", () => {
    const tanH = Math.tan((40 * Math.PI) / 360);
    expect(fitDistance(1, 1, 40, 3)).toBeCloseTo(1.3 / tanH, 9);
    expect(fitDistance(1, 1, 40, 0)).toBeCloseTo(1.3 / (tanH * 0.2), 9);
  });
});
