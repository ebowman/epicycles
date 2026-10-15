import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { besselJ, besselJp } from "../src/math/bessel";
import { BODIES } from "../src/math/bodies";
import { chainArms, evaluateChain, geocentricChain, keplerFourier } from "../src/math/expansion";
import { position, solveKepler } from "../src/math/kepler";
import { parseChain } from "../src/data/loader";

const g = (n: string) => JSON.parse(readFileSync(new URL(`./golden/${n}`, import.meta.url), "utf8"));
const maxErr: Record<string, number> = {};
function close(name: string, a: number, b: number, rel = 1e-9, abs = 1e-12): void {
  const err = Math.abs(a - b);
  maxErr[name] = Math.max(maxErr[name] ?? 0, err);
  expect(err, `${name}: ${a} vs ${b}`).toBeLessThanOrEqual(abs + rel * Math.abs(b));
}

describe("kepler", () => {
  it("solveKepler matches golden and residual", () => {
    for (const s of g("kepler.json")) {
      const E = solveKepler(s.M, s.e);
      close("solveKepler", E, s.E, 1e-12, 1e-12);
      // E - M is periodic in M: residual of the unwrapped-branch solution against M itself
      expect(Math.abs(E - s.e * Math.sin(E) - s.M)).toBeLessThan(1e-13 * Math.max(1, Math.abs(s.M)));
    }
  });
  it("rejects bad eccentricity", () => {
    expect(() => solveKepler(1, 1)).toThrow();
    expect(() => solveKepler(1, -0.1)).toThrow();
  });
  it("positions match golden", () => {
    for (const s of g("positions.json")) {
      const r = position(BODIES[s.body]!, s.t);
      for (let i = 0; i < 3; i++) close("position", r[i]!, s.r[i]);
    }
  });
});

describe("bessel", () => {
  it("J_n and J_n' match scipy to 1e-12 absolute", () => {
    const rows = g("bessel.json");
    expect(rows.length).toBeGreaterThan(100);
    for (const s of rows) {
      close("besselJ", besselJ(s.n, s.x), s.J, 0, 1e-13);
      close("besselJp", besselJp(s.n, s.x), s.Jp, 0, 1e-13);
    }
  });
  it("edge cases", () => {
    expect(besselJ(0, 0)).toBe(1);
    expect(besselJ(3, 0)).toBe(0);
    expect(() => besselJ(1, -1)).toThrow();
  });
});

describe("expansion", () => {
  it("keplerFourier matches golden", () => {
    for (const s of g("kepler_fourier.json")) {
      const c = keplerFourier(s.e, s.N);
      expect(c.re.length).toBe(s.c.length);
      for (let i = 0; i < s.c.length; i++) close("keplerFourier", c.re[i]!, s.c[i], 1e-9, 1e-12);
    }
  });
  it("geocentric chain evaluation matches Python and truth", () => {
    for (const s of g("geocentric.json")) {
      const chain = geocentricChain(BODIES[s.body]!, BODIES.EARTH!, s.Nb, s.Ne);
      const r = evaluateChain(chain, s.t);
      const arms = chainArms(chain, s.t);
      const last = arms.length - 3;
      for (let i = 0; i < 3; i++) {
        close("chain", r[i]!, s.r[i]);
        close("armsLast", arms[last + i]!, s.arms_last[i]);
        close("armsMid", arms[3 * Math.floor(chain.length / 2) + i]!, s.arms_mid[i], 1e-9, 1e-11);
      }
      if (s.Nb >= 20 && s.body === "MARS") {
        const a = position(BODIES[s.body]!, s.t), b = position(BODIES.EARTH!, s.t);
        for (let i = 0; i < 3; i++) close("truthVsChain", r[i]!, a[i]! - b[i]!, 0, 1e-6);
      }
    }
  });
  it("full Mars N=5 chain matches field by field", () => {
    const d = g("mars_N5_chain.json");
    const chain = geocentricChain(BODIES.MARS!, BODIES.EARTH!, 5, 5);
    expect(chain.length).toBe(d.epicycles.length);
    for (let i = 0; i < chain.length; i++) {
      const c = chain[i]!, p = d.epicycles[i];
      expect(c.source).toBe(p.source);
      expect(c.k).toBe(p.k);
      close("chain.radius", c.radius, p.radius);
      close("chain.omega", c.omega, p.omega);
      close("chain.phase", c.phase, p.phase, 1e-9, 1e-9);
      for (let j = 0; j < 3; j++) {
        close("chain.normal", c.normal[j]!, p.normal[j]);
        close("chain.u", c.u[j]!, p.u[j]);
        close("chain.v", c.v[j]!, p.v[j]);
      }
    }
    // loader round trip of the same structure
    expect(parseChain(d).epicycles.length).toBe(chain.length);
  });
  it("reports max discrepancies", () => {
    console.log("MAX GOLDEN DISCREPANCIES", JSON.stringify(maxErr));
  });
});
