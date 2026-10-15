import { describe, expect, it } from "vitest";
import { readFileSync } from "node:fs";
import { parseChain } from "../src/data/loader";
import { EARTH } from "../src/math/bodies";
import { evaluateChain, geocentricChain } from "../src/math/expansion";
import { position } from "../src/math/kepler";
import {
  BODY_ELEMENTS, BODY_IDS, NAFF_WINDOW, circleCountFor, fourierN, nMaxFor, parseParams, presetDefaults,
  ptolemaicChain, selectTopK, studyRows,
} from "../src/sim/presets";
import { AppState } from "../src/sim/state";
import { FORMULAS } from "../src/ui/info";
import katex from "katex";

const dataDir = new URL("../public/data/", import.meta.url);
const naffJson = () => JSON.parse(readFileSync(new URL("perturbed_ceres.json", dataDir), "utf8"));

describe("preset loader smoke", () => {
  for (const p of ["ptolemy", "fourier", "ceres-naff"] as const) {
    it(`preset ${p} resolves to a valid state`, () => {
      const s = parseParams(`?preset=${p}`);
      expect(BODY_IDS).toContain(s.body);
      expect(s.mode).toBe(s.body === "naff" ? "naff" : "kepler");
      expect(Number.isInteger(s.N)).toBe(true);
      expect(s.N).toBeGreaterThanOrEqual(1);
      expect(s.N).toBeLessThanOrEqual(nMaxFor(s.body));
      expect(Number.isFinite(s.t)).toBe(true);
      if (s.mode === "naff") {
        expect(s.t).toBeGreaterThanOrEqual(NAFF_WINDOW[0]);
        expect(s.t).toBeLessThanOrEqual(NAFF_WINDOW[1]);
      } else {
        const chain = geocentricChain(BODY_ELEMENTS[s.body as "mars"], EARTH, s.N, s.N);
        expect(chain.length).toBeGreaterThan(0);
      }
      expect(presetDefaults(p).body).toBe(s.body);
    });
  }
  it("ptolemy is Mars N=1 near the 2003 opposition with error shown; fourier animates", () => {
    const s = parseParams("?preset=ptolemy");
    expect([s.body, s.N, s.showError]).toEqual(["mars", 1, true]);
    expect(Math.abs(s.t - 2452880)).toBeLessThan(30);
    expect(parseParams("?preset=fourier").fourier).toBe(true);
    expect(fourierN(0)).toBe(1);
    expect(fourierN(1e6)).toBe(20);
    expect(fourierN(2)).toBeGreaterThan(1);
  });
});

describe("Ptolemaic model", () => {
  it("is 2 circles (k=+1 body + k=+1 Earth) with ~33 deg peak longitude error over 60 yr", () => {
    const c = ptolemaicChain(geocentricChain(BODY_ELEMENTS.mars, EARTH, 1, 1));
    expect(c.map((e) => e.source).sort()).toEqual(["body", "earth"]);
    expect(parseParams("?preset=ptolemy").ptolemaic).toBe(true);
    let mx = 0;
    for (let t = 2451545; t < 2451545 + 60 * 365.25; t += 3) {
      const a = evaluateChain(c, t), pm = position(BODY_ELEMENTS.mars, t), pe = position(EARTH, t);
      let d = Math.atan2(a[1], a[0]) - Math.atan2(pm[1] - pe[1], pm[0] - pe[0]);
      d = Math.abs(Math.atan2(Math.sin(d), Math.cos(d))) * 180 / Math.PI;
      mx = Math.max(mx, d);
    }
    expect(mx).toBeGreaterThan(25);
    expect(mx).toBeLessThan(40);
  });
});

describe("URL param parsing", () => {
  it("defaults", () => {
    const s = parseParams("");
    expect([s.body, s.N, s.t, s.paused, s.preset, s.zoom, s.info]).toEqual(["mars", 20, 2452880, false, null, null, true]);
  });
  it("explicit params, clamping and junk", () => {
    const s = parseParams("?body=ceres&N=7&t=2451545&paused=1&speed=9999&zoom=opp&info=0&circles=0");
    expect([s.body, s.N, s.t, s.paused, s.speed, s.zoom, s.info, s.showCircles]).toEqual(["ceres", 7, 2451545, true, 200, "opposition", false, false]);
    const j = parseParams("?body=pluto&N=abc&t=&preset=zzz");
    expect([j.body, j.N, j.t, j.preset]).toEqual(["mars", 20, 2452880, null]);
    expect(parseParams("?body=mars&N=1000").N).toBe(60);
    expect(parseParams("?body=phaethon&N=1000").N).toBe(300);
    expect(parseParams("?body=mars&N=-3").N).toBe(1);
  });
  it("naff body clamps t into the study window", () => {
    const s = parseParams("?body=naff&t=2452880");
    expect(s.t).toBe(NAFF_WINDOW[0]);
    expect(parseParams("?body=naff&N=500").N).toBe(200);
  });
  it("explicit params override a preset", () => {
    const s = parseParams("?preset=ptolemy&N=5&body=eros");
    expect([s.body, s.N]).toEqual(["eros", 5]);
  });
});

describe("N range per body", () => {
  it("limits", () => {
    expect(nMaxFor("mars")).toBe(60);
    expect(nMaxFor("ceres")).toBe(60);
    expect(nMaxFor("eros")).toBe(60);
    expect(nMaxFor("phaethon")).toBe(300);
    expect(nMaxFor("naff")).toBe(200);
    expect(circleCountFor("phaethon", 300)).toBe(1202);
    expect(circleCountFor("naff", 200)).toBe(401);
  });
  it("AppState honours per-body max without breaking the default", () => {
    const s = new AppState();
    expect(s.nMax).toBe(60);
    s.setNMax(300); expect(s.setN(300)).toBe(true); expect(s.N).toBe(300);
    s.setNMax(60); expect(s.N).toBe(60);
  });
  it("phaethon N=300 chain is finite (up to 1202 circles)", () => {
    const c = geocentricChain(BODY_ELEMENTS.phaethon, EARTH, 300, 300);
    expect(c.length).toBeGreaterThan(900);
    expect(c.length).toBeLessThanOrEqual(1202);
    expect(c.every((e) => Number.isFinite(e.radius) && Number.isFinite(e.phase))).toBe(true);
  });
});

describe("NAFF", () => {
  const chain = parseChain(naffJson());
  it("perturbed_ceres.json loads with 401 circles", () => {
    expect(chain.epicycles.length).toBe(401);
    expect(chain.meta.t_ref).toBe(2461200.5);
  });
  it("selectTopK keeps the 2K+1 largest by radius, sorted descending", () => {
    const top = selectTopK(chain.epicycles, 10);
    expect(top.length).toBe(21);
    for (let i = 1; i < top.length; i++) expect(top[i]!.radius).toBeLessThanOrEqual(top[i - 1]!.radius);
    const cut = top[top.length - 1]!.radius;
    const rest = chain.epicycles.filter((e) => !top.includes(e));
    expect(rest.every((e) => e.radius <= cut)).toBe(true);
  });
  it("selectTopK clamps and handles unsorted input", () => {
    expect(selectTopK(chain.epicycles, 1000).length).toBe(401);
    expect(selectTopK(chain.epicycles, 0).length).toBe(3);
    expect(selectTopK(chain.epicycles, NaN).length).toBe(3);
    const shuffled = [...chain.epicycles].reverse();
    expect(selectTopK(shuffled, 5)[0]!.radius).toBe(chain.epicycles[0]!.radius);
  });
  it("study table rows", () => {
    const rows = studyRows(chain.study);
    expect(rows.find((r) => r.label === "NAFF K=200")?.circles).toBe("401");
    expect(rows.some((r) => r.label.startsWith("Kepler"))).toBe(true);
    expect(studyRows(null)).toEqual([]);
  });
});

describe("KaTeX panel formulas", () => {
  it("all render without error", () => {
    for (const tex of Object.values(FORMULAS)) {
      expect(() => katex.renderToString(tex, { throwOnError: true })).not.toThrow();
    }
  });
});

describe("drawn circle count (info panel / HUD source of truth = chain length)", () => {
  it("analytic chain length equals circleCountFor; Ptolemaic N=1 draws 2 (not 6)", () => {
    for (const n of [1, 3, 7]) {
      expect(geocentricChain(BODY_ELEMENTS.mars, EARTH, n, n).length).toBe(circleCountFor("mars", n));
    }
    const full = geocentricChain(BODY_ELEMENTS.mars, EARTH, 1, 1);
    expect(circleCountFor("mars", 1)).toBe(6);
    expect(ptolemaicChain(full).length).toBe(2);
  });
});
