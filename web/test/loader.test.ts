import { describe, expect, it } from "vitest";
import { readFileSync, readdirSync } from "node:fs";
import { parseChain } from "../src/data/loader";
import { BODIES } from "../src/math/bodies";
import { evaluateChain } from "../src/math/expansion";
import { position } from "../src/math/kepler";

const dataDir = new URL("../public/data/", import.meta.url);
const read = (f: string) => JSON.parse(readFileSync(new URL(f, dataDir), "utf8"));
const clone = () => structuredClone(read("mars_N5.json"));

describe("parseChain valid", () => {
  it("parses every file in public/data (except index)", () => {
    for (const f of readdirSync(dataDir)) {
      if (f === "index.json" || !f.endsWith(".json")) continue;
      const c = parseChain(read(f));
      expect(c.epicycles.length).toBeGreaterThan(0);
    }
  });
  it("naff chain has naff sources and study key", () => {
    const c = parseChain(read("perturbed_ceres.json"));
    expect(c.epicycles.every((e) => e.source === "naff")).toBe(true);
    expect(c.study).toBeDefined();
  });
  it("mars_N60 chain matches geocentric Kepler truth < 1e-9 AU", () => {
    const c = parseChain(read("mars_N60.json"));
    let worst = 0;
    for (const dt of [-8000, -1234.5, 0, 55.5, 687, 4000, 9000]) {
      const t = c.meta.t_ref + dt;
      const r = evaluateChain(c.epicycles, t);
      const a = position(BODIES.MARS!, t), b = position(BODIES.EARTH!, t);
      for (let i = 0; i < 3; i++) worst = Math.max(worst, Math.abs(r[i]! - (a[i]! - b[i]!)));
    }
    console.log("mars_N60 vs truth max err (AU):", worst);
    expect(worst).toBeLessThan(1e-9);
  });
});

describe("parseChain malformed", () => {
  const bad = (mut: (d: any) => void, re: RegExp) => {
    const d = clone();
    mut(d);
    expect(() => parseChain(d)).toThrow(re);
  };
  it("rejects non-object", () => expect(() => parseChain(5)).toThrow(/\(root\)/));
  it("bad schema_version", () => bad((d) => (d.schema_version = 2), /schema_version/));
  it("missing meta", () => bad((d) => delete d.meta, /meta/));
  it("bad units", () => bad((d) => (d.meta.units.length = "km"), /meta\.units\.length/));
  it("missing t_ref", () => bad((d) => delete d.meta.t_ref, /meta\.t_ref/));
  it("epicycles not array", () => bad((d) => (d.epicycles = {}), /epicycles/));
  it("bad source", () => bad((d) => (d.epicycles[3].source = "moon"), /epicycles\[3\]\.source/));
  it("negative radius", () => bad((d) => (d.epicycles[0].radius = -1), /epicycles\[0\]\.radius/));
  it("string omega", () => bad((d) => (d.epicycles[2].omega = "x"), /epicycles\[2\]\.omega/));
  it("NaN phase (null)", () => bad((d) => (d.epicycles[1].phase = null), /epicycles\[1\]\.phase/));
  it("normal wrong length", () => bad((d) => (d.epicycles[0].normal = [0, 0]), /epicycles\[0\]\.normal/));
  it("u not unit", () => bad((d) => (d.epicycles[4].u = [2, 0, 0]), /epicycles\[4\]\.u.*unit/));
  it("non-integer k", () => bad((d) => (d.epicycles[0].k = 1.5), /epicycles\[0\]\.k/));
  it("unsorted radii", () =>
    bad((d) => {
      d.epicycles.reverse();
    }, /sorted/));
});
