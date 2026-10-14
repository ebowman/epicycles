/** Load and validate epicycle-chain JSON (docs/epicycle-schema.md, version 1). */
import { Epicycle, Source } from "../math/expansion";
import { Vec3 } from "../math/kepler";

export interface ChainMeta {
  body: string | null;
  earth: string | null;
  N_body: number | null;
  N_earth: number | null;
  t_ref: number;
  units: { length: string; time: string; angle: string };
}

export interface LoadedChain {
  meta: ChainMeta;
  epicycles: Epicycle[];
  study?: unknown;
}

const UNIT_TOL = 1e-6;

function fail(path: string, msg: string): never {
  throw new Error(`Invalid epicycle chain: ${path}: ${msg}`);
}

function isObj(x: unknown): x is Record<string, unknown> {
  return typeof x === "object" && x !== null && !Array.isArray(x);
}

function num(x: unknown, path: string): number {
  if (typeof x !== "number" || !Number.isFinite(x)) fail(path, "must be a finite number");
  return x;
}

function strOrNull(x: unknown, path: string): string | null {
  if (x === null || x === undefined) return null;
  if (typeof x !== "string") fail(path, "must be a string or null");
  return x;
}

function intOrNull(x: unknown, path: string): number | null {
  if (x === null || x === undefined) return null;
  if (typeof x !== "number" || !Number.isInteger(x)) fail(path, "must be an integer or null");
  return x;
}

function vec3(x: unknown, path: string): Vec3 {
  if (!Array.isArray(x) || x.length !== 3) fail(path, "must be an array of length 3");
  return [num(x[0], `${path}[0]`), num(x[1], `${path}[1]`), num(x[2], `${path}[2]`)];
}

function unit(v: Vec3, path: string): Vec3 {
  const n = Math.hypot(v[0], v[1], v[2]);
  if (Math.abs(n - 1) > UNIT_TOL) fail(path, `must have unit length (got ${n})`);
  return v;
}

export function parseChain(json: unknown): LoadedChain {
  if (!isObj(json)) fail("(root)", "must be an object");
  if (json.schema_version !== 1) fail("schema_version", `must be 1 (got ${JSON.stringify(json.schema_version)})`);
  const m = json.meta;
  if (!isObj(m)) fail("meta", "must be an object");
  const u = m.units;
  if (!isObj(u)) fail("meta.units", "must be an object");
  if (u.length !== "AU") fail("meta.units.length", 'must be "AU"');
  if (u.time !== "day") fail("meta.units.time", 'must be "day"');
  if (u.angle !== "rad") fail("meta.units.angle", 'must be "rad"');
  const tRef = num(m.t_ref, "meta.t_ref");
  const meta: ChainMeta = {
    body: strOrNull(m.body, "meta.body"),
    earth: strOrNull(m.earth, "meta.earth"),
    N_body: intOrNull(m.N_body, "meta.N_body"),
    N_earth: intOrNull(m.N_earth, "meta.N_earth"),
    t_ref: tRef,
    units: { length: "AU", time: "day", angle: "rad" },
  };
  if (!Array.isArray(json.epicycles)) fail("epicycles", "must be an array");
  const epicycles: Epicycle[] = json.epicycles.map((raw: unknown, i: number): Epicycle => {
    const p = `epicycles[${i}]`;
    if (!isObj(raw)) fail(p, "must be an object");
    const radius = num(raw.radius, `${p}.radius`);
    if (radius < 0) fail(`${p}.radius`, "must be >= 0");
    const src = raw.source;
    if (src !== "body" && src !== "earth" && src !== "naff") {
      fail(`${p}.source`, `must be "body", "earth" or "naff" (got ${JSON.stringify(src)})`);
    }
    return {
      radius,
      omega: num(raw.omega, `${p}.omega`),
      phase: num(raw.phase, `${p}.phase`),
      normal: unit(vec3(raw.normal, `${p}.normal`), `${p}.normal`),
      u: unit(vec3(raw.u, `${p}.u`), `${p}.u`),
      v: unit(vec3(raw.v, `${p}.v`), `${p}.v`),
      tRef,
      source: src as Source,
      k: intOrNull(raw.k, `${p}.k`),
    };
  });
  for (let i = 1; i < epicycles.length; i++) {
    if (epicycles[i]!.radius > epicycles[i - 1]!.radius) {
      fail(`epicycles[${i}].radius`, "epicycles must be sorted by radius descending");
    }
  }
  const out: LoadedChain = { meta, epicycles };
  if ("study" in json) out.study = json.study;
  return out;
}

export async function loadChain(url: string): Promise<LoadedChain> {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`Failed to fetch ${url}: HTTP ${resp.status}`);
  return parseChain(await resp.json());
}
