/** Pure (DOM-free) body/preset/URL-parameter logic for the SPA. */
import type { Elements } from "../math/kepler";
import type { Epicycle } from "../math/expansion";
import { CERES, EROS, MARS, PHAETHON } from "../math/bodies";

export type BodyId = "mars" | "ceres" | "eros" | "phaethon" | "naff";
export type Mode = "kepler" | "naff";
export type PresetId = "ptolemy" | "fourier" | "ceres-naff";

export const BODY_IDS: readonly BodyId[] = ["mars", "ceres", "eros", "phaethon", "naff"];
export const BODY_LABEL: Record<BodyId, string> = {
  mars: "Mars", ceres: "Ceres", eros: "Eros", phaethon: "Phaethon", naff: "Perturbed Ceres (NAFF)",
};
export const BODY_ELEMENTS: Record<Exclude<BodyId, "naff">, Elements> = {
  mars: MARS, ceres: CERES, eros: EROS, phaethon: PHAETHON,
};
/** Heliocentric orbit drawn for NAFF mode: the unperturbed osculating Ceres ellipse (context only). */
export const NAFF_HELIO_ELEMENTS: Elements = CERES;

/** Max N (analytic) or K (naff) per body. Phaethon (e=0.89) needs N up to ~300. */
export function nMaxFor(body: BodyId): number {
  switch (body) {
    case "phaethon": return 300;
    case "naff": return 200;
    default: return 60;
  }
}
/** Number of drawn circles: analytic chain has 2(2N+1) terms (body + Earth; 4N rotating circles plus two static k=0 offsets); NAFF top-K has 2K+1. */
export function circleCountFor(body: BodyId, n: number): number {
  return body === "naff" ? 2 * n + 1 : 2 * (2 * n + 1);
}
export const BIG_N_WARN = 60;

export const NAFF_EPOCH_JD = 2461200.5;
export const NAFF_SPAN_YEARS = 200;
/** Study window [start, end] (JD): the 200-yr integration starts at the epoch. */
export const NAFF_WINDOW: readonly [number, number] = [NAFF_EPOCH_JD, NAFF_EPOCH_JD + NAFF_SPAN_YEARS * 365.25];
/** Interior (middle 80 %) sub-window where the NAFF window taper does not bias the fit. */
export const NAFF_INTERIOR: readonly [number, number] = [
  NAFF_WINDOW[0] + 0.1 * (NAFF_WINDOW[1] - NAFF_WINDOW[0]), NAFF_WINDOW[0] + 0.9 * (NAFF_WINDOW[1] - NAFF_WINDOW[0]),
];
export const DEFAULT_T = 2452880; // Mars opposition, 2003-08-28
export const NAFF_DEFAULT_T = 0.5 * (NAFF_WINDOW[0] + NAFF_WINDOW[1]);

/**
 * Top-K selection for the NAFF chain: the 2K+1 circles of largest radius
 * (K positive + K negative frequencies + the mean), i.e. what an analytic order-K chain would hold.
 * Returned sorted by radius descending; K is clamped to [1, floor((len-1)/2)].
 */
export function selectTopK(chain: readonly Epicycle[], K: number): Epicycle[] {
  const kMax = Math.max(1, Math.floor((chain.length - 1) / 2));
  const k = Math.min(kMax, Math.max(1, Math.round(Number.isFinite(K) ? K : 1)));
  const sorted = [...chain].sort((a, b) => b.radius - a.radius);
  return sorted.slice(0, Math.min(sorted.length, 2 * k + 1));
}

export interface InitialState {
  body: BodyId;
  mode: Mode;
  N: number;
  t: number;
  speed: number;
  paused: boolean;
  preset: PresetId | null;
  zoom: "opposition" | null;
  showCircles: boolean;
  showError: boolean;
  info: boolean;
  fourier: boolean; // run the 1 -> 20 animation on start
  ptolemaic: boolean; // N=1 only: keep just the k=+1 deferent and k=+1 epicycle (2 circles, paper model A)
}

const isBody = (s: string | null): s is BodyId => s !== null && (BODY_IDS as readonly string[]).includes(s);

/** Defaults implied by a preset (before explicit URL overrides). */
export function presetDefaults(p: PresetId | null): Partial<InitialState> {
  switch (p) {
    case "ptolemy": return { body: "mars", N: 1, t: DEFAULT_T, showError: true, ptolemaic: true };
    case "fourier": return { body: "mars", N: 1, t: DEFAULT_T, fourier: true };
    case "ceres-naff": return { body: "naff", N: 20, t: NAFF_DEFAULT_T };
    default: return {};
  }
}

export function clampToWindow(t: number, body: BodyId): number {
  if (body !== "naff") return t;
  return Math.min(NAFF_WINDOW[1], Math.max(NAFF_WINDOW[0], t));
}

/** Parse URL query (?preset=&body=&N=&t=&speed=&paused=1&zoom=opp&circles=0&error=1&info=0) into a valid state. Never throws. */
export function parseParams(search: string): InitialState {
  const q = new URLSearchParams(search);
  const num = (k: string): number | undefined => {
    const raw = q.get(k);
    if (raw === null || raw.trim() === "") return undefined;
    const v = Number(raw);
    return Number.isFinite(v) ? v : undefined;
  };
  const pr = q.get("preset");
  const preset: PresetId | null = pr === "ptolemy" || pr === "fourier" || pr === "ceres-naff" ? pr : null;
  const d = presetDefaults(preset);
  const bq = q.get("body");
  const body: BodyId = isBody(bq) ? bq : d.body ?? "mars";
  const mode: Mode = body === "naff" ? "naff" : "kepler";
  const nmax = nMaxFor(body);
  const N = Math.min(nmax, Math.max(1, Math.round(num("N") ?? d.N ?? 20)));
  const t = clampToWindow(num("t") ?? (body === "naff" ? NAFF_DEFAULT_T : d.t ?? DEFAULT_T), body);
  const sp = num("speed") ?? 30;
  return {
    body, mode, N, t,
    speed: Math.min(200, Math.max(1, sp)),
    paused: q.get("paused") === "1",
    preset,
    zoom: q.get("zoom") === "opp" ? "opposition" : null,
    showCircles: q.get("circles") !== "0",
    showError: q.has("error") ? q.get("error") === "1" : d.showError ?? false,
    info: q.get("info") !== "0",
    fourier: q.has("fourier") ? q.get("fourier") === "1" : d.fourier ?? false,
    ptolemaic: q.has("ptolemaic") ? q.get("ptolemaic") === "1" : d.ptolemaic ?? false,
  };
}

/** Paper model A: of the N=1 chain keep only the k=+1 body circle (deferent) and k=+1 Earth circle (epicycle). */
export function ptolemaicChain(chain: readonly Epicycle[]): Epicycle[] {
  return chain.filter((c) => c.k === 1);
}

export const FOURIER_N_MAX = 20;
export const FOURIER_SECONDS = 8;
/** N for the Fourier animation after `elapsed` seconds (1 -> 20, then holds). */
export function fourierN(elapsed: number): number {
  if (!(elapsed > 0)) return 1;
  return 1 + Math.floor((FOURIER_N_MAX - 1) * Math.min(1, elapsed / FOURIER_SECONDS));
}

/** Formatting for the study error tables (pure; consumed by the info panel). */
export interface StudyRow { label: string; circles: string; interior: string; full: string }
const sci = (x: unknown): string => (typeof x === "number" && Number.isFinite(x) ? x.toExponential(2) : "n/a");
export function studyRows(study: unknown): StudyRow[] {
  const rows: StudyRow[] = [];
  if (typeof study !== "object" || study === null) return rows;
  const s = study as Record<string, any>;
  const naff = s.naff_errors;
  if (naff && typeof naff === "object") {
    for (const k of Object.keys(naff).sort((a, b) => Number(a) - Number(b))) {
      const r = naff[k] ?? {};
      rows.push({ label: `NAFF K=${k}`, circles: String(r.n_circles ?? 2 * Number(k) + 1), interior: sci(r.interior_max), full: sci(r.full_max) });
    }
  }
  const other = s.kepler_other;
  if (other && typeof other === "object") {
    for (const k of Object.keys(other)) {
      const r = other[k] ?? {};
      rows.push({ label: `Kepler (${k.replace(/^kepler_/, "")})`, circles: "-", interior: sci(r.interior_max), full: sci(r.full_max) });
    }
  }
  const osc = s.kepler_osculating_chain_errors;
  if (osc && typeof osc === "object") {
    const keys = Object.keys(osc).sort((a, b) => Number(a) - Number(b));
    if (keys.length) {
      const r = osc[keys[keys.length - 1]!] ?? {};
      rows.push({ label: "Kepler (osculating)", circles: "-", interior: sci(r.interior_max), full: sci(r.full_max) });
    }
  }
  return rows;
}
