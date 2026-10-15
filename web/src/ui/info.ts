/** Explanatory side panel (KaTeX formulas, bundled locally; no CDN). */
import katex from "katex";
import "katex/dist/katex.min.css";
import type { StudyRow } from "../sim/presets";

export const FORMULAS = {
  geo: String.raw`\mathbf r_{\mathrm{geo}}(t)=\mathbf r_{\mathrm{body}}(t)-\mathbf r_{\oplus}(t)`,
  fourier: String.raw`z(t)=a\sum_{k=-N}^{N}c_k\,e^{\,ik\,M(t)},\qquad M(t)=n\,(t-\tau)`,
  bessel: String.raw`c_{k}=\frac1k\Big[\tfrac{1+\sqrt{1-e^{2}}}{e}J_k(ke)-J_{k+1}(ke)\Big],\quad c_{-k}=\frac1k\Big[\tfrac{e}{1+\sqrt{1-e^{2}}}J_k(ke)-J_{k+1}(ke)\Big],\quad c_0=-\tfrac32 e`,
  naff: String.raw`\mathbf r(t)\approx\sum_{j}A_j\,e^{\,i\omega_j t}`,
} as const;

export function renderFormula(tex: string, display = true): string {
  return katex.renderToString(tex, { displayMode: display, throwOnError: false, output: "html" });
}

function section(parent: HTMLElement, title: string, text: string, tex?: string): void {
  const s = document.createElement("section");
  const h = document.createElement("h3"); h.textContent = title;
  const p = document.createElement("p"); p.textContent = text;
  s.append(h, p);
  if (tex) { const d = document.createElement("div"); d.className = "formula"; d.innerHTML = renderFormula(tex); s.append(d); }
  parent.append(s);
}

export interface InfoPanel {
  setMode(note: string, rows: StudyRow[], windowText: string): void;
}

export function buildInfo(root: HTMLElement): InfoPanel {
  root.replaceChildren();
  const title = document.createElement("h2"); title.textContent = "How to read this";
  const mode = document.createElement("section"); mode.id = "modenote"; mode.setAttribute("aria-live", "polite");
  root.append(title, mode);
  section(root, "Geocentric difference",
    "Seen from Earth, a body sits at the vector difference of two heliocentric Kepler orbits. The right view draws exactly this path (green) and the epicycle chain that reproduces it.",
    FORMULAS.geo);
  section(root, "Fourier series = epicycles",
    "A Kepler orbit is periodic in the mean anomaly M, so z = x + iy has a Fourier series. Each term is a circle of radius a|c_k| turning at k times the mean motion, riding on the previous one. Truncating at |k| <= N gives 2N rotating circles per orbit plus a static offset (the k=0 term); for body and Earth together N=1 is 4 rotating circles plus an offset. The circular two-term approximation (paper model A, the Ptolemaic checkbox) keeps only the k=+1 deferent and k=+1 epicycle.",
    FORMULAS.fourier);
  section(root, "Bessel coefficients",
    "The coefficients are closed-form in Bessel functions of the first kind, so no numerical fitting is needed. Their decay is set by the eccentricity e: Mars (e = 0.09) needs a few circles, Phaethon (e = 0.89) needs hundreds. The bound in the paper certifies an order growing like (1-e)^{-3/2}; the actual minimal order is smaller and saturates as e approaches 1.",
    FORMULAS.bessel);
  section(root, "Why perturbed orbits need NAFF",
    "With planetary perturbations the motion is no longer periodic but almost periodic: a sum of incommensurate frequencies. Fitting them from a finite time window (NAFF, Hann taper) gives a finite circle list that is accurate inside the window but not beyond it.",
    FORMULAS.naff);
  const p = document.createElement("p"); p.className = "paperlink";
  p.append("Full derivation: ");
  const a = document.createElement("a");
  a.href = "https://ebowman.github.io/epicycles/paper.pdf"; a.target = "_blank"; a.rel = "noopener";
  a.textContent = "Geocentric epicycle expansions of Keplerian orbits (PDF)";
  const r = document.createElement("a");
  r.href = "https://github.com/ebowman/epicycles"; r.target = "_blank"; r.rel = "noopener"; r.textContent = "source repository";
  p.append(a, " | ", r);
  root.append(p);
  return {
    setMode(note, rows, windowText) {
      mode.replaceChildren();
      const h = document.createElement("h3"); h.textContent = "Current mode";
      const n = document.createElement("p"); n.textContent = note;
      mode.append(h, n);
      if (rows.length) {
        const w = document.createElement("p"); w.textContent = windowText;
        const det = document.createElement("details"); det.open = true;
        const sum = document.createElement("summary"); sum.textContent = "Study error table (max error, AU)";
        const tb = document.createElement("table");
        tb.createCaption().textContent = "Reconstruction error vs the N-body truth";
        const head = tb.createTHead().insertRow();
        for (const c of ["Model", "Circles", "Interior", "Full"]) {
          const th = document.createElement("th"); th.scope = "col"; th.textContent = c; head.append(th);
        }
        const body = tb.createTBody();
        for (const r of rows) {
          const tr = body.insertRow();
          for (const v of [r.label, r.circles, r.interior, r.full]) tr.insertCell().textContent = v;
        }
        det.append(sum, tb);
        mode.append(w, det);
      }
    },
  };
}
