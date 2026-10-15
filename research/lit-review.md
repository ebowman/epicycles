# Literature review: Newtonian orbits as epicycle (Fourier) systems

Bead: `epicycles-3im.1`. Compiled 2026-10-14.

**Hypothesis under review.** Ptolemaic epicycles are, mathematically, Fourier (more precisely, quasi-periodic / Bohr almost-periodic) series. Any bounded classical Newtonian trajectory seen from Earth can therefore be written to arbitrary precision as a system of epicycles, through an explicit transform from the heliocentric model:

- `r_geo = r_body - r_earth`;
- each Kepler orbit has an exact Fourier series in its mean anomaly with Bessel-function coefficients;
- perturbed motion is quasi-periodic and can be recovered numerically with Laskar's NAFF.

**Verification conventions.**
- **Verified** means the bibliographic data was checked on 2026-10-14 against one of: the Crossref API (`api.crossref.org/works/<doi>`), the arXiv API, the publisher or host page, or (for VSOP87) the ADS bibcode as given in the IMCCE distribution README.
- A sample of DOIs was also resolved through `https://doi.org/`, and each returned an HTTP 302 redirect.
- **UNVERIFIED** entries could not be checked against a primary record. They are kept here as leads and are left out of `refs.bib`.
- **Paywalled** marks entries described from the abstract or from secondary sources only, as the bead requires.

Citation keys in square brackets refer to `research/refs.bib`.

---

## A. Epicycles as Fourier / quasi-periodic series: history and philosophy

**A1. N. R. Hanson (1960), "The Mathematical Power of Epicyclical Astronomy," *Isis* 51(2):150–158. [Hanson1960]**
- DOI: https://doi.org/10.1086/348869 (Crossref verified).
- *Paywalled*: the freely circulating PDF is an image-only scan, so this entry is based on the abstract and secondary summaries.
- **What it shows.** Hanson argues that the deferent-and-epicycle device is computationally universal. Enough circles can reproduce essentially any observed path; secondary accounts, including Wikipedia, credit him with constructions of polygonal paths such as squares and triangles. He also argues that the decline of Ptolemaic astronomy was not caused by weak mathematics.
- **Relation to our hypothesis.** This is the classic statement of the "epicycles can draw anything" thesis, and the first item any reviewer will cite. It is kinematic and representational only. It does not derive epicycles from Newtonian dynamics.

**A2. G. Gallavotti (2001), "Quasi periodic motions from Hipparchus to Kolmogorov," *Rend. Lincei Mat. Appl.* 12:125–152. [Gallavotti2001]**
- arXiv: https://arxiv.org/abs/chao-dyn/9907004 (verified). The journal reference comes from arXiv; the Sapienza IRIS repository dates the issue 2002.
- **What it shows.** This is the closest scholarly prior statement of the thesis. Gallavotti:
  - rewrites a quasi-periodic Fourier series with frequencies ω₁, ω₂ as a deferent with nested epicycles ordered by amplitude (his eqs. 13–15);
  - expands the Kepler ellipse for Mars in powers of *e*, giving "a deferent and two epicycles" at O(e) and more epicycles at higher order (his eqs. 17–20);
  - argues that perturbed motions near integrable ones are quasi-periodic (Kolmogorov) and so are, again, epicycles. Over long times they are only approximated (Poincaré non-convergence, Arnold diffusion).
- He credits Schiaparelli (1874) as the earliest Fourier reading of epicycles (see A9).
- He does **not** write the geocentric difference `r_body - r_earth` explicitly, does not use Bessel coefficients, and does not treat 3D.
- **Relation to our hypothesis.** The conceptual core ("epicycles = quasi-periodic Fourier series, and KAM explains why they work") is already published here. Our paper must cite it prominently and position itself as an explicit, quantitative and geocentric realisation of it.

**A3. Yu. I. Manin (2018), "Time and periodicity from Ptolemy to Schroedinger: paradigm shifts vs continuity," arXiv:1812.00364. [Manin2018]**
- https://arxiv.org/abs/1812.00364 (verified).
- **What it shows.** Manin offers a history-of-mathematics essay that traces one continuous thread, the mathematics of periodicity, from Ptolemaic epicycles to Fourier analysis and quantum amplitudes.
- **Relation to our hypothesis.** It confirms that the "epicycles → Fourier" continuity is a standard theme among mathematicians. It contains no orbital computation.

**A4. D. Boccaletti (2001), "From the epicycles of the Greeks to Kepler's ellipse – the breakdown of the circle paradigm," arXiv:physics/0107009. [Boccaletti2001]**
- https://arxiv.org/abs/physics/0107009 (verified).
- **What it shows.** A historical account of the composition-of-circles paradigm from Hipparchus to Kepler, and of why the *Astronomia Nova* broke it.
- **Relation to our hypothesis.** Useful historical context. It frames Kepler as a break with circles, whereas our framing is that the ellipse is still an infinite epicycle sum.

**A5. E. Scholz (2022), "From heliocentrism to epicycles: a commentary on pre-Ptolemaic astronomy," arXiv:2208.02137. [Scholz2022]**
- https://arxiv.org/abs/2208.02137 (verified).
- **What it shows.** Scholz shows that translating heliocentric uniform circular orbits to an Earth-based observer, using classical geometric means, leads naturally to epicyclic constructions. He suggests that Aristarchus's heliocentrism and the invention of epicycles may be historically linked.
- **Relation to our hypothesis.** This is exactly our transform `r_geo = r_body - r_earth`, restricted to circular orbits (one deferent plus one epicycle). Our contribution extends it to eccentric, inclined and perturbed orbits.

**A6. C. C. Carman (2010), "La refutabilidad del sistema de epiciclos y deferentes de Ptolomeo," *Principia* 14(2). [Carman2010]**
- DOI: https://doi.org/10.5007/1808-1711.2010v14n2p211 (Crossref verified).
- **What it shows.** A philosophical analysis of whether the Ptolemaic system is falsifiable, given that epicycles can fit any curve. It works in the tradition of the "Homer Simpson orbit" demonstration (A10).
- **Relation to our hypothesis.** Supplies the main conceptual objection a referee will raise: "if epicycles can fit anything, the representation is empty." Our answer has to be that, for Newtonian orbits, the epicycle **frequencies and amplitudes are fixed by physics** (n_body, n_earth, e, i, Bessel coefficients). They are not free fit parameters.

**A7. C. C. Carman (2015), "El sistema de epiciclos y deferentes de Ptolomeo: una reconstrucción," *Metatheoria* 6(1):47–72. [Carman2015]**
- DOI: https://doi.org/10.48160/18532330me6.246 (Crossref verified).
- **What it shows.** A structuralist (metatheoretical) reconstruction of the Ptolemaic theory.
- **Relation to our hypothesis.** Philosophical background only.

**A8. G. Acosta, E. Smith, O. Kosheleva, V. Kreinovich (2019), "Epicycles are almost as good as trigonometric series: general system-based analysis," *Appl. Math. Sci.* 13(16):769–773. [Acosta2019]**
- DOI: https://doi.org/10.12988/ams.2019.97103 (Crossref verified). Also issued as a UTEP CS technical report.
- **What it shows.** A short note arguing, on general approximation grounds, that epicycle systems are nearly as expressive as trigonometric series.
- **Relation to our hypothesis.** Weak venue, but it shows the "epicycles ≈ Fourier" claim keeps being re-derived independently.

**A9. G. V. Schiaparelli (1874), memoir cited by Gallavotti as the earliest Fourier reading of epicycles. UNVERIFIED.**
- Gallavotti (A2, footnote 10) cites "a memory of 1874 of G. Schiaparelli."
- The 1874 lecture we could identify is "Le sfere omocentriche di Eudosso, di Callippo e di Aristotele," read 26 Nov 1874 and printed 1877. We could not confirm that it, or "I precursori di Copernico" (1873), contains a Fourier interpretation.
- Lead only; omitted from the .bib.

**A10. S. Ginnobili & C. C. Carman (2008), "Deferentes, epiciclos y adaptaciones" (origin of the 1000-epicycle "Homer Simpson orbit"). UNVERIFIED.**
- The work is referenced in a GeoGebra-journal article and in Carman's later papers, but we could not retrieve the venue or pagination.
- **Relation to our hypothesis.** This is the origin of the viral "epicycles draw Homer Simpson" animation, computed with a DFT, which Mathologer (F1) popularised.
- Omitted from the .bib.

---

## B. Ptolemaic models versus Keplerian motion (quantitative)

**B1. O. Neugebauer (1975), *A History of Ancient Mathematical Astronomy*, Springer. [Neugebauer1975]**
- DOI: https://doi.org/10.1007/978-3-642-61910-6 (Crossref verified).
- **What it shows.** The standard reference for the ancient models. It covers Apollonius's equivalence of the eccentric and epicyclic models, and the Ptolemaic planetary models in full detail.
- **Relation to our hypothesis.** The authority for what Ptolemy actually did, as opposed to the Fourier caricature: the equant, eccentric deferents, and frequencies tied to the solar year.

**B2. N. M. Swerdlow & O. Neugebauer (1984), *Mathematical Astronomy in Copernicus's De Revolutionibus*, Springer. [SwerdlowNeugebauer1984]**
- DOI: https://doi.org/10.1007/978-1-4613-8262-1 (Crossref verified).
- **What it shows.** Copernicus's models still use epicycles (epicyclets replace the equant), and the heliocentric change of frame alone does not remove circles.
- **Relation to our hypothesis.** Historical anchor showing that frame change and circle-composition are independent ideas. This is the point our transform makes explicit.

**B3. J. Evans (1998), *The History and Practice of Ancient Astronomy*, Oxford UP. [Evans1998]**
- DOI: https://doi.org/10.1093/oso/9780195095395.001.0001 (Crossref verified).
- **What it shows.** A textbook treatment of the Ptolemaic models, with hands-on constructions. It includes the equivalence of the epicycle and the Earth's heliocentric orbit for the superior planets.
- **Relation to our hypothesis.** The pedagogical baseline. Our 3D visualisation should be framed against it.

**B4. Ptolemy, *Almagest*, trans. G. J. Toomer (1984; Princeton reissue 1998). [Toomer1998]**
- DOI: https://doi.org/10.1515/9780691213361 (Crossref verified; Crossref dates this DOI record 1999).
- **What it shows.** The primary source: Ptolemy's actual models and parameters.
- **Relation to our hypothesis.** Needed for any "Ptolemy versus explicit Fourier epicycles" comparison.

**B5. J. Evans (1984), "On the function and the probable origin of Ptolemy's equant," *Am. J. Phys.* 52(12):1080–1089. [Evans1984]**
- DOI: https://doi.org/10.1119/1.13764 (Crossref verified).
- **What it shows.** The equant reproduces Kepler's second law to first order in eccentricity.
- **Relation to our hypothesis.** This is an AJP precedent for "Ptolemy approximates Kepler." The equant is a non-uniform device; in Fourier language it is equivalent to the first harmonic correction. Our paper can state this equivalence explicitly.

**B6. J. Evans (1988), "The division of the Martian eccentricity from Hipparchos to Kepler: a history of the approximations to Kepler motion," *Am. J. Phys.* 56(11):1009–1024. [Evans1988]**
- DOI: https://doi.org/10.1119/1.15383 (Crossref verified).
- **What it shows.** Traces successive approximations to Keplerian motion through O(e) and O(e²), comparing the eccentric, the equant and the bisected eccentricity.
- **Relation to our hypothesis.** The closest AJP-style treatment of "historic models as truncations of Kepler motion." It stays in the heliocentric or planet-centred frame and does not use Fourier/Bessel language.

**B7. C. C. Carman & G. L. Recio (2019), "Ptolemaic planetary models and Kepler's laws," *Arch. Hist. Exact Sci.* 73(1):39–124. [CarmanRecio2019]**
- DOI: https://doi.org/10.1007/s00407-018-0219-x (Crossref verified).
- *Paywalled*: described from the abstract.
- **What it shows.** Compares Ptolemaic longitude models with Keplerian models (first two laws) under optimised parameters. It separates structural differences from parameter differences in predictive power.
- **Relation to our hypothesis.** The most rigorous modern quantitative Ptolemy-versus-Kepler comparison. A new paper must not claim to be the first such comparison.

**B8. I. Rushkin (2015), "Optimizing the Ptolemaic model of planetary and solar motion," arXiv:1502.01967. [Rushkin2015]**
- https://arxiv.org/abs/1502.01967 (verified).
- **What it shows.** Ptolemy's equal eccentric/equant shifts reproduce the eccentricity only to first order. A 5/3 ratio would improve accuracy, and the best attainable error is quadratic in *e*.
- **Relation to our hypothesis.** A concrete error-order analysis of the classical model. Our Bessel expansion gives the general all-orders version.

---

## C. Fourier–Bessel (Kapteyn) expansion of Keplerian motion

**C1. F. W. Bessel (1824/1826), "Untersuchung des Theils der planetarischen Störungen, welche aus der Bewegung der Sonne entsteht," *Abh. Königl. Akad. Wiss. Berlin* (1824, publ. 1826), pp. 1–52. UNVERIFIED.**
- The title, pages and dates come from the secondary account in Dutka (C2); the primary text was not inspected.
- **What it shows.** It introduces the Bessel functions Jₙ as the Fourier coefficients of the solution of Kepler's equation, e.g. E − M = Σ (2/n) Jₙ(ne) sin nM.
- **Relation to our hypothesis.** This is the origin of "a Kepler orbit is an exact Fourier series in mean anomaly with Bessel coefficients," and it is fully classical.
- Omitted from the .bib; cite via Dutka (1995).

**C2. J. Dutka (1995), "On the early history of Bessel functions," *Arch. Hist. Exact Sci.* 49(2):105–134. [Dutka1995]**
- DOI: https://doi.org/10.1007/BF00376544 (Crossref verified).
- **What it shows.** A history of how Bessel functions arose from Kepler's equation and planetary perturbations. Contributors include Lagrange, Bessel and others.
- **Relation to our hypothesis.** The citable secondary source for C1.

**C3. J. Dutka (1997), "A note on 'Kepler's equation'," *Arch. Hist. Exact Sci.* 51(1):59–65. [Dutka1997]**
- DOI: https://doi.org/10.1007/BF00376451 (Crossref verified).
- **What it shows.** Historical notes on Kepler's equation and its series solutions.
- **Relation to our hypothesis.** Context for the Lagrange power series in *e* and its convergence. The Laplace limit e ≈ 0.6627 applies to that series. It does **not** apply to the Fourier series in M, which converges for all e < 1. This distinction should be made explicit in our paper.

**C4. W. Kapteyn (1893), "Recherches sur les fonctions de Fourier–Bessel," *Ann. Sci. ÉNS* 10:91–122. [Kapteyn1893]**
- DOI: https://doi.org/10.24033/asens.385 (Crossref verified).
- **What it shows.** Introduces and analyses series of the form Σ aₙ Jₙ(nz), now called "Kapteyn series," including their convergence domain.
- **Relation to our hypothesis.** The Fourier coefficients of position and velocity on a Kepler orbit are Kapteyn-type, so this supplies the convergence theory for our epicycle amplitudes.

**C5. G. N. Watson (1944), *A Treatise on the Theory of Bessel Functions*, 2nd ed., Cambridge UP; 1st ed. 1922. [Watson1944]**
- The book has no DOI. Its existence and edition were verified through Crossref records of contemporary reviews, e.g. *Math. Gazette* review of the 2nd edition, https://doi.org/10.1017/S0025557200206384.
- **What it shows.** Chapter 17 treats Kapteyn series and their application to Kepler's problem, including convergence.
- **Relation to our hypothesis.** The standard reference for the convergence and error bounds of the Bessel/Kapteyn epicycle amplitudes.

**C6. R. C. Tautz, I. Lerche, D. Dominici (2011), "Methods for summing general Kapteyn series," *J. Phys. A* 44:385202. [Tautz2011]**
- DOI: https://doi.org/10.1088/1751-8113/44/38/385202 (Crossref verified).
- **What it shows.** Modern techniques for closed-form summation of Kapteyn series. Kepler-problem series appear as examples.
- **Relation to our hypothesis.** Useful for deriving tail bounds and closed forms of the truncation error.

**C7. R. Borghi (2024), "On the Bessel solution of Kepler's equation," *Mathematics* 12(1):154. [Borghi2024]**
- DOI: https://doi.org/10.3390/math12010154 (Crossref verified); arXiv:2312.01437.
- **What it shows.** Re-examines the Lagrange–Bessel series solution of Kepler's equation using Stieltjes-series theory, and shows the field is still active.
- **Relation to our hypothesis.** Recent and citable evidence on convergence acceleration and truncation behaviour of the Bessel series.

**C8. P. Colwell (1993), *Solving Kepler's Equation over Three Centuries*, Willmann-Bell. [Colwell1993]**
- The book has no DOI. Verified through Crossref's record of E. M. Standish's review, *J. Hist. Astron.* 26 (1995):175–177, https://doi.org/10.1177/002182869502600213.
- **What it shows.** A comprehensive survey of solution methods for Kepler's equation, including the Fourier–Bessel series and the Laplace limit.
- **Relation to our hypothesis.** A one-stop citation for the classical facts we rely on.

**C9. C. D. Murray & S. F. Dermott (1999), *Solar System Dynamics*, Cambridge UP. [MurrayDermott1999]**
- DOI: https://doi.org/10.1017/CBO9781139174817 (Crossref verified; Crossref lists online date 2000).
- **What it shows.** Chapter 2 gives the standard expansions of r, the true anomaly and related quantities in mean anomaly, with Bessel-function coefficients. Later chapters cover the disturbing-function expansions of perturbed motion.
- **Relation to our hypothesis.** The modern textbook source for the heliocentric half of our transform.

**C10. D. Brouwer & G. M. Clemence (1961), *Methods of Celestial Mechanics*, Academic Press. [BrouwerClemence1961]**
- Verified through Crossref chapter records under the prefix 10.1016/B978-1-4832-0075-0, e.g. https://doi.org/10.1016/B978-1-4832-0075-0.50003-4. No book-level DOI was found.
- **What it shows.** The classical treatment of elliptic-motion expansions in Bessel functions and of perturbation theory, including Hansen's method.
- **Relation to our hypothesis.** The classical authority for the expansions.

**C11. G. Pucacco (2012), "Normal forms for the epicyclic approximations of the perturbed Kepler problem," *New Astron.* 17(5):475–482. [Pucacco2012]**
- DOI: https://doi.org/10.1016/j.newast.2012.01.001 (Crossref verified); arXiv:1303.2925.
- **What it shows.** Expands the Kepler orbit in a rotating (Hill) frame as Fourier series in mean anomaly. The first order is a circle plus a second circle at twice the frequency, offset by −3ae/2; the paper calls these "epicyclic expansions." Normal forms extend this to perturbed cases. In the unperturbed case the results coincide with the Kepler-equation expansion.
- **Relation to our hypothesis.** A modern, peer-reviewed use of the word "epicycle" for exactly the Fourier terms we use, though in a co-rotating frame and in the plane. This is strong evidence that the Kepler-as-epicycles expansion is known in the dynamics community. Note that the "epicyclic approximation" of galactic dynamics is a related usage.

---

## D. Newtonian orbit → epicycles, derived explicitly or numerically

**D1. Q. M. Sugon Jr., S. Bragais, D. J. McNamara (2008), "Copernicus's epicycles from Newton's gravitational force law via linear perturbation theory in geometric algebra," arXiv:0807.2708. [Sugon2008]**
- https://arxiv.org/abs/0807.2708 (verified). No journal publication was found.
- **What it shows.** Treats the orbit as a perturbed circle, solves the linearised (Hill-type) equation with an exponential Fourier series, and identifies the harmonics with Copernicus's elements: 0th = eccentric, 1st = deferent, 2nd = epicycle.
- **Relation to our hypothesis.** A direct prior derivation of "Newton ⇒ epicycles." It is first-order, heliocentric and planar. Our transform is all-orders, geocentric and 3D.

**D2. J. A. Peralta, A. Calles, E. Yépez (2003), "El análisis de Fourier de las trayectorias planetarias y el modelo copernicano del sistema solar," *Rev. Mex. Fís.* 49(3):283–289. [Peralta2003]**
- https://rmf.smf.mx/ojs/index.php/rmf/article/view/3165. Verified via the journal's OJS record as returned by search; no DOI. The journal's citation metadata lists the last page incorrectly as "283–0"; the PDF runs to p. 289.
- **What it shows.** Integrates Newton's equations (Verlet) in Cartesian coordinates and applies an FFT to the time series. The Fourier terms "naturally generate the deferent and epicycles of the Copernican model for each planet."
- **Relation to our hypothesis.** **The most direct prior art for the numerical half of our programme:** Newtonian integration → FFT → epicycles. It is heliocentric/Copernican, planar, unperturbed, and uses FFT rather than NAFF. It has no analytic Bessel transform and no geocentric 3D treatment. It must be cited.

**D3. Universidad Distrital (Bogotá), "Aproximación didáctica a los modelos que describen el movimiento de los planetas." UNVERIFIED.**
- Repository record: https://repository.udistrital.edu.co/items/a3c6c368-d4bd-4718-97c2-fc491499d938 (returned HTTP 403; author and year not confirmed).
- **What it shows.** Per a search-result abstract, it applies Verlet and DFT to planetary trajectories and to Halley's comet, following Peralta et al.
- **Relation to our hypothesis.** A possible prior numerical treatment of a small-body (comet) orbit as epicycles.
- Omitted from the .bib.

**D4. F. Paillusson & M. Booth (2020), "Physics and science: the art of taking a stance about undecidable questions," arXiv:2005.08754. [Paillusson2020]**
- https://arxiv.org/abs/2005.08754 (verified).
- **What it shows.** Uses computed heliocentric trajectories and their geocentric coordinate transforms to argue kinematic equivalence (underdetermination) between the two world systems.
- **Relation to our hypothesis.** Shows the geocentric transform is used as an illustration in the philosophy of science. No Fourier or epicycle decomposition is given.

**D5. "epa058" (n.d.), *Ptolemy-Fourier-Solar-System*, GitHub repository. [epa058]**
- https://github.com/epa058/Ptolemy-Fourier-Solar-System (verified to exist, accessed 2026-10-14; not peer reviewed).
- **What it shows.** Computes geocentric planetary trajectories from Keplerian elements via Kepler's equation and reconstructs them as sums of rotating circles (FFT peaks plus least-squares). It is 2D only; 3D inclination and real ephemerides are listed as future work.
- **Relation to our hypothesis.** Hobbyist prior art for the geocentric numerical decomposition. It shows that the "geocentric Fourier-epicycle reconstruction" idea is in circulation, so novelty cannot rest on that idea alone.

---

## E. Almost-periodic functions and frequency analysis (NAFF)

**E1. H. Bohr (1925), "Zur Theorie der fast periodischen Funktionen. I," *Acta Math.* 45:29–127. [Bohr1925]**
- DOI: https://doi.org/10.1007/BF02395468 (Crossref verified).
- **What it shows.** Founds the theory of almost-periodic functions. These are the uniform closures of trigonometric polynomials with arbitrary real frequencies, which generalises Fourier series to non-commensurate frequencies.
- **Relation to our hypothesis.** The correct mathematical home for the claim. A geocentric orbit `r_body − r_earth` combines two (generally incommensurate) mean motions, so it is not periodic and has no Fourier series in the strict sense. It is quasi-periodic and therefore Bohr almost-periodic. "Epicycles = Fourier" should be stated as "epicycles = Bohr almost-periodic expansion."

**E2. J. Laskar (1990), "The chaotic motion of the solar system: a numerical estimate of the size of the chaotic zones," *Icarus* 88(2):266–291. [Laskar1990]**
- DOI: https://doi.org/10.1016/0019-1035(90)90084-M (Crossref verified).
- **What it shows.** Introduces numerical frequency analysis (NAFF) on secular solar-system solutions and shows the inner solar system is chaotic.
- **Relation to our hypothesis.** The origin of NAFF. It also sets the **limit** of the hypothesis: chaotic motion is quasi-periodic only over finite times.

**E3. J. Laskar, C. Froeschlé, A. Celletti (1992), "The measure of chaos by the numerical analysis of the fundamental frequencies. Application to the standard mapping," *Physica D* 56(2–3):253–269. [Laskar1992]**
- DOI: https://doi.org/10.1016/0167-2789(92)90028-L (Crossref verified).
- **What it shows.** Formalises NAFF and its accuracy, which is much better than plain FFT for quasi-periodic signals.
- **Relation to our hypothesis.** The methods reference for our perturbed-case demonstration.

**E4. J. Laskar (1993), "Frequency analysis for multi-dimensional systems. Global dynamics and diffusion," *Physica D* 67(1–3):257–281. [Laskar1993]**
- DOI: https://doi.org/10.1016/0167-2789(93)90210-R (Crossref verified).
- **What it shows.** Frequency map analysis for multi-degree-of-freedom systems, and its use to detect diffusion.
- **Relation to our hypothesis.** Supports using frequency drift as a diagnostic of where an epicycle representation breaks down.

**E5. J. Laskar (1999), "Introduction to frequency map analysis," in *Hamiltonian Systems with Three or More Degrees of Freedom*, Springer, pp. 134–150. [Laskar1999]**
- DOI: https://doi.org/10.1007/978-94-011-4673-9_13 (Crossref verified).
- **What it shows.** A tutorial exposition of NAFF and frequency maps.
- **Relation to our hypothesis.** The pedagogical NAFF reference.

**E6. J. Laskar (2003), "Frequency map analysis and quasiperiodic decompositions," arXiv:math/0305364. [Laskar2003]**
- https://arxiv.org/abs/math/0305364 (verified).
- **What it shows.** Convergence theory for NAFF quasi-periodic approximations, with error estimates in terms of window and data length.
- **Relation to our hypothesis.** Provides the error bounds for the numerical (perturbed) half of the paper.

**E7. Y. N. Fu & J. Laskar (2019), "Frequency analysis and representation of slowly diffusing planetary solutions," *A&A* 628:A84. [FuLaskar2019]**
- DOI: https://doi.org/10.1051/0004-6361/201730492 (Crossref verified); arXiv:1701.07108.
- **What it shows.** Notes that planetary ephemerides are traditionally represented as finite sums of periodic (or Poisson) terms. It then extends the representation to solutions with slowly drifting frequencies.
- **Relation to our hypothesis.** The state of the art for "how far can a sum of epicycles go" on long timescales. Over Myr spans the frequencies themselves drift, and a fixed epicycle set fails.

**E8. P. Bretagnon & G. Francou (1988), "Planetary theories in rectangular and spherical variables. VSOP87 solutions," *A&A* 202:309. [Bretagnon1988]**
- ADS bibcode 1988A&A...202..309B. Bibliographic data verified from the IMCCE distribution README, https://ftp.imcce.fr/pub/ephem/planets/vsop87/README. The journal has no DOI for this volume.
- **What it shows.** The VSOP87 analytical planetary theory: heliocentric coordinates as long sums of Poisson terms t^α · cos(A + Bt), fitted to DE200.
- **Relation to our hypothesis.** The operational embodiment of "planetary motion = epicycles." With α = 0, each term is literally a uniform circular motion, and VSOP87 has thousands of them. These are heliocentric, so the geocentric version is VSOP87(planet) − VSOP87(Earth). This is the strongest evidence that the representation itself is not novel in practice.

**E9. S. M. Kudryavtsev (2016), "Analytical series representing DE431 ephemerides of terrestrial planets," *MNRAS* 456(4):4015–4019. [Kudryavtsev2016]**
- DOI: https://doi.org/10.1093/mnras/stv2892 (Crossref verified).
- **What it shows.** Spectral (frequency-analysis) fitting of a modern numerical ephemeris (DE431) to Poisson series over long spans.
- **Relation to our hypothesis.** A modern counterpart of E8 obtained by frequency analysis of a numerical ephemeris. It is methodologically the same as our proposed NAFF step, again heliocentric.

---

## F. Popular and educational

**F1. Mathologer (B. Polster) (2018), "Epicycles, complex Fourier series and Homer Simpson's orbit," YouTube. [Mathologer2018]**
- https://www.youtube.com/watch?v=qS4H6PEcCCA (title verified; air date 6 Jul 2018 per TheTVDB listing).
- **What it shows.** The complex-Fourier-series explanation of drawing arbitrary curves with epicycles. It is built on the Ginnobili–Carman animation and discusses the Moon's orbit.
- **Relation to our hypothesis.** Establishes that "epicycles = complex Fourier series" is widely known popular folklore.

**F2. 3Blue1Brown (G. Sanderson) (2019), "But what is a Fourier series? From heat flow to drawing with circles (DE4)," YouTube. [3b1b2019]**
- https://www.youtube.com/watch?v=r6sGWTCMz2k (title verified; the 2019 year is from the DE series chronology, not from the page metadata).
- **What it shows.** A visual derivation of complex Fourier coefficients as rotating vectors that draw arbitrary closed curves.
- **Relation to our hypothesis.** The visualisation idiom our SPA would extend from 2D closed curves to 3D quasi-periodic orbits.

**F3. Wikipedia, "Deferent and epicycle" (accessed 2026-10-14). [WikiDeferent]**
- https://en.wikipedia.org/wiki/Deferent_and_epicycle.
- **What it shows.** States that any smooth curve can be approximated to arbitrary accuracy with enough epicycles, citing Hanson. It also discusses the Fourier connection.
- **Relation to our hypothesis.** Documents that the general claim is encyclopedia-level common knowledge.

---

## Novelty assessment

### What is classical (established, citable, and not ours to claim)
1. **Kepler motion has an exact Fourier series in mean anomaly with Bessel coefficients.**
   - Lagrange and Bessel (1824); Kapteyn series; Watson ch. 17; Brouwer & Clemence; Murray & Dermott ch. 2.
   - The series converges for all e < 1. The Laplace limit concerns the power series in *e*, not this series.
2. **A heliocentric circular orbit seen from Earth is a deferent plus an epicycle.**
   - Apollonius's equivalence; Copernicus; Neugebauer; Evans; Scholz 2022.
3. **Ptolemy's devices are truncations of Kepler motion at O(e).**
   - Evans 1984 and 1988 (AJP); Rushkin 2015; Carman & Recio 2019.
4. **Quasi-periodic motions are epicycles, and KAM/Kolmogorov explains why near-integrable dynamics looks epicyclic.**
   - Stated and argued explicitly by Gallavotti (2001), including the Kepler ellipse as a deferent plus epicycle series.
5. **Planetary positions in practice are represented as huge sums of uniform circular motions,** obtained analytically or by frequency analysis of numerical ephemerides.
   - VSOP87; Laskar's NAFF; Kudryavtsev 2016; Fu & Laskar 2019.
6. **Newtonian orbit → (FFT or Fourier) → deferent and epicycles, demonstrated directly.**
   - Peralta, Calles & Yépez 2003 (numerical, heliocentric); Sugon et al. 2008 (analytic, first order); Pucacco 2012 ("epicyclic expansions" in Hill's frame).

### What is folklore (true in spirit, usually stated loosely or slightly wrongly)
- **"Epicycles are Fourier series."** This appears in Hanson, Mathologer, 3Blue1Brown, Wikipedia and Acosta et al. Three caveats are routinely glossed over:
  - **(a)** A geocentric orbit involves at least two incommensurate frequencies (n_body and n_earth). It is therefore quasi-periodic, and the right object is a Bohr almost-periodic expansion, not a Fourier series.
  - **(b)** Ptolemy's real models use the equant and eccentrics and are not harmonic-lattice expansions.
  - **(c)** "Epicycles can fit anything" makes the representation empirically vacuous (Carman 2010) **unless** the frequencies and amplitudes are fixed by dynamics.
- **"Any trajectory to arbitrary precision."** On a finite time window this is trivial: Weierstrass approximation by trigonometric polynomials applies to any continuous curve. Beyond that window it holds only for bounded quasi-periodic motion. It fails globally for unbound (hyperbolic or parabolic) trajectories and holds only for finite times on chaotic ones (Laskar 1990; Fu & Laskar 2019).

### What we did *not* find, and so what a new paper can credibly contribute
We found no source that does the following together:

- **(i)** Writes down the **explicit closed-form 3D geocentric epicycle expansion** of an arbitrary Keplerian body. This is r_geo(t) = Σₖ Aₖ(e_b, i_b, ω_b, Ω_b) e^{ik n_b t} − Σⱼ Bⱼ(e_E) e^{ij n_E t}, with Bessel/Kapteyn amplitudes. For each body it has the "ellipse = two counter-rotating circles in eccentric anomaly, re-expanded in mean anomaly" structure. Inclination is handled as a 3D sum of uniformly rotating vectors (circles in tilted planes), not as a planar approximation.
- **(ii)** Gives **rigorous truncation-error bounds** (Bessel/Kapteyn tail estimates) and the number of epicycles needed for a given angular accuracy as a function of *e*, for both the body and the Earth term.
- **(iii)** Demonstrates on a **real perturbed small body** (for example an NEA such as Apophis, or Ceres) that **NAFF applied to a geocentric N-body ephemeris** recovers the analytic epicycles, plus perturbation lines at frequency combinations, with quantified residuals and a clear statement of the time window and the frequency-drift limits.
- **(iv)** **Delimits precisely** where the "arbitrary precision" claim holds:
  - bounded quasi-periodic motion: yes;
  - chaotic motion: finite-time only;
  - unbound orbits: no.

The closest prior items each cover only a piece:
- Gallavotti covers the concept and the KAM framing, but is planar, has no explicit geocentric transform and no Bessel coefficients.
- Peralta et al. 2003 and the epa058 repository cover the numerical FFT pipeline, but are planar, unperturbed, without error bounds, and either heliocentric or not peer reviewed.
- Pucacco and Sugon cover the analytic epicycle expansion, but in a rotating or heliocentric frame and at first order.
- VSOP87 and Kudryavtsev have the representation, but heliocentric and not framed as epicycles.

### Candid bottom line
- The **mathematical isomorphism is not new.** A paper claiming to *discover* that epicycles equal Fourier analysis, or that Newtonian orbits are epicycle systems, would be rejected as known (Gallavotti 2001, Bessel 1824, VSOP87).
- The **specific, explicit, quantitative geocentric-3D synthesis is not, to our search, in the literature.** That covers the closed-form transform, error bounds, the NAFF validation on a real perturbed body, and the precise validity domain. It would make a solid **expository and pedagogical research article**.

### Recommendation: how to frame the paper
- **Venue.** Frame it as a *synthesis-with-new-results* paper for *American Journal of Physics* or *European Journal of Physics* (alternatively *Celestial Mechanics and Dynamical Astronomy* education/notes, or *Archive for History of Exact Sciences* if the historical comparison is expanded).
- **Working title.** Something like *"Ptolemy's epicycles as the geocentric Fourier–Bessel expansion of Newtonian orbits: an explicit transform, error bounds, and a frequency-analysis demonstration."*
- **Opening.** State up front that the epicycle/Fourier correspondence is classical (Hanson, Gallavotti, Bessel/Kapteyn, VSOP87).
- **Claimed contributions.** Claim exactly items (i)–(iv) above, plus the open-source 3D interactive visualisation as a pedagogical contribution.
- **Answer to the "epicycles fit anything" objection (Carman).** In the Newtonian construction, every epicycle radius and rate is *determined* by (a, e, i, ω, Ω, M₀) of the body and of the Earth. Epicycles become a *derived coordinate representation*, not a fitting model.
- **Phrasing.** Use "almost-periodic (Bohr)" rather than "Fourier" wherever two or more incommensurate frequencies appear.
- **Overreach to avoid.** Avoid "any classical trajectory": restrict the claim to bounded quasi-periodic motion and finite-time approximation of chaotic motion.

---

## Source count

| | Count |
|---|---|
| Total sources | 46 |
| Verified (in `refs.bib`) | 42 |
| UNVERIFIED (A9 Schiaparelli 1874, A10 Ginnobili & Carman 2008, C1 Bessel 1824 primary, D3 Universidad Distrital thesis; omitted from the .bib) | 4 |
