# Derivation: Newtonian orbits seen from Earth as epicycle systems

Bead `epicycles-3im.2`. Input: `research/lit-review.md` (notation and framing), `research/refs.bib` (keys in square brackets). All numerical statements below were produced by scripts run on 2026-10-14 (numpy/scipy); the "Claims to verify" in Part 6 are phrased so a pytest can reproduce them.

**Framing (from the lit review, binding).**
- The heliocentric Fourier–Bessel (Kapteyn) expansion of Kepler motion is classical [Kapteyn1893; Watson1944; BrouwerClemence1961; MurrayDermott1999; Dutka1995]. Part 1 is background, re-derived for self-containedness and to fix notation, not claimed as new.
- The epicycle/quasi-periodic correspondence is classical [Hanson1960; Gallavotti2001].
- A geocentric orbit carries at least two incommensurate frequencies, so the right object is a Bohr almost-periodic expansion [Bohr1925], not a Fourier series in one period.
- The new material is the explicit closed-form 3D geocentric assembly (Parts 2–3), the rigorous truncation bounds with the N-needed table (Part 4), and the precise domain of validity (Parts 5 and 7).

---

## Notation

| Symbol | Meaning |
|---|---|
| $a, e$ | semi-major axis, eccentricity of a Keplerian orbit; $s \equiv \sqrt{1-e^2}$ |
| $M, E$ | mean and eccentric anomaly; Kepler's equation $M = E - e\sin E$ |
| $n$ | mean motion; $M = n(t-t_p)$ |
| $z = x+iy$ | position in the perifocal plane, $x$ toward periapsis, origin at the focus |
| $c_k$, $k\in\mathbb Z$ | complex Fourier coefficient: $z = a\sum_{k\in\mathbb Z} c_k e^{ikM}$ (all $c_k$ real) |
| $J_k$ | Bessel function of the first kind, order $k$ |
| $\rho(e)$ | $e\,e^{s}/(1+s)$, the Kapteyn decay ratio, $<1$ for $0\le e<1$ |
| $\Omega, i, \omega$ | longitude of node, inclination, argument of periapsis |
| $\mathsf R$ | $\mathsf R_z(\Omega)\mathsf R_x(i)\mathsf R_z(\omega)$, perifocal to ecliptic |
| $\mathbf P,\mathbf Q,\mathbf N$ | columns of $\mathsf R$: periapsis direction, in-plane quadrature direction, orbit normal |
| $A_k\in\mathbb C^3$ | vector amplitude of harmonic $k\ge1$: term $\mathrm{Re}(A_k e^{ikM})$ |
| $b$, $\oplus$ | subscripts for the observed body and for Earth |
| $N$ | truncation order: keep $\lvert k\rvert\le N$ |
| $\varepsilon_N(a,e)$ | rigorous sup-norm truncation bound (Part 4) |
| $\boldsymbol\omega$ | frequency vector of a quasi-periodic motion (Part 5) |

---

## Part 1. Planar Kepler orbit as a Fourier series in the mean anomaly

### 1.1 Setup
In the perifocal frame, with the focus at the origin,
$$x = a(\cos E - e),\qquad y = a\sqrt{1-e^2}\,\sin E,\qquad M = E - e\sin E .$$
The map $M\mapsto E$ is a smooth, strictly increasing bijection of the circle for $0\le e<1$, so $z(M)$ is real-analytic and $2\pi$-periodic in $M$. It therefore has a Fourier series, and the series converges uniformly (the decay rate is made explicit in 1.4).

### 1.2 Ellipse as two counter-rotating circles in $E$
Using $\cos E=\tfrac12(e^{iE}+e^{-iE})$ and $\sin E=\tfrac1{2i}(e^{iE}-e^{-iE})$:
$$\frac{z}{a} = -e + \frac{1+s}{2}\,e^{iE} + \frac{1-s}{2}\,e^{-iE}. \tag{1.1}$$
This is already an epicycle system in the *eccentric* anomaly: a prograde circle of radius $(1+s)/2$, a retrograde circle of radius $(1-s)/2$, and an offset $-e$. The remaining work is to re-expand $e^{\pm iE}$ in the uniform variable $M$. This is the classical step [Bessel 1824 as cited in Dutka1995; Watson1944 ch. 17].

### 1.3 Coefficients
**Lemma.** For integers $m\neq0$ and $k\neq0$,
$$\frac1{2\pi}\int_0^{2\pi} e^{imE}\,e^{-ikM}\,dM = \frac mk\,J_{k-m}(ke).$$

*Proof.* Since $e^{-ikM}dM=d(e^{-ikM})/(-ik)$, integrating by parts (the boundary term vanishes because $e^{imE}e^{-ikM}$ is $2\pi$-periodic in $M$) gives
$$\frac1{2\pi}\int e^{imE}e^{-ikM}dM=\frac{1}{2\pi\,ik}\int e^{-ikM}\,(im)\,e^{imE}dE=\frac mk\cdot\frac1{2\pi}\int_0^{2\pi}e^{-i[(k-m)E-ke\sin E]}dE=\frac mk\,J_{k-m}(ke),$$
where $M=E-e\sin E$ was inserted, using Bessel's integral $J_n(x)=\frac1{2\pi}\int_0^{2\pi}e^{i(n\tau-x\sin\tau)}d\tau$ (real, even in the sign of the exponent). $\square$

Apply it to (1.1), with $m=+1$ for the prograde circle and $m=-1$ for the retrograde one.

**Coefficients for $k\ge1$.**
$$c_k=\frac{(1+s)J_{k-1}(ke)-(1-s)J_{k+1}(ke)}{2k}
=\frac1k\Big[J_k'(ke)+\frac{s}{e}J_k(ke)\Big], \tag{1.2}$$
using $J_{k-1}-J_{k+1}=2J_k'$ and $J_{k-1}+J_{k+1}=2J_k/x$ (at $x=ke$ this is $2J_k(ke)/(ke)$, which gives the $1/e$ after the $1/(2k)$ factor... explicitly, $J_{k-1}(ke)+J_{k+1}(ke)=\tfrac{2}{e}J_k(ke)$).

**Coefficients for $-k$, $k\ge1$** (apply the lemma at $-k$, using $J_{-n}(-x)=J_n(x)$):
$$c_{-k}=\frac{(1-s)J_{k-1}(ke)-(1+s)J_{k+1}(ke)}{2k}
=\frac1k\Big[J_k'(ke)-\frac{s}{e}J_k(ke)\Big]. \tag{1.3}$$

**Equivalent form without derivatives** (using $J_k'(ke)=J_k(ke)/e-J_{k+1}(ke)$, which follows from the two recurrences):
$$c_k=\frac1k\Big[\frac{1+s}{e}J_k(ke)-J_{k+1}(ke)\Big],\qquad c_{-k}=\frac1k\Big[\frac{1-s}{e}J_k(ke)-J_{k+1}(ke)\Big]. \tag{1.4}$$

**The $k=0$ term.** Averaging (1.1) over $M$: $\langle\cos E\rangle_M=-e/2$ (because $dM=(1-e\cos E)dE$ gives $\frac1{2\pi}\int\cos E(1-e\cos E)dE=-e/2$), and $\langle\sin E\rangle_M=0$. So
$$c_0 = -\tfrac{3e}{2}\quad(\text{i.e. the mean of } x/a \text{ is } -3e/2 \text{ and the mean of } y \text{ is } 0). \tag{1.5}$$

### 1.4 Real form (the form in the task brief)
Write $z/a=c_0+\sum_{k\ge1}\big[\alpha_k\cos kM+i\beta_k\sin kM\big]$ with $\alpha_k=c_k+c_{-k}$ and $\beta_k=c_k-c_{-k}$. From (1.2)–(1.3):
$$\frac{x}{a}=-\frac{3e}{2}+\sum_{k\ge1}\frac2k J_k'(ke)\cos kM,\qquad
\frac{y}{a}=\frac{2\sqrt{1-e^2}}{e}\sum_{k\ge1}\frac1k J_k(ke)\sin kM. \tag{1.6}$$
**Independent derivation of (1.6)** (trigonometric route, as cross-check). The cosine coefficient of $x/a=\cos E-e$ is $\frac2\pi\int_0^\pi(\cos E-e)\cos kM\,dM$. Integrating by parts ($\sin kM$ vanishes at $M=0,\pi$) gives $\frac{2}{\pi k}\int_0^\pi\sin E\,\sin(kE-ke\sin E)\,dE=\frac2kJ_k'(ke)$, using $J_k'(x)=\frac1\pi\int_0^\pi\sin\tau\sin(k\tau-x\sin\tau)d\tau$. For $y$, integrating by parts ($y=0$ at $E=0,\pi$) gives $\frac{2s}{\pi k}\int_0^\pi\cos E\cos(kE-ke\sin E)dE=\frac{2s}{k}\cdot\frac{J_{k-1}+J_{k+1}}{2}=\frac{2s}{ke}J_k(ke)$. Both agree with (1.6) and with the brief. (Numerical confirmation in Part 4.)

### 1.5 Checks and small-$e$ behaviour
- **$e\to0$.** $J_1'(0)=\tfrac12$, $J_1(ke)/e\to \tfrac12$ for $k=1$, and $J_k(ke)=O(e^k)$. So $c_1\to1$, $c_{-1}\to0$ ($c_{-1}=O(e^2)$), all other $c_k\to0$, $c_0\to0$. Thus $z/a\to e^{iM}$: a single circle of radius $a$, as required.
- **Leading terms** (series in $e$, checked numerically): $c_0=-\tfrac32e$, $c_1=1-\tfrac12e^2+\dots$, $c_2=\tfrac12e+\dots$, $c_3=\tfrac38e^2+\dots$, $c_{-1}=\tfrac18e^2+\dots$, $c_{-2}=\tfrac1{24}e^3+\dots$. So $z/a\approx-\tfrac32e+e^{iM}+\tfrac e2e^{2iM}+O(e^2)$, the familiar first-order-in-$e$ result.
- **Periapsis sum rule.** At $M=0$, $E=0$, $z=a(1-e)$, so $-\tfrac{3e}2+\sum_{k\ge1}\frac2kJ_k'(ke)=1-e$ (checked numerically, $e=0.3$).

### 1.6 One term equals one epicycle
The term $a c_k e^{ikM}$, $M=n(t-t_p)$, is a uniform circular motion with
- radius $a|c_k|$,
- angular rate $kn$ (counter-clockwise for $k>0$, clockwise for $k<0$),
- phase $\arg c_k$ at $t=t_p$ (here $0$ or $\pi$; all $c_k$ are real).

The $k=0$ term is a fixed offset of the centre of the whole system (the "eccentric" displacement $-\tfrac32 ea$ along the apsidal line). Hence the orbit is exactly a deferent-and-epicycle hierarchy with an infinite family of circles centred on one point. Geometrically each circle is concentric; the "epicycle on a deferent" picture is the same sum grouped as nested circles.

### 1.7 Convergence and decay of $|c_k|$
**Kapteyn's inequality** [Kapteyn1893; Watson1944 §8.7]: for $0<z\le1$ and $\nu\ge0$,
$$|J_\nu(\nu z)|\le\Big(\frac{z\,e^{\sqrt{1-z^2}}}{1+\sqrt{1-z^2}}\Big)^{\nu}.$$
With $z=e$ this gives
$$|J_k(ke)|\le\rho(e)^k,\qquad \rho(e)=\frac{e\,e^{\sqrt{1-e^2}}}{1+\sqrt{1-e^2}}<1\quad(0\le e<1). \tag{1.7}$$
$\rho$ is increasing with $\rho(0)=0$ and $\rho(1)=1$. (The asymptotic form $J_k(ke)\sim\rho^k/\sqrt{2\pi k}\,(1-e^2)^{-1/4}$ for fixed $e$, large $k$, is sharper by $k^{-1/2}$ but is not needed for a rigorous bound.) Near $e=1$, with $\delta=1-e$: $\ln\rho\approx-(2\delta)^{3/2}/3$.

**Bound on $c_k$.** Since $\rho$ is increasing, $J_{k+1}(ke)=J_{k+1}\big((k+1)\tfrac{ke}{k+1}\big)$ obeys $|J_{k+1}(ke)|\le\rho(ke/(k+1))^{k+1}\le\rho^{k+1}$. From (1.4):
$$|c_k|\le\frac{\rho^k}{k}\Big[\frac{1+s}{e}+\rho\Big],\qquad|c_{-k}|\le\frac{\rho^k}{k}\Big[\frac{1-s}{e}+\rho\Big]. \tag{1.8}$$
Adding,
$$|c_k|+|c_{-k}|\le\frac{2\rho^k}{k}\Big[\frac1e+\rho\Big]. \tag{1.9}$$
(The bounds (1.7)–(1.9) were verified numerically for $k\le399$ and $e\in\{0.0167,0.0934,0.25,0.6,0.89,0.99\}$.)

**Convergence.** Because $\rho<1$, $\sum_k|c_k|<\infty$ for every $0\le e<1$, so the series converges absolutely and uniformly in $M$ (Weierstrass M-test). The constant $1/e$ in (1.9) is harmless: for small $e$, $\rho^k\sim(1.36e)^k$ cancels it for $k\ge2$, and the $k=1$ term is bounded by the true value $c_1\le1$.

Note: the Laplace limit $0.6627$ concerns the power series of $E$ in powers of $e$ (a singularity at complex $e$ of modulus $0.6627$ in the variable that uses $1+e^2$). It does **not** limit this Fourier series, which converges for all $e<1$ (lit review, Novelty assessment, item 1).

---

## Part 2. Three dimensions: inclined orbit as sums of circles in the orbital plane

### 2.1 Rotation to the ecliptic frame
$\mathbf r=\mathsf R\,(x,y,0)^{\mathsf T}$ with $\mathsf R=\mathsf R_z(\Omega)\mathsf R_x(i)\mathsf R_z(\omega)$. Its columns are
$$\mathbf P=\begin{pmatrix}\cos\Omega\cos\omega-\sin\Omega\cos i\sin\omega\\ \sin\Omega\cos\omega+\cos\Omega\cos i\sin\omega\\ \sin i\sin\omega\end{pmatrix},\quad
\mathbf Q=\begin{pmatrix}-\cos\Omega\sin\omega-\sin\Omega\cos i\cos\omega\\ -\sin\Omega\sin\omega+\cos\Omega\cos i\cos\omega\\ \sin i\cos\omega\end{pmatrix},\quad
\mathbf N=\begin{pmatrix}\sin\Omega\sin i\\ -\cos\Omega\sin i\\ \cos i\end{pmatrix}.$$
They are orthonormal with $\mathbf P\times\mathbf Q=\mathbf N$. Then $\mathbf r=x\mathbf P+y\mathbf Q$.

### 2.2 Vector amplitude of each harmonic
Since $x\mathbf P+y\mathbf Q=\mathrm{Re}\big[z\,(\mathbf P-i\mathbf Q)\big]$ and $\mathrm{Re}[c\,e^{-ikM}(\mathbf P-i\mathbf Q)]=\mathrm{Re}[\bar c\,e^{ikM}(\mathbf P+i\mathbf Q)]$,
$$\mathbf r(M)=a\,c_0\,\mathbf P+\sum_{k\ge1}\mathrm{Re}\big(A_k e^{ikM}\big),\qquad
A_k=a\big[c_k(\mathbf P-i\mathbf Q)+c_{-k}(\mathbf P+i\mathbf Q)\big]\in\mathbb C^3, \tag{2.1}$$
($c_{-k}$ real). Equivalently $A_k=a[\alpha_k\mathbf P-i\beta_k\mathbf Q]$ with $\alpha_k=c_k+c_{-k}=\tfrac2kJ_k'(ke)$, $\beta_k=c_k-c_{-k}=\tfrac{2s}{ke}J_k(ke)$.

### 2.3 Each vector term is an ellipse in the orbital plane
$\mathrm{Re}(A_ke^{ikM})=a\alpha_k\cos(kM)\,\mathbf P+a\beta_k\sin(kM)\,\mathbf Q$. This is an ellipse with centre at the origin, semi-axes $a|\alpha_k|$ along $\mathbf P$ and $a|\beta_k|$ along $\mathbf Q$, in the plane with normal $\mathbf N$, traversed at rate $kn$.

### 2.4 Ellipse as two counter-rotating circles
$$\mathrm{Re}(A_ke^{ikM})=\underbrace{a|c_k|\big[\cos(kM+\phi_k)\mathbf P+\sin(kM+\phi_k)\mathbf Q\big]}_{\text{prograde circle: rate }+kn\text{ about }\mathbf N}
+\underbrace{a|c_{-k}|\big[\cos(kM-\phi_{-k})\mathbf P-\sin(kM-\phi_{-k})\mathbf Q\big]}_{\text{retrograde circle: rate }-kn\text{ about }\mathbf N}, \tag{2.2}$$
with $\phi_{\pm k}=\arg c_{\pm k}\in\{0,\pi\}$. Radii $a|c_k|$ and $a|c_{-k}|$; semi-major axis of the ellipse $=a(|c_k|+|c_{-k}|)$, semi-minor axis $=a\,\big||c_k|-|c_{-k}|\big|$; the ellipse is circular iff $c_{-k}=0$. The orientation of the major axis is at angle $\tfrac12(\phi_k+\phi_{-k})$ from $\mathbf P$.

### 2.5 Remark: planar case
Because the Kepler orbit is planar, **all** the circles in (2.2), for every $k$, lie in the single orbital plane $\perp\mathbf N$. So for one body the decomposition is simply the planar statement of 1.6 with $\mathbf P,\mathbf Q$ as the basis; "in-plane" the terms $c_ke^{ikM}$ with $k>0$ and $k<0$ are already counter-rotating circles about $\mathbf N$, with no further work. The inclination matters only once two orbital planes are combined in the geocentric difference (Part 3): the body's circles then live in a plane with normal $\mathbf N_b$ and Earth's in one with normal $\mathbf N_\oplus$ (Earth's plane is the ecliptic, $\mathbf N_\oplus\approx\hat z$ up to a slowly varying tilt). The geocentric epicycle system is therefore a set of circles in two (generally different) planes, giving ecliptic latitude without extra machinery.

---

## Part 3. Geocentric transform

### 3.1 Two epicycle families, no cross terms
With $\mathbf r_{\rm geo}(t)=\mathbf r_b(t)-\mathbf r_\oplus(t)$ and $M_b=n_b(t-t_b)$, $M_\oplus=n_\oplus(t-t_\oplus)$, absorbing the perihelion times into the amplitudes ($A_k\to A_ke^{-ikn_bt_b}$ etc.),
$$\mathbf r_{\rm geo}(t)=\underbrace{a_bc_0^{(b)}\mathbf P_b-a_\oplus c_0^{(\oplus)}\mathbf P_\oplus}_{\text{constant offset}}+\sum_{k\ge1}\mathrm{Re}\big(A^{(b)}_ke^{ikn_bt}\big)-\sum_{j\ge1}\mathrm{Re}\big(A^{(\oplus)}_je^{ijn_\oplus t}\big). \tag{3.1}$$
The spectrum is exactly the union of the two harmonic ladders $\{kn_b\}\cup\{jn_\oplus\}$, because the difference is linear. There are **no** combination frequencies $kn_b\pm jn_\oplus$ in the Cartesian vector. Combination frequencies appear only in nonlinear functions of $\mathbf r_{\rm geo}$ (ecliptic longitude, distance, angular rate).

### 3.2 Almost-periodic, not periodic
Both ladders are uniformly (absolutely) summable by (1.9). So (3.1) is a uniformly convergent trigonometric series with two base frequencies, which is a Bohr almost-periodic function [Bohr1925] with frequency module $\mathbb Z n_b+\mathbb Z n_\oplus$. If $n_b/n_\oplus$ is irrational, $\mathbf r_{\rm geo}$ is **not periodic**, so it has no Fourier series in a single period. Example: for Mars $n_b/n_\oplus\approx0.5317$ (mean elements), and the synodic period is a derived quantity of the lattice, not an exact period. If $n_b/n_\oplus=p/q$ is rational, the motion is exactly periodic with period $2\pi q/n_b$, and (3.1) is an ordinary Fourier series. Motion is quasi-periodic on a 2-torus in the generic irrational case [Gallavotti2001].

### 3.3 The classical Ptolemaic model as the leading truncation
Keep only $k=1$ for both bodies, neglect $e$ ($c_1\to1$, $c_0,c_{-1}\to0$), and work in the ecliptic plane with complex notation:
$$z_{\rm geo}\approx a_be^{i\theta_b}-a_\oplus e^{i\theta_\oplus},\qquad\theta=nt+\text{const}.$$
- **Outer planet** (Mars, Jupiter, Saturn). The deferent is the body's $k=1$ circle, radius $a_b$, period $2\pi/n_b$. The epicycle is Earth's $k=1$ term, radius $a_\oplus$, period $2\pi/n_\oplus=1$ year. The epicycle-to-deferent ratio is $a_\oplus/a_b$: $0.656$ for Mars, versus Ptolemy's $39;30/60\approx0.658$ (Almagest values quoted from memory; cf. [Toomer1998; Evans1998; Neugebauer1975], to be verified at paper stage). The vector from the epicycle centre to the planet is $-\mathbf r_\oplus$, i.e. it is always parallel to the Earth-to-Sun direction: the classical rule that the epicycle radius stays parallel to the mean-Sun line.
- **Inner planet** (Mercury, Venus). The roles swap: the deferent is Earth's term with a sign flip, $-a_\oplus e^{i\theta_\oplus}$, i.e. the Sun's geocentric circle of period 1 year, and the epicycle is the planet's $k=1$ circle of radius $a_b$ and its own period. Ratio $a_b/a_\oplus$: $0.723$ for Venus, versus Ptolemy's $43;10/60\approx0.72$ (same caveat).
- By Apollonius's equivalence [Neugebauer1975; Evans1998], which is the algebraic commutativity of the sum $a_be^{i\theta_b}-a_\oplus e^{i\theta_\oplus}$, the two choices of which circle is deferent give the same geocentric curve. Which is called the deferent is a matter of nesting order.
- Sign and orientation: the minus sign on Earth's term is the half-turn phase shift of the epicycle (Sun's geocentric position is $-\mathbf r_\oplus$).

### 3.4 Eccentricity, the eccentric, and the equant
With the $k=0,1,2$ terms of Part 1 (first order in $e$),
$$\frac za\approx-\tfrac{3e}{2}+e^{iM}+\tfrac e2e^{2iM}.$$
**Equant.** Place the central body at $O=0$, the deferent centre at $C=-ea$, and the equant point at $Q=-2ea$ (so $C$ bisects $OQ$, Ptolemy's bisection). A point on the deferent circle of radius $a$ about $C$, whose direction from $Q$ advances uniformly at angle $M$, is $P=Q+\rho e^{iM}$ with $\rho\approx a+\mathrm{Re}(\bar d e^{iM})$, $d=C-Q=ea$. Therefore
$$P-O\approx-2ea+a e^{iM}+\tfrac{ea}{2}\big(e^{2iM}+1\big)=a\Big[-\tfrac{3e}{2}+e^{iM}+\tfrac e2e^{2iM}\Big]+O(e^2),$$
which coincides with the first-order Kepler series. So **the equant model reproduces the $k=0,1,2$ terms to first order in $e$**, and the missing pieces are $c_{-1}=O(e^2)$, $c_{-2}=O(e^3)$ and the $O(e^2)$ corrections to $c_1$, $c_2$. This is the known result that Ptolemy's equant (and the bisected eccentricity) captures Kepler's equation to $O(e)$ [Evans1984; Evans1988; CarmanRecio2019; Rushkin2015]. The $k=2$ circle is, in Ptolemy's language, the equant's contribution. Note the equant model is a different parametrisation of the same three leading terms; it is not itself a finite harmonic lattice expansion (lit review, folklore item b).

---

## Part 4. Truncation error

### 4.1 Single body
Let $z_N=a\sum_{|k|\le N}c_ke^{ikM}$ ($2N$ circles plus the offset). By (1.9),
$$\sup_{M}|z-z_N|\le a\sum_{k>N}(|c_k|+|c_{-k}|)\le 2a\Big[\frac1e+\rho\Big]\sum_{k>N}\frac{\rho^k}{k}\le\varepsilon_N(a,e)\equiv a\,\frac{2(1/e+\rho)\,\rho^{N+1}}{(N+1)(1-\rho)}. \tag{4.1}$$
The final step bounds $\sum_{k>N}\rho^k/k\le\frac1{N+1}\sum_{k>N}\rho^k=\frac{\rho^{N+1}}{(N+1)(1-\rho)}$. Since $z(M)$ is exactly $2\pi$-periodic in $M$, the supremum over $M$ equals the supremum over all $t$: **the bound holds uniformly for all time** for an unperturbed bound orbit. For $e\to0$ the bound is $O(\rho^{N+1}/e)\to0$ for $N\ge1$.

**Scaling.** $N_{\min}\sim\ln(1/\mathrm{tol})/\ln(1/\rho)$ up to logarithmic factors; and near $e\to1$, $\ln(1/\rho)\approx(2(1-e))^{3/2}/3$, so $N\propto(1-e)^{-3/2}$. High eccentricity is costly.

### 4.2 Geocentric
By (3.1) and the triangle inequality, keeping $|k|\le N_b$ for the body and $|k|\le N_\oplus$ for Earth, the 3D rotation $\mathsf R$ being an isometry,
$$\sup_t|\mathbf r_{\rm geo}-\mathbf r_{\rm geo,N}|\le\varepsilon_{N_b}(a_b,e_b)+\varepsilon_{N_\oplus}(a_\oplus,e_\oplus). \tag{4.2}$$
A position error $\varepsilon$ translates into an angular error $\lesssim\varepsilon/\Delta_{\min}$, where $\Delta_{\min}$ is the minimum Earth–body distance (e.g. about $0.37$ AU for Mars at a perihelic opposition, about $a_b-a_\oplus$ in general for an outer planet). Accuracy requirements in arcseconds must be set against $\Delta_{\min}$, not $a$.

### 4.3 Numerical results (scripts in the session scratchpad, not in the repo)
**(a) Coefficient-formula check.** Direct Newton solve of Kepler's equation versus the series (1.2)–(1.5), $20001$ points in $M\in[0,2\pi]$:

| $e$ | $N$ | $\max|z_N-z_{\rm true}|/a$ |
|---|---|---|
| $0.0167$ | $60$ | $9.6\times10^{-16}$ |
| $0.1$ | $60$ | $1.0\times10^{-15}$ |
| **$0.5$** | **$60$** | **$1.06\times10^{-14}$** |
| $0.9$ | $60$ | $3.8\times10^{-3}$ (truncation, not formula error; see table below) |

So the coefficient formulas are verified to $\sim10^{-14}$ at $e=0.5$, $N=60$.

**(b) Kapteyn inequality (1.7) and the coefficient bounds (1.8)** hold for all tested $k\le399$, $e\in\{0.0167,0.0934,0.25,0.6,0.89,0.99\}$ (to floating-point tolerance).

**(c) Number of harmonics $N$ needed** so that $\sup_M|z-z_N|/a<\mathrm{tol}$ (error on a grid of $3001$ points on $[0,\pi]$; the error is symmetric under $M\to-M$). "Observed" is the smallest $N$ with the true sup-error below tol; "bound" is the smallest $N$ with $\varepsilon_N/a$ below tol from (4.1). The bound was checked to hold for every $N$ from $0$ up to $800$ (6000 for $e=0.89$), with a $10^{-13}$ allowance for roundoff.

| $e$ | $\rho(e)$ | $N$ obs, tol $10^{-6}$ | $N$ obs, tol $10^{-12}$ | $N$ bound, $10^{-6}$ | $N$ bound, $10^{-12}$ |
|---|---|---|---|---|---|
| $0.0167$ (Earth) | $0.0227$ | $4$ | $7$ | $4$ | $8$ |
| $0.0934$ (Mars) | $0.1267$ | $6$ | $12$ | $7$ | $13$ |
| $0.25$ | $0.3345$ | $10$ | $22$ | $12$ | $24$ |
| $0.6$ | $0.7418$ | $33$ | $75$ | $43$ | $87$ |
| $0.89$ | $0.9644$ | $230$ | $575$ | $350$ | $712$ |

The bound overestimates $N$ by a factor of at most $\approx1.5$ in this range, so it is practically usable. Raw (untruncated-at-first-order) behaviour for small $N$: sup-error for $N=1$ is $8.5\times10^{-3}$ ($e=0.0167$), $5.0\times10^{-2}$ ($0.0934$), $0.148$ ($0.25$), $0.43$ ($0.6$), $0.73$ ($0.89$), in units of $a$.

**(d) Ptolemaic and truncated geocentric Mars.** Planar (coplanar) approximation, J2000 mean elements for Mars and Earth, mean motions as constants (Standish-type values), time span 60 years (28.1 synodic periods) at 0.5 d steps. Reference: the exact two-body (non-interacting Kepler) geocentric longitude. Maximum longitude error and number of retrograde intervals (true count: **28**):

| Truncation | max longitude error (deg) | RMS (deg) | retrograde loops |
|---|---|---|---|
| A: $k=1$ for both (circular deferent + epicycle) | $33.2$ | $9.7$ | $28$ |
| B: $k\in\{0,1\}$ both (eccentric deferent, circular epicycle) | $9.40$ | $2.7$ | $28$ |
| C: Mars $k\in\{0,1,2\}$ (equant level), Earth $k=1$ | $5.05$ | $0.97$ | $28$ |
| D: both $k\in\{0,1,2\}$ | $0.557$ | $0.159$ | $28$ |
| E: both $\lvert k\rvert\le2$ | $0.728$ | $0.182$ | $28$ |
| F: both $\lvert k\rvert\le5$ | $6\times10^{-4}$ | $1\times10^{-4}$ | $28$ |
| G: both $\lvert k\rvert\le10$ | $<10^{-4}$ | $<10^{-4}$ | $28$ |

Observations. (i) Even the circular N=1 model, A, produces all the retrograde loops (28 in 60 years, matching truth), which is the qualitative Ptolemaic achievement. (ii) Its worst longitude error is $\sim33^\circ$, reached near perihelic oppositions where $\Delta\approx0.37$ AU versus $0.52$ AU in the circular model, which is why Ptolemy needed the eccentric and equant. (iii) The error is not monotone in the truncation set for small sets (E is slightly worse than D in sup norm, because $c_{-1}$, $c_{-2}$ shift errors that previously partly cancelled); monotone decrease sets in once $N\gtrsim5$. (iv) The angular errors are consistent with (4.2) divided by $\Delta_{\min}$.

---

## Part 5. The perturbed N-body case

### 5.1 What survives from the Kepler case
The planets' mutual perturbations break the exactness of Parts 1–3 but, for the Solar System's regular part, not their form. In action–angle variables $(\boldsymbol I,\boldsymbol\theta)$ of a near-integrable Hamiltonian $H=H_0(\boldsymbol I)+\mu H_1(\boldsymbol I,\boldsymbol\theta)$:
- **Lagrange–Laplace secular theory** gives linear precession and node-regression eigenfrequencies $g_i$, $s_i$ in addition to the mean motions $n_j$ (the secular system is linear at first order). This replaces the fixed $(\omega,\Omega)$ of Part 2 by slowly rotating frames, i.e. extra frequencies $g_i,s_i$ in the spectrum [MurrayDermott1999 ch. 7; BrouwerClemence1961].
- **KAM theory** guarantees that for sufficiently small $\mu$ and non-resonant, Diophantine frequencies, a large measure of initial conditions lie on invariant tori; motion on them is quasi-periodic [Gallavotti2001 for the epicycle reading; Kolmogorov, Arnold, Moser, historically; these three are not in `refs.bib` and are flagged for the bibliography].
- On such a torus, with $d$ rationally independent frequencies $\boldsymbol\omega=(\omega_1,\dots,\omega_d)$, each Cartesian coordinate has the form
$$\mathbf r(t)=\sum_{\mathbf m\in\mathbb Z^d}\mathbf A_{\mathbf m}e^{i\langle\mathbf m,\boldsymbol\omega\rangle t},\qquad\sum_{\mathbf m}|\mathbf A_{\mathbf m}|<\infty,$$
for analytic Hamiltonians (absolute summability follows from analyticity of the torus embedding), a Bohr almost-periodic function [Bohr1925; Laskar1993; Laskar2003]. Each term is a uniformly rotating vector, i.e. an epicycle. The frequencies are **not** integer multiples of one fundamental frequency, so this is not a Fourier series in a period. They are integer combinations of $\boldsymbol\omega$, which for the Solar System are the $n_j$, $g_i$, $s_i$ and a few resonance-related arguments [Laskar1990; Laskar2003].
- The geocentric vector subtracts two such expansions, so the two-family structure (3.1) generalises: a union of the body's and Earth's quasi-periodic spectra, still free of cross terms in the Cartesian vector.

### 5.2 Constructing it: NAFF
Laskar's numerical analysis of fundamental frequencies [Laskar1990; Laskar1992; Laskar1993; Laskar1999] extracts $\{\mathbf A_{\mathbf m},\langle\mathbf m,\boldsymbol\omega\rangle\}$ from a numerically integrated or ephemeris time series $f(t)=\mathbf r(t)$ (one complex coordinate, or one per component) iteratively: maximise $\big|\langle f,e^{i\nu t}\rangle_\chi\big|$ over $\nu$ (with a window $\chi$, typically Hann), subtract the term found, repeat, with Gram-Schmidt orthogonalisation of successive basis functions. The analytic coefficients of Parts 1–3 provide a prior against which NAFF output can be validated (the planned demonstration on a real perturbed body, contribution (iii) of the lit review). Reported precision for quasi-periodic signals is far better than the Fourier resolution, with errors decreasing like a high power of $1/T$ for a Hann window [Laskar1993; check the exact exponent in the paper], while two lines closer than $\sim2\pi/T$ cannot be separated.

### 5.3 Limits (explicit)
1. **Chaos.** The Solar System's inner planets are chaotic with a Lyapunov time of about 5 Myr [Laskar1990]. Quasi-periodicity then holds only approximately on a finite window; frequencies drift slowly [FuLaskar2019]. "Arbitrary precision for all time" is false.
2. **Finite windows.** For any perturbed orbit the representation is of the form $\sum A_{\mathbf m}e^{i\langle\mathbf m,\boldsymbol\omega\rangle t}$ with slowly varying $A,\boldsymbol\omega$; for fixed $T$ it can be made accurate to the NAFF resolution but extrapolation outside the fitted window degrades with the drift rate.
3. **Frequency resolution.** Distinct lines cannot be resolved if separated by less than $\sim2\pi/T$. Frequency combinations that are very close (near-resonances) cannot be separated by a short ephemeris; the number of terms needed grows with $T$.
4. **Near-resonant/non-Diophantine regions.** KAM tori are destroyed near resonances; the measure of surviving tori $\to1$ only as $\mu\to0$. For the planets $\mu\sim10^{-3}$ and the Great Inequality (5:2 Jupiter-Saturn near-resonance) is a long-period term, not a breakdown, but it is the largest perturbation line.
5. **Unbound orbits.** For $e\ge1$ the orbit is not bounded and $z(M)$ is not periodic; the Fourier–Bessel series of Part 1 does not exist (hyperbolic Kepler motion has no $2\pi$-periodic $M$). Quasi-periodicity is a property of bounded motion. A comet on a hyperbolic flyby can be matched by trigonometric polynomials on any finite window (Weierstrass) but with no dynamical meaning and no extrapolation. Orbits with $e\to1^-$ converge, but with $N\propto(1-e)^{-3/2}$.
6. **Earth is not a Kepler point either.** Earth's orbital elements vary secularly (e.g. $e_\oplus$ between $\sim0.003$ and $\sim0.06$ over $10^5$ years) and the Earth-Moon barycentre is what moves on the Kepler-like orbit; the Earth's centre has an additional monthly term of amplitude $\approx4.67\times10^3$ km. The Sun is not at the barycentre either (displacement $\sim0.005$ AU due to Jupiter). All of these are additional frequency lines, not conceptual obstacles.

---

## Part 6. Claims to verify numerically (for pytest)

Write $z_N$ for the truncated series of Part 1 and $z_{\rm true}$ for the direct Newton solve of Kepler's equation (tolerance $10^{-14}$ in $E$). Unless stated, "error" means sup over a uniform grid of at least $10^4$ points in $M\in[0,2\pi]$.

1. **Coefficient formulas.** For $e=0.5$, $N=60$: $\max|z_N-z_{\rm true}|/a<10^{-12}$ (measured $1.06\times10^{-14}$). Also for $e\in\{0.0167,0.1\}$, $N=60$: $<10^{-12}$ (measured $\sim10^{-15}$).
2. **Derivative-free form.** Forms (1.2)/(1.3) and (1.4) agree to $10^{-14}$ for all $1\le k\le100$, $e\in\{0.0167,0.0934,0.25,0.6,0.89\}$; and the real form (1.6) agrees with $\alpha_k=c_k+c_{-k}$, $\beta_k=c_k-c_{-k}$.
3. **Periapsis sum rule.** $-\tfrac{3e}2+\sum_{k=1}^{K}\tfrac2kJ_k'(ke)\to1-e$ ($e=0.3$, $K=200$, tolerance $10^{-12}$).
4. **Circle at $e=0$.** $c_0=c_{-k}=0$ for $k\ge1$, $c_1=1$, $c_k=0$ for $k\ge2$ (exact); $z_N(M)=e^{iM}$.
5. **Small-$e$ coefficients.** $c_0=-3e/2$ exactly; $c_2/e\to\tfrac12$, $c_{-1}/e^2\to\tfrac18$, $1-c_1\approx e^2/2$ at $e=10^{-3}$, relative tolerance $10^{-3}$.
6. **Kapteyn bound.** $|J_k(ke)|\le\rho(e)^k$ and bounds (1.8) for $1\le k\le399$, $e\in\{0.0167,0.0934,0.25,0.6,0.89,0.99\}$ (allow $10^{-15}$ absolute slack).
7. **Truncation bound.** $\max_M|z-z_N|/a\le\varepsilon_N(1,e)+10^{-13}$ for all $0\le N\le400$ and $e\in\{0.0167,0.0934,0.25,0.6,0.89\}$ (use $N$ up to 1200 for $0.89$).
8. **N-needed table.** For the $e$ values above, the smallest $N$ with sup-error $<10^{-6}$ ($<10^{-12}$) is $4,6,10,33,230$ ($7,12,22,75,575$); the bound-based $N$ is $4,7,12,43,350$ ($8,13,24,87,712$). The measured $N$ is $\le$ bound-based $N$ in every case.
9. **3D rotation.** $\mathsf R$ from 2.1 is orthogonal with $\det=+1$ and $\mathbf P\times\mathbf Q=\mathbf N$; for random $(\Omega,i,\omega)$ and $e=0.5$, $N=60$, $\mathbf r_N=a c_0\mathbf P+\sum_{k=1}^N\mathrm{Re}(A_ke^{ikM})$ matches $\mathsf R(x,y,0)^{\mathsf T}$ to $10^{-12}$.
10. **Two-circle decomposition (2.2).** For random $k\le10$, $e=0.5$, the sum of the prograde and retrograde circles in (2.2) equals $\mathrm{Re}(A_ke^{ikM})$ to $10^{-14}$; the semi-axes of the point set are $a(|c_k|+|c_{-k}|)$ and $a\,||c_k|-|c_{-k}||$; every term lies in the plane through $\mathbf P,\mathbf Q$ (dot product with $\mathbf N$ is zero to $10^{-14}$).
11. **Geocentric sum and bound.** For Mars and Earth (elements of Part 4.3d, full 3D allowed), $\sup_t|\mathbf r_{\rm geo}-\mathbf r_{{\rm geo},N}|\le\varepsilon_N(a_b,e_b)+\varepsilon_N(a_\oplus,e_\oplus)$ for $N=1,\dots,10$.
12. **Ptolemaic N=1, Mars (planar test of 4.3d).** Over 60 years at 0.5 d steps: circular $k=1$ model A has exactly 28 retrograde intervals (equal to the exact two-body count) and max geocentric longitude error $\le34^\circ$ (measured $33.2^\circ$); model B (k in {0,1}) $\le9.5^\circ$ (measured $9.40^\circ$); model C (equant-level) $\le5.2^\circ$ (measured $5.05^\circ$); model F ($|k|\le5$) $\le0.001^\circ$ (measured $6\times10^{-4}$).
13. **Cartesian spectrum has no cross terms (NOT YET MEASURED; thresholds are proposals).** The FFT of $\mathbf r_{\rm geo}(t)$ (Mars − Earth, long window) has lines only near $kn_b$ and $jn_\oplus$, none near $kn_b\pm jn_\oplus$ above the leakage floor (use a Hann window and $\ge$ 200 yr; threshold $10^{-6}$ relative amplitude).
14. **Rational ratio gives periodicity (control; NOT YET MEASURED).** For $n_b/n_\oplus=1/2$ (synthetic), $\mathbf r_{\rm geo}(t+4\pi/n_\oplus)=\mathbf r_{\rm geo}(t)$ to $10^{-12}$; for an irrational ratio such as $1/\sqrt2$ the minimum over $t'\in[T,10T]$ of $\max_t|\mathbf r(t+t')-\mathbf r(t)|$ is bounded away from $0$ at the scale set by the amplitudes (qualitative).

---

## Part 7. Caveats and flagged issues (hypothesis qualifications)

1. **"Fourier series" is wrong for the geocentric motion.** One period does not exist unless $n_b/n_\oplus$ is rational. The correct statement is Bohr almost-periodic with frequency module $\mathbb Zn_b+\mathbb Zn_\oplus$ (Part 3.2). The heliocentric Kepler orbit alone *is* an honest Fourier series in $M$.
2. **"Arbitrary precision" holds only for bounded, unperturbed Kepler orbits, uniformly for all time.** For perturbed motion it holds on finite windows or on KAM tori with the stated resolution limit $\sim2\pi/T$; for chaotic motion only on finite time; for unbound ($e\ge1$) orbits the series does not exist (Part 5.3).
3. **The Kepler series is not new.** Bessel (1824), Kapteyn (1893) and Watson (1922/1944) have the expansion; Gallavotti (2001) has the KAM/epicycle reading; Peralta et al. (2003) have the numerical FFT pipeline (lit review). Parts 1 and 5.1 must be presented as background.
4. **Epicycles are not a fitting model in this construction.** Every radius and rate is fixed by $(a,e,i,\omega,\Omega,M_0)$ of body and Earth. This answers the "epicycles fit anything" objection [Carman2010] only for the dynamically derived expansion, not for Ptolemy's own free parameters.
5. **Ptolemaic models are not literally truncations of (3.1).** Ptolemy's equant and eccentric (Part 3.4) reproduce the $k=0,1,2$ terms only to first order in $e$ and are not harmonic-lattice expansions [Evans1984; Evans1988]. The pure N=1 circle model has errors of tens of degrees for Mars (Part 4.3d), so "Ptolemy ≈ leading Fourier truncation" is true qualitatively (retrograde loops) but not quantitatively (about $33^\circ$ at worst).
6. **Truncation error is not monotone** in the choice of retained terms at small $N$ (Part 4.3d: E vs D); the bound (4.1) is monotone, the realised longitude error need not be.
7. **Angular accuracy depends on $\Delta_{\min}$,** not on $a$: (4.2) bounds Cartesian error; the angular error scales as $1/\Delta$, and diverges for close approaches.
8. **High eccentricity is expensive:** $N\propto(1-e)^{-3/2}$; $e=0.89$ needs 575 harmonics ($\approx1150$ circles) for $10^{-12}$. Epicycle counts are not "a handful" except for planet-like $e$.
9. **Reductions assumed in Part 4.3(d):** coplanar planar approximation, mean-element Kepler orbits without perturbations, no light time, aberration, precession or nutation. These change the numbers quoted but not the structure of the argument. The numerical results are a test of the truncation machinery against an exact two-body reference, not an ephemeris comparison.
10. **Non-rigorous inputs.** (a) Kapteyn's inequality is quoted from Watson §8.7 and checked numerically, not re-proved here. (b) The Almagest ratios 39;30/60 and 43;10/60 are from memory. (c) The Hann-window NAFF error exponent is quoted from memory. (d) Kolmogorov/Arnold/Moser and Lagrange–Laplace have no entry in `refs.bib`; add if cited in the paper.
11. **Typo-prone point corrected during derivation.** The `x/a` coefficient is $\tfrac2kJ_k'(ke)$ (not $\tfrac2{k^2}$); the $k=0$ term is $-\tfrac{3e}2$. Both match the task brief and numerics.
