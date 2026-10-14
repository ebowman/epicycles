/** Element sets ported from epicycles/bodies.py (identical numbers). */
import { Elements, MU_SUN } from "./kepler";

export const J2000 = 2451545.0;
const rad = (d: number): number => (d * Math.PI) / 180;

function planet(name: string, a: number, e: number, I: number, L: number,
                varpi: number, node: number, epoch = J2000): Elements {
  // Python's (L - varpi) % 360.0 has the sign of the divisor
  let m = (L - varpi) % 360.0;
  if (m < 0) m += 360.0;
  return { a, e, i: rad(I), Omega: rad(node), omega: rad(varpi - node), M0: rad(m),
           epoch, mu: MU_SUN, name };
}

function asteroid(name: string, a: number, e: number, I: number, node: number,
                  w: number, ma: number, epoch = 2461200.5): Elements {
  return { a, e, i: rad(I), Omega: rad(node), omega: rad(w), M0: rad(ma),
           epoch, mu: MU_SUN, name };
}

export const EARTH = planet("Earth", 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0);
export const MARS = planet("Mars", 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891);
export const JUPITER = planet("Jupiter", 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909);
export const SATURN = planet("Saturn", 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448);

export const CERES = asteroid("Ceres", 2.765552595034094, 0.07969229514816586, 10.58802780183462,
                              80.24862682043221, 73.29421453021587, 274.4193463761342);
export const EROS = asteroid("Eros", 1.458243716760167, 0.2228779627700761, 10.82854410314273,
                             304.2679713350896, 178.9181319135911, 62.51145501986792);
export const PHAETHON = asteroid("Phaethon", 1.271464620755913, 0.8896722844020535, 22.31052728999002,
                                 265.0988060323288, 322.3001684957205, 301.4858238742286);

export const BODIES: Record<string, Elements> = {
  EARTH, MARS, JUPITER, SATURN, CERES, EROS, PHAETHON,
};
