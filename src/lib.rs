// pga — Projective Geometric Algebra Cl(3,0,1): 16-component multivectors,
// blades and motors for 3D Euclidean geometry (plane-based / dual PGA).
//
// P(R*_{3,0,1}) = Cl(3,0,1), with generators {e1, e2, e3, e0} where e1^2 = e2^2 =
// e3^2 = +1 and e0^2 = 0.  The degenerate metric is what makes the translation
// half of a rigid motion multiplicative rather than affine, so a motor is an
// even-grade versor and applying it is one product, M X M~.
//
// Basis blades are indexed by bitmask (bit 0 = e1, bit 1 = e2, bit 2 = e3, bit 3 =
// e0), so a basis product is a XOR plus an inversion count — see gp_blade.
//
// I = e123 e0 is nilpotent (I^2 = 0), so the Hodge dual is NOT multiplication by
// I: it is the blade-wise complement DUAL_DST, converted into this basis order
// from PGA4CS Table 4 (bivector.net), and dual(dual(x)) = (-1)^grade(x) . x.

pub mod multivector;
pub mod primitives;

pub use multivector::*;
pub use primitives::*;
