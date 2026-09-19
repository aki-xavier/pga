// pga — Projective Geometric Algebra Cl(3,0,1): 16-component multivectors,
// blades and motors for 3D Euclidean geometry (plane-based / dual PGA).
//
// Consumed by the simu crate as a sibling path dependency.  Its test suite
// (tests/multivector.rs) cross-checks the whole 16x16 product table against an
// independent reference implementation.
//
// P(R*_{3,0,1}) = Cl(3,0,1) — the "plane-based" (dual) PGA: the algebra of
// Euclidean 3D point/line/plane geometry.
//
// Basis generators {e1, e2, e3, e0} with e1^2 = e2^2 = e3^2 = +1, e0^2 = 0,
// generators anticommute pairwise.  Every basis blade is indexed by its
// bitmask (bit 0 = e1, bit 1 = e2, bit 2 = e3, bit 3 = e0), so the geometric
// product of two basis blades is a XOR/sign computation with an extra metric
// factor when both factors contain the null generator e0:
//
//   [0] grade 0: 1                        [8]  grade 1: e0    (plane at infinity)
//   [1][2][4]   grade 1: e1, e2, e3       [9][10][12] grade 2: e1e0, e2e0, e3e0 (ideal lines)
//   [3][5][6]   grade 2: e12, e13, e23    [11][13][14] grade 3: e12e0, e13e0, e23e0
//   [7] grade 3: e123 (the origin point)  [15] grade 4: I = e123e0 (null, I^2 = 0)
//
// Geometric semantics (the reason for the degenerate metric):
//   plane = grade-1 vector  n + d e0         (n unit normal, distance d)
//   line  = grade-2 blade   direction + moment (Euclidean, squares negative)
//   point = grade-3 blade   e123 + x e032 + y e013 + z e021
//   meet  = outer product (plane ^ plane = line, plane ^ line = point)
//   join  = regressive product (point v point = line, three points = plane)
//   motor = even-grade versor (rotation + translation), M X M~ applies it.
//
// The pseudoscalar I = e123 e0 squares to 0 (nilpotent), so the Hodge dual is
// NOT multiplication by I; it is the blade-wise complement below (PGA4CS
// Table 4 / bivector.net, converted to this basis order), and
// dual(dual(x)) = (-1)^grade(x) . x.

pub mod multivector;
pub mod primitives;

pub use multivector::*;
pub use primitives::*;
