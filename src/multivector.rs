// multivector.rs — the 16-component multivector, its operators, the Hodge dual,
// and the motor group (dual quaternions of SE(3)). Port of multivector.v.

use crate::primitives::point;
use std::fmt;

pub const NUM_COMPONENTS: usize = 16;
pub const NUM_GRADES: usize = 5;

/// Metric of the four generators in bit order (e1, e2, e3, e0).
pub const GENERATOR_METRIC: [f64; 4] = [1.0, 1.0, 1.0, 0.0];

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Multivector {
    pub values: [f64; NUM_COMPONENTS],
}

// --- constructors ------------------------------------------------------------

/// mv_zero returns the zero multivector.
pub fn mv_zero() -> Multivector {
    Multivector::default()
}

/// mv_scalar returns a scalar multivector (only the grade-0 component).
pub fn mv_scalar(s: f64) -> Multivector {
    let mut m = Multivector::default();
    m.values[0] = s;
    m
}

/// mv_vector builds a plane (grade-1) from coefficients (a, b, c, d): the
/// plane a x + b y + c z + d = 0.
pub fn mv_vector(a: f64, b: f64, c: f64, d: f64) -> Multivector {
    let mut m = Multivector::default();
    m.values[1] = a;
    m.values[2] = b;
    m.values[4] = c;
    m.values[8] = d;
    m
}

/// mv_bivector builds a bivector from the six components in the order
/// (e12, e13, e23, e1e0, e2e0, e3e0).
pub fn mv_bivector(b12: f64, b13: f64, b23: f64, b01: f64, b02: f64, b03: f64) -> Multivector {
    let mut m = Multivector::default();
    m.values[3] = b12;
    m.values[5] = b13;
    m.values[6] = b23;
    m.values[9] = b01;
    m.values[10] = b02;
    m.values[12] = b03;
    m
}

/// mv_trivector builds a trivector from the four components in the order
/// (e123, e12e0, e13e0, e23e0) — useful for raw point-like blades.
pub fn mv_trivector(t123: f64, t120: f64, t130: f64, t230: f64) -> Multivector {
    let mut m = Multivector::default();
    m.values[7] = t123;
    m.values[11] = t120;
    m.values[13] = t130;
    m.values[14] = t230;
    m
}

/// pseudoscalar returns I = e123 e0 (squares to 0).
pub fn pseudoscalar() -> Multivector {
    let mut m = Multivector::default();
    m.values[15] = 1.0;
    m
}

/// e123 returns I3 = e1 e2 e3, the Euclidean pseudoscalar (squares to -1).
pub fn e123() -> Multivector {
    let mut m = Multivector::default();
    m.values[7] = 1.0;
    m
}

/// e1 is a basis vector (grade 1). Fresh value per call.
pub fn e1() -> Multivector {
    mv_vector(1.0, 0.0, 0.0, 0.0)
}

pub fn e2() -> Multivector {
    mv_vector(0.0, 1.0, 0.0, 0.0)
}

pub fn e3() -> Multivector {
    mv_vector(0.0, 0.0, 1.0, 0.0)
}

/// e0 is the plane at infinity.
pub fn e0() -> Multivector {
    mv_vector(0.0, 0.0, 0.0, 1.0)
}

// --- component access --------------------------------------------------------

/// popcount counts the set bits of a small non-negative integer.
///
/// Public only because the ported test suite (which V keeps inside the module)
/// uses it.
pub fn popcount(x: u32) -> u32 {
    let mut n = 0;
    let mut v = x;
    while v > 0 {
        n += v & 1;
        v >>= 1;
    }
    n
}

impl Multivector {
    /// grade returns the grade-g projection.
    pub fn grade(&self, g: usize) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            if popcount(i as u32) as usize == g {
                res.values[i] = self.values[i];
            }
        }
        res
    }

    /// scalar_part returns the grade-0 component.
    pub fn scalar_part(&self) -> f64 {
        self.values[0]
    }

    /// vector_part returns the four grade-1 components (e1, e2, e3, e0).
    pub fn vector_part(&self) -> [f64; 4] {
        [
            self.values[1],
            self.values[2],
            self.values[4],
            self.values[8],
        ]
    }

    /// bivector_part returns the six grade-2 components in the order
    /// (e12, e13, e23, e1e0, e2e0, e3e0).
    pub fn bivector_part(&self) -> [f64; 6] {
        [
            self.values[3],
            self.values[5],
            self.values[6],
            self.values[9],
            self.values[10],
            self.values[12],
        ]
    }

    /// pseudoscalar_part returns the grade-4 component.
    pub fn pseudoscalar_part(&self) -> f64 {
        self.values[15]
    }

    /// is_zero reports whether all components are approximately zero.
    pub fn is_zero(&self) -> bool {
        for v in self.values {
            if v.abs() > 1e-10 {
                return false;
            }
        }
        true
    }

    /// vmax returns the maximum absolute component.
    pub fn vmax(&self) -> f64 {
        let mut mx = 0.0;
        for v in self.values {
            let a = v.abs();
            if a > mx {
                mx = a;
            }
        }
        mx
    }

    // --- operators -----------------------------------------------------------

    pub fn add(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            res.values[i] = self.values[i] + o.values[i];
        }
        res
    }

    pub fn sub(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            res.values[i] = self.values[i] - o.values[i];
        }
        res
    }

    pub fn scale(&self, s: f64) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            res.values[i] = self.values[i] * s;
        }
        res
    }

    pub fn div_scalar(&self, s: f64) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            res.values[i] = self.values[i] / s;
        }
        res
    }

    pub fn neg(&self) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            res.values[i] = -self.values[i];
        }
        res
    }

    /// approx_eq performs approximate equality (allclose, atol=1e-6).
    ///
    /// V calls this `eq`; the name is taken by the language's own comparison
    /// trait here, and `PartialEq` is the exact comparison.
    pub fn approx_eq(&self, o: Multivector) -> bool {
        for i in 0..NUM_COMPONENTS {
            if (self.values[i] - o.values[i]).abs() > 1e-6 {
                return false;
            }
        }
        true
    }

    /// copy returns a fresh multivector with the same components. (Rust values
    /// are copied by assignment; this exists for the port's call sites.)
    pub fn copy(&self) -> Multivector {
        *self
    }

    /// str renders non-zero components by grade.
    pub fn str(&self) -> String {
        let mut parts: Vec<String> = Vec::new();
        for g in 0..NUM_GRADES {
            for i in 0..NUM_COMPONENTS {
                if popcount(i as u32) as usize != g {
                    continue;
                }
                let v = self.values[i];
                if v.abs() <= 1e-10 {
                    continue;
                }
                let name = blade_name(i);
                if name == "1" {
                    parts.push(format!("{v:.4}"));
                } else {
                    parts.push(format!("{v:+.4}*{name}"));
                }
            }
        }
        if parts.is_empty() {
            return "Multivector(0)".to_string();
        }
        format!("Multivector({})", parts.join(" "))
    }

    // --- algebra -------------------------------------------------------------

    /// gp computes the geometric product.  Basis-blade masks multiply by XOR;
    /// the sign counts inversions over the generator list; a shared e0 factor
    /// makes the term vanish because e0^2 = 0.
    pub fn gp(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for ma in 0..NUM_COMPONENTS {
            let a = self.values[ma];
            if a == 0.0 {
                continue;
            }
            for mb in 0..NUM_COMPONENTS {
                let b = o.values[mb];
                if b == 0.0 {
                    continue;
                }
                let (dst, coeff) = gp_blade(ma as u32, mb as u32);
                if coeff == 0.0 {
                    continue;
                }
                res.values[dst] += coeff * a * b;
            }
        }
        res
    }

    /// op computes the outer product a ^ b (the PGA meet): zero for blades
    /// sharing a basis factor, else +/- the union blade.
    pub fn op(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for ma in 0..NUM_COMPONENTS {
            let a = self.values[ma];
            if a == 0.0 {
                continue;
            }
            for mb in 0..NUM_COMPONENTS {
                let b = o.values[mb];
                if b == 0.0 {
                    continue;
                }
                if ma & mb != 0 {
                    continue;
                }
                let (dst, coeff) = gp_blade(ma as u32, mb as u32);
                res.values[dst] += coeff * a * b;
            }
        }
        res
    }

    /// ip computes the Hestenes fat-dot inner product.
    ///
    /// A|B = sum_{r,s>=1} < <A>_r <B>_s >_|r-s| ; scalar (grade-0) terms are zero.
    pub fn ip(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for ga in 1..NUM_GRADES {
            let a_g = self.grade(ga);
            if a_g.is_zero() {
                continue;
            }
            for gb in 1..NUM_GRADES {
                let b_g = o.grade(gb);
                if b_g.is_zero() {
                    continue;
                }
                let prod = a_g.gp(b_g);
                res = res.add(prod.grade((gb as i64 - ga as i64).unsigned_abs() as usize));
            }
        }
        res
    }

    /// lc is the left contraction A _| B: sums of <_A_g _B_h>_(h-g) for g <= h.
    pub fn lc(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for ga in 1..NUM_GRADES {
            let a_g = self.grade(ga);
            if a_g.is_zero() {
                continue;
            }
            for gb in ga..NUM_GRADES {
                let b_g = o.grade(gb);
                if b_g.is_zero() {
                    continue;
                }
                res = res.add(a_g.gp(b_g).grade(gb - ga));
            }
        }
        res
    }

    /// rc is the right contraction A |_ B: sums of <_A_g _B_h>_(g-h) for g >= h.
    pub fn rc(&self, o: Multivector) -> Multivector {
        let mut res = Multivector::default();
        for ga in 1..NUM_GRADES {
            let a_g = self.grade(ga);
            if a_g.is_zero() {
                continue;
            }
            for gb in 1..=ga {
                let b_g = o.grade(gb);
                if b_g.is_zero() {
                    continue;
                }
                res = res.add(a_g.gp(b_g).grade(ga - gb));
            }
        }
        res
    }

    /// reverse applies the reversal involution: grade-k blade * (-1)^(k(k-1)/2).
    pub fn reverse(&self) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            // signed on purpose: V's int arithmetic gives k(k-1)/2 = 0 for the
            // scalar component, where an unsigned k - 1 would underflow
            let k = popcount(i as u32) as i32;
            if (k * (k - 1) / 2) % 2 != 0 {
                res.values[i] = -self.values[i];
            } else {
                res.values[i] = self.values[i];
            }
        }
        res
    }

    /// grade_involution flips the sign of odd-grade components.
    pub fn grade_involution(&self) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            if !popcount(i as u32).is_multiple_of(2) {
                res.values[i] = -self.values[i];
            } else {
                res.values[i] = self.values[i];
            }
        }
        res
    }

    /// conjugate is the Clifford conjugate (reverse + grade involution).
    pub fn conjugate(&self) -> Multivector {
        self.reverse().grade_involution()
    }

    // --- Hodge dual and regressive joins -------------------------------------

    /// dual is the PGA Hodge / Poincare complement: blade-wise with the signs of
    /// PGA4CS Table 4 (bivector.net), converted to the (e1,e2,e3,e0) bit order.
    /// dual(dual(x)) = (-1)^grade(x) . x.
    pub fn dual(&self) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            if self.values[i] != 0.0 {
                res.values[DUAL_DST[i]] += DUAL_SIGN[i] * self.values[i];
            }
        }
        res
    }

    /// undual returns (-1)^grade(x) . dual(x), the inverse of dual.
    pub fn undual(&self) -> Multivector {
        let mut res = Multivector::default();
        for i in 0..NUM_COMPONENTS {
            if self.values[i] != 0.0 {
                let sgn = DUAL_SIGN[i];
                if !popcount(DUAL_DST[i] as u32).is_multiple_of(2) {
                    // (-1)^grade of the result: odd grades flip
                    res.values[DUAL_DST[i]] += -sgn * self.values[i];
                } else {
                    res.values[DUAL_DST[i]] += sgn * self.values[i];
                }
            }
        }
        res
    }

    /// meet returns the intersection of two blades: the outer product (native
    /// PGA incidence), e.g. meet(plane, plane) = line, meet(plane, line) = point.
    pub fn meet(&self, o: Multivector) -> Multivector {
        self.op(o)
    }

    /// join returns the union of two blades via the regressive product
    /// (self* ^ o*)*, e.g. join(point, point) = line, join(line, point) = plane.
    pub fn join(&self, o: Multivector) -> Multivector {
        self.dual().op(o.dual()).dual()
    }

    // --- norms and inverses --------------------------------------------------

    /// norm returns sqrt(|<self.reverse(self)>_0|); lines and rotors have norm 1,
    /// null blades (ideal lines, points) have norm 0.
    pub fn norm(&self) -> f64 {
        let s = self.gp(self.reverse()).values[0];
        s.abs().sqrt()
    }

    /// normalized returns the unit-norm multivector (or zero if degenerate).
    pub fn normalized(&self) -> Multivector {
        let n = self.norm();
        if n < 1e-12 {
            return mv_zero();
        }
        self.div_scalar(n)
    }

    /// inverse returns A^-1 = rev(A) (s - p I)/s^2 with A rev(A) = s + p I — exact
    /// for blades and motors (unit elements give rev(A)); null elements (points,
    /// ideal lines) panic.
    pub fn inverse(&self) -> Multivector {
        let prod = self.gp(self.reverse());
        let s = prod.values[0];
        let p = prod.values[15];
        if s.abs() < 1e-12 {
            panic!("inverse: degenerate null element");
        }
        // (s + p I)^-1 = (s - p I) / s^2 since I^2 = 0
        let inv = self
            .reverse()
            .scale(s)
            .sub(self.reverse().gp(pseudoscalar()).scale(p));
        inv.div_scalar(s * s)
    }

    // --- motion (the even subalgebra: dual quaternions of SE(3)) -------------

    /// apply returns M . v . M~ for a versor M; rotates/translates points, lines
    /// and planes.
    pub fn apply(&self, v: Multivector) -> Multivector {
        self.gp(v).gp(self.reverse())
    }

    /// exp exponentiates a scalar + bivector to the motor group: exp(s + B) =
    /// e^s exp(B).  Simple Euclidean bivectors use cos/sin, ideal (nilpotent)
    /// bivectors truncate at 1 + B, and general screw bivectors use the axis-pair
    /// closed form.
    pub fn exp(&self) -> Multivector {
        if !self.grade(1).is_zero() || !self.grade(3).is_zero() || !self.grade(4).is_zero() {
            panic!("exp: only scalar + bivector are supported");
        }
        let scale = self.values[0].exp();
        let b = self.grade(2);
        if b.is_zero() {
            return mv_scalar(scale);
        }
        let b2 = b.gp(b);
        let s = b2.values[0];
        if s.abs() < 1e-12 {
            // nilpotent: B^2 = 0 (pure translation / ideal line)
            return mv_scalar(scale).add(b.scale(scale));
        }
        let p = b2.values[15];
        let u = (-s).sqrt();
        let v = -p / (2.0 * u);
        let perp = b.gp(pseudoscalar()).div_scalar(u); // theta_hat _| = B.I / u
        let hat = b.sub(perp.scale(v)).div_scalar(u);
        let cu = scale * u.cos();
        let su = scale * u.sin();
        mv_scalar(cu)
            .add(hat.scale(su))
            .gp(mv_scalar(1.0).add(perp.scale(v)))
    }

    /// log returns the bivector B with exp(B) = self for a unit motor.  Pure
    /// translations return their ideal line; the -1 motor panics.
    pub fn log(&self) -> Multivector {
        let a = self.values[0];
        let b = self.grade(2);
        let c = self.values[15];
        let sin_u = (1.0 - a * a).max(0.0).sqrt();
        if sin_u < 1e-9 {
            if b.is_zero() {
                if a < 0.0 {
                    panic!("log: the -1 motor has no unique logarithm");
                }
                return mv_zero();
            }
            return b; // pure translation
        }
        let u = sin_u.atan2(a);
        let v = -c / sin_u;
        let perp = b.gp(pseudoscalar()).div_scalar(sin_u); // sin_u . theta_hat _| = B.I
        let beta = v * u.cos();
        let hat = b.sub(perp.scale(beta)).div_scalar(sin_u);
        hat.scale(u).add(perp.scale(v))
    }

    /// to_matrix returns the equivalent 4x4 homogeneous transform [R|t], flattened
    /// row-major into 16 components (row r, col c is at index 4*r + c).
    pub fn to_matrix(&self) -> [f64; 16] {
        let origin_t = self.apply(point(0.0, 0.0, 0.0));
        let px = self.apply(point(1.0, 0.0, 0.0)).coords();
        let py = self.apply(point(0.0, 1.0, 0.0)).coords();
        let pz = self.apply(point(0.0, 0.0, 1.0)).coords();
        let tx = origin_t.coords();
        [
            px[0] - tx[0],
            py[0] - tx[0],
            pz[0] - tx[0],
            tx[0],
            px[1] - tx[1],
            py[1] - tx[1],
            pz[1] - tx[1],
            tx[1],
            px[2] - tx[2],
            py[2] - tx[2],
            pz[2] - tx[2],
            tx[2],
            0.0,
            0.0,
            0.0,
            1.0,
        ]
    }

    /// blade_grade returns the grade of a pure blade, or -1 when the multivector
    /// holds mixed grades or is zero.
    ///
    /// V keeps this helper module-private and unused; it is kept here for parity.
    #[allow(dead_code)]
    fn blade_grade(&self) -> i32 {
        let mut g: i32 = -1;
        for i in 0..NUM_COMPONENTS {
            if self.values[i] != 0.0 {
                let gi = popcount(i as u32) as i32;
                if g != -1 && gi != g {
                    return -1;
                }
                g = gi;
            }
        }
        g
    }
}

impl fmt::Display for Multivector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.str())
    }
}

/// gp_blade returns the product of two basis blades as (result mask, coeff).
///
/// Public only because the ported test suite (which V keeps inside the module)
/// cross-checks it against its own reference implementation.
pub fn gp_blade(ma: u32, mb: u32) -> (usize, f64) {
    let mut swaps = 0u32;
    for i in 0..4 {
        if mb & (1 << i) != 0 {
            swaps += popcount(ma >> (i + 1));
        }
    }
    // common basis factors square to their metric; the null generator e0
    // squares to 0, so any shared e0 factor annihilates the term.
    let mut coeff = 1.0;
    for (i, m) in GENERATOR_METRIC.iter().enumerate() {
        if (ma & mb) & (1 << i) != 0 {
            coeff *= m;
        }
    }
    let mut sign = 1.0;
    if !swaps.is_multiple_of(2) {
        sign = -1.0;
    }
    ((ma ^ mb) as usize, sign * coeff)
}

/// rotor returns the rotation rotor by `angle` radians about the axis through
/// the origin in direction `axis`: R = exp(-angle/2 . L) with L = unit line.
pub fn rotor(axis: [f64; 3], angle: f64) -> Multivector {
    let len2 = axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2];
    if len2 < 1e-18 {
        panic!("rotor: zero axis");
    }
    let inv_len = 1.0 / len2.sqrt();
    let u = mv_vector(axis[0] * inv_len, axis[1] * inv_len, axis[2] * inv_len, 0.0);
    // unit line through the origin: u _| I3
    let line = u.lc(e123());
    let half = angle / 2.0;
    mv_scalar(half.cos()).sub(line.scale(half.sin()))
}

/// translator returns the translator T = exp(-d/2 . e0 ^ t) = 1 - (e0 ^ t)/2
/// (the ideal line e0 ^ t is nilpotent, so the series truncates).
pub fn translator(displacement: [f64; 3]) -> Multivector {
    let t = mv_vector(displacement[0], displacement[1], displacement[2], 0.0);
    mv_scalar(1.0).sub(e0().op(t).scale(0.5))
}

/// motor builds the rigid motion M = T . R (rotate then translate).
pub fn motor(axis: [f64; 3], angle: f64, displacement: [f64; 3]) -> Multivector {
    translator(displacement).gp(rotor(axis, angle))
}

/// motor_identity returns the identity motor.
pub fn motor_identity() -> Multivector {
    mv_scalar(1.0)
}

/// interpolate slerps between motors m1 and m2: M(t) = m1 . exp(t . log(m1~ . m2)).
pub fn interpolate(m1: Multivector, m2: Multivector, t: f64) -> Multivector {
    let rel = m1.reverse().gp(m2);
    m1.gp(rel.log().scale(t).exp())
}

/// blade_name returns the symbol of a basis blade, sorted by generator
/// (e1, e2, e3, e0).
fn blade_name(i: usize) -> &'static str {
    match i {
        0 => "1",
        1 => "e1",
        2 => "e2",
        4 => "e3",
        8 => "e0",
        3 => "e12",
        5 => "e13",
        6 => "e23",
        9 => "e1e0",
        10 => "e2e0",
        12 => "e3e0",
        7 => "e123",
        11 => "e12e0",
        13 => "e13e0",
        14 => "e23e0",
        15 => "I",
        _ => "?",
    }
}

/// The Hodge complement: blade -> (complement slot, sign).  Converted from
/// PGA4CS Table 4 / bivector.net: 1<->I, e0<->e123, e1<->e032, e2<->e013,
/// e3<->e021, e01<->e23, e02<->e31, e03<->e12, e032->-e1, e013->-e2,
/// e021->-e3, e123->-e0.
const DUAL_DST: [usize; 16] = [
    15, // 1    -> I
    14, // e1   -> e032 (= e23e0 slot)
    13, // e2   -> e013 (= e13e0 slot)
    12, // e12  -> e03 (= e3e0 slot)
    11, // e3   -> e021 (= e12e0 slot)
    10, // e13  -> -e02 (= e2e0 slot)
    9,  // e23  -> e01 (= e1e0 slot)
    8,  // e123 -> -e0
    7,  // e0   -> e123
    6,  // e01  -> e23
    5,  // e02  -> -e31 (= e13 slot, negated)
    4,  // e021 -> -e3
    3,  // e03  -> e12
    2,  // e013 -> -e2
    1,  // e032 -> -e1
    0,  // I    -> 1
];

const DUAL_SIGN: [f64; 16] = [
    1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
];
