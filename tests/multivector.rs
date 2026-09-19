// The geometric product is cross-checked against an independent reference
// (ref_blade_gp) that folds blades back into generators and bubble-sorts them with
// the generator metric [1, 1, 1, 0] (e0 squares to 0) — so a wrong sign or a wrong
// slot in the 16x16 table is caught against a different derivation rather than
// against itself.

use pga::*;

fn almost(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn mv_almost(a: Multivector, b: Multivector) -> bool {
    for i in 0..16 {
        if (a.values[i] - b.values[i]).abs() > 1e-9 {
            return false;
        }
    }
    true
}

/// Both blades expand into generator index lists, are bubble-sorted (each inversion
/// flips the sign), and each duplicated generator folds into its metric coefficient.
fn ref_blade_gp(ma: u32, mb: u32) -> (usize, f64) {
    let mut gens: Vec<u32> = Vec::new();
    for i in 0..4 {
        if ma & (1 << i) != 0 {
            gens.push(i);
        }
    }
    for i in 0..4 {
        if mb & (1 << i) != 0 {
            gens.push(i);
        }
    }
    let mut sign = 1.0;
    let n = gens.len();
    for i in 0..n.saturating_sub(1) {
        for j in 0..n - 1 - i {
            if gens[j] > gens[j + 1] {
                gens.swap(j, j + 1);
                sign = -sign;
            }
        }
    }
    let mut coeff = 1.0;
    let mut mask = 0u32;
    for g in gens {
        if mask & (1 << g) != 0 {
            coeff *= GENERATOR_METRIC[g as usize];
            mask &= !(1 << g);
        } else {
            mask |= 1 << g;
        }
    }
    (mask as usize, sign * coeff)
}

fn blade(mask: usize) -> Multivector {
    let mut m = Multivector::default();
    m.values[mask] = 1.0;
    m
}

#[test]
fn basis_metric() {
    assert!(mv_almost(e1().gp(e1()), mv_scalar(1.0)));
    assert!(mv_almost(e2().gp(e2()), mv_scalar(1.0)));
    assert!(mv_almost(e3().gp(e3()), mv_scalar(1.0)));
    assert!(e0().gp(e0()).is_zero());
    // euclidean bivectors square to -1, ideal lines are null
    assert!(mv_almost(
        mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0).gp(mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        mv_scalar(-1.0)
    ));
    assert!(mv_bivector(0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
        .gp(mv_bivector(0.0, 0.0, 0.0, 1.0, 0.0, 0.0))
        .is_zero());
    // I3 squares to -1, the pseudoscalar I squares to 0
    assert!(mv_almost(e123().gp(e123()), mv_scalar(-1.0)));
    assert!(mv_almost(pseudoscalar().gp(pseudoscalar()), mv_zero()));
}

#[test]
fn gp_reference_crosscheck() {
    for ma in 0..16u32 {
        for mb in 0..16u32 {
            let (dst, coeff) = gp_blade(ma, mb);
            let (rdst, rcoeff) = ref_blade_gp(ma, mb);
            assert_eq!(dst, rdst, "blade product mask for {ma}x{mb}");
            assert!(
                almost(coeff, rcoeff),
                "blade product coefficient for {ma}x{mb}: {coeff} vs {rcoeff}"
            );
        }
    }
}

#[test]
fn gp_properties() {
    let a = mv_scalar(2.0)
        .add(e1())
        .add(mv_bivector(0.5, -1.0, 0.25, 0.1, -0.3, 0.2));
    let b = mv_vector(1.0, -2.0, 0.5, 0.25).add(mv_bivector(0.0, 0.25, -0.5, 0.2, 0.1, 0.4));
    let c = e3().add(pseudoscalar()).add(mv_scalar(-1.5));
    assert!(mv_almost(a.add(b).gp(c), a.gp(c).add(b.gp(c))));
    assert!(mv_almost(a.gp(b).gp(c), a.gp(b.gp(c))));
}

#[test]
fn outer_product() {
    assert!(e1().op(e1()).is_zero());
    assert!(mv_almost(e1().op(e2()).op(e3()), e123()));
    assert!(mv_almost(e1().op(e2()).op(e3()).op(e0()), pseudoscalar()));
    assert!(e1().op(mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)).is_zero());
}

#[test]
fn reverse_involution() {
    assert!(mv_almost(e0().reverse(), e0()));
    assert!(mv_almost(e1().reverse(), e1()));
    assert!(mv_almost(
        mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0).reverse(),
        mv_bivector(-1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ));
    assert!(mv_almost(e123().reverse(), e123().neg()));
    assert!(mv_almost(pseudoscalar().reverse(), pseudoscalar()));
    assert!(mv_almost(e1().grade_involution(), e1().neg()));
    assert!(mv_almost(
        mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0).grade_involution(),
        mv_bivector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ));
}

#[test]
fn dual() {
    // spot values from the bivector.net table converted to this basis order
    assert!(mv_almost(e1().dual(), e2().op(e3()).op(e0())));
    assert!(mv_almost(e2().dual(), e1().op(e3()).op(e0())));
    assert!(mv_almost(e3().dual(), e1().op(e2()).op(e0())));
    assert!(mv_almost(e0().dual(), e123()));
    assert!(mv_almost(e123().dual(), e0().neg()));
    assert!(mv_almost(e1().op(e2()).dual(), e3().op(e0())));
    assert!(mv_almost(pseudoscalar().dual(), mv_scalar(1.0)));
    assert!(mv_almost(mv_scalar(1.0).dual(), pseudoscalar()));
    // dual^2 = (-1)^grade x and undual inverts it
    for i in 0..16 {
        let x = blade(i);
        let dd = x.dual().dual();
        if popcount(i as u32).is_multiple_of(2) {
            assert!(mv_almost(dd, x), "dual^2 for blade {i}");
        } else {
            assert!(mv_almost(dd, x.neg()), "dual^2 for blade {i}");
        }
        assert!(
            mv_almost(x.undual().dual(), x),
            "undual inverts dual for {i}"
        );
    }
}

#[test]
fn meet() {
    // planes meet in lines, plane/lines meet in points
    assert!(mv_almost(
        plane([1.0, 0.0, 0.0], 0.0).meet(plane([0.0, 1.0, 0.0], 0.0)),
        e1().op(e2())
    ));
    // point ^ plane gives the oriented (n.x + d) . I
    assert!(almost(point(2.0, 0.0, 0.0).meet(e1()).values[15], 2.0));
    assert!(almost(point(0.0, 0.0, 1.0).meet(e3()).values[15], 1.0));
    assert!(almost(point(0.0, 0.0, 0.0).meet(e1()).values[15], 0.0));
    assert!(mv_almost(e1().meet(e1()), mv_zero()));
}

#[test]
fn join() {
    // the join of two points is a line through them
    let l = point(0.0, 0.0, 0.0).join(point(1.0, 0.0, 0.0));
    assert!(almost(point_line_dist(point(1.0, 0.0, 0.0), l), 0.0));
    assert!(almost(point_line_dist(point(2.0, 0.0, 0.0), l), 0.0));
    assert!(almost(point_line_dist(point(0.0, 1.0, 0.0), l), 1.0));
}

#[test]
fn rotor_translator_motor() {
    let r = rotor([0.0, 0.0, 1.0], std::f64::consts::PI / 2.0);
    let got = r.apply(point(1.0, 0.0, 0.0)).coords();
    assert!(got[0].abs() < 1e-6 && (got[1] - 1.0).abs() < 1e-6);
    // rotation about x sends (0,0,1) to (0,-1,0)
    let rx = rotor([1.0, 0.0, 0.0], std::f64::consts::PI / 2.0);
    let got2 = rx.apply(point(0.0, 0.0, 1.0)).coords();
    assert!(got2[0].abs() < 1e-6 && (got2[1] + 1.0).abs() < 1e-6);
    let tt = translator([1.0, 0.0, 0.0])
        .apply(point(2.0, 1.0, 0.0))
        .coords();
    assert!((tt[0] - 3.0).abs() < 1e-6 && (tt[1] - 1.0).abs() < 1e-6);
    // motor: rotate about x then translate z
    let m = motor([1.0, 0.0, 0.0], std::f64::consts::PI / 2.0, [0.0, 0.0, 1.0]);
    let mc = m.apply(point(1.0, 0.0, 0.0)).coords();
    assert!((mc[0] - 1.0).abs() < 1e-6 && (mc[2] - 1.0).abs() < 1e-6);
    // motors are unit: m m~ has scalar 1
    assert!(almost(m.gp(m.reverse()).values[0], 1.0));
}

#[test]
fn exp_log_roundtrip() {
    let b = mv_bivector(0.5, -0.2, 0.3, 0.1, 0.4, -0.2);
    assert!(mv_almost(b.exp().log(), b), "screw bivector");
    let br = mv_bivector(0.3, 0.0, 0.0, 0.0, 0.0, 0.0);
    assert!(mv_almost(br.exp().log(), br), "pure rotation");
    let bt = mv_bivector(0.0, 0.0, 0.0, 0.4, -0.2, 0.3);
    assert!(mv_almost(bt.exp().log(), bt), "pure translation");
    // exp of a translator's generator equals the translator itself
    let t = translator([1.0, 0.0, 0.0]);
    assert!(mv_almost(
        mv_bivector(0.0, 0.0, 0.0, 0.5, 0.0, 0.0).exp(),
        t
    ));
}

#[test]
fn inverse() {
    let pi = plane([0.3, -0.4, 0.5], 1.0);
    assert!(mv_almost(pi.gp(pi.inverse()), mv_scalar(1.0)));
    let l = line_from_points(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0));
    assert!(mv_almost(l.gp(l.inverse()), mv_scalar(1.0)));
    let r = rotor([1.0, 2.0, 3.0], 1.234);
    assert!(mv_almost(r.gp(r.inverse()), mv_scalar(1.0)));
    assert!(mv_almost(r.inverse(), r.reverse()));
}

#[test]
fn interpolate_midpoint() {
    let mid = interpolate(
        motor_identity(),
        motor([0.0, 0.0, 1.0], std::f64::consts::PI / 2.0, [0.0, 0.0, 0.0]),
        0.5,
    );
    let got = mid.apply(point(1.0, 0.0, 0.0)).coords();
    assert!((got[0] - (std::f64::consts::PI / 4.0).cos()).abs() < 1e-6);
    assert!((got[1] - (std::f64::consts::PI / 4.0).sin()).abs() < 1e-6);
}

#[test]
fn str_renders_non_zero_components() {
    let s = e1().add(e2().op(e3())).str();
    assert!(s.contains("e1"));
    assert!(s.contains("e23"));
    assert!(mv_zero().str().contains('0'));
}
