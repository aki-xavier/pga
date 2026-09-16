// primitives.rs — geometry-layer tests, ported from primitives_test.v:
// coordinates, distances, angles, projections and the PGA join/meet
// constructions on real primitive data.

use pga::*;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

#[test]
fn point_plane_incidence() {
    assert!(close(point(1.0, 2.0, 3.0).coords()[0], 1.0, 1e-9));
    assert!(close(point(1.0, 2.0, 3.0).coords()[1], 2.0, 1e-9));
    assert!(close(point(1.0, 2.0, 3.0).coords()[2], 3.0, 1e-9));
    // signed distances (plane z + 2 = 0, i.e. z = -2)
    let pi = plane([0.0, 0.0, 1.0], 2.0);
    assert!(close(point_plane_dist(point(0.0, 0.0, 5.0), pi), 7.0, 1e-9));
    assert!(close(point_plane_dist(point(0.0, 0.0, 0.0), pi), 2.0, 1e-9));
    // a point on the plane gives 0
    assert!(close(
        point_plane_dist(point(1.0, -1.0, -2.0), pi),
        0.0,
        1e-9
    ));
}

#[test]
fn distances() {
    // 3-4-5 point distance
    assert!(close(
        point_dist(point(0.0, 0.0, 0.0), point(3.0, 4.0, 0.0)),
        5.0,
        1e-9
    ));
    // point-line: (1,0,3) to the z-axis
    let lz = line_from_points(point(0.0, 0.0, 0.0), point(0.0, 0.0, 1.0));
    assert!(close(point_line_dist(point(1.0, 0.0, 3.0), lz), 1.0, 1e-9));
    assert!(close(
        point_line_dist(point(2.0, 2.0, 10.0), lz),
        8.0f64.sqrt(),
        1e-9
    ));
    // line-line: skew lines
    let l1 = line_from_points(point(1.0, 2.0, 0.0), point(1.0, 2.0, 1.0));
    let l2 = line_from_points(point(0.0, 0.0, 0.0), point(0.0, 0.0, 1.0));
    assert!(close(line_line_dist(l1, l2), 5.0f64.sqrt(), 1e-9));
    // intersecting lines give 0
    let ix = line_from_points(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0));
    let iy = line_from_points(point(0.0, 0.0, 0.0), point(0.0, 1.0, 0.0));
    assert!(close(line_line_dist(ix, iy), 0.0, 1e-9));
}

#[test]
fn angles() {
    // planes: 45 deg between x+y=0 and x=0, 90 deg between z=0 and x=0
    assert!(close(
        plane_angle(plane([1.0, 1.0, 0.0], 0.0), plane([1.0, 0.0, 0.0], 0.0)),
        std::f64::consts::PI / 4.0,
        1e-9
    ));
    assert!(close(
        plane_angle(plane([0.0, 0.0, 1.0], 0.0), plane([1.0, 0.0, 0.0], 0.0)),
        std::f64::consts::PI / 2.0,
        1e-9
    ));
    // lines: the x-axis and the z-axis are 90 deg apart; parallel lines give 0
    let lx = line_from_points(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0));
    let lz = line_from_points(point(0.0, 0.0, 0.0), point(0.0, 0.0, 1.0));
    assert!(close(line_angle(lx, lz), std::f64::consts::PI / 2.0, 1e-9));
    let lx2 = line_from_points(point(0.0, 0.0, 0.0), point(2.0, 0.0, 0.0));
    assert!(close(line_angle(lx, lx2), 0.0, 1e-9));
}

#[test]
fn join_constructions() {
    // line through two points: passes them and its direction is their offset
    let l = line_from_points(point(0.0, 0.0, 0.0), point(1.0, 0.0, 0.0));
    let d = l.ideal_direction();
    assert!(d[0].abs() > 0.9); // direction along x (scale-normalised below)
                               // collinear points join to the same line
    assert!(close(
        point_line_dist(
            point(1.0, 2.0, 0.0),
            line_from_points(point(0.0, 2.0, 0.0), point(3.0, 2.0, 0.0))
        ),
        0.0,
        1e-9
    ));
    // plane from three points: (0,0,0),(1,0,0),(0,1,0) span the plane z=0
    let pl = point(0.0, 0.0, 0.0)
        .join(point(1.0, 0.0, 0.0))
        .join(point(0.0, 1.0, 0.0));
    assert!(close(point_plane_dist(point(0.3, 0.7, 0.0), pl), 0.0, 1e-9));
    // the join's orientation is arbitrary: test the unsigned distance
    assert!(close(
        point_plane_dist(point(0.0, 0.0, 5.0), pl).abs(),
        5.0,
        1e-9
    ));
}

#[test]
fn projections() {
    let pi = plane([0.0, 0.0, 1.0], 0.0);
    let q = project_point_onto_plane(point(1.0, 2.0, 3.0), pi);
    assert!(close(point_plane_dist(q, pi), 0.0, 1e-9));
    assert!(close(point_dist(q, point(1.0, 2.0, 0.0)), 0.0, 1e-9));
    // onto a line: (1,0,3) projects to (0,0,3) on the z-axis
    let lz = line_from_points(point(0.0, 0.0, 0.0), point(0.0, 0.0, 1.0));
    let q2 = project_point_onto_line(point(1.0, 0.0, 3.0), lz);
    assert!(close(point_line_dist(q2, lz), 0.0, 1e-9));
    assert!(close(point_dist(q2, point(0.0, 0.0, 3.0)), 0.0, 1e-9));
}

#[test]
fn reflect() {
    let rf = point(1.0, 2.0, 3.0).reflect([0.0, 0.0, 1.0]);
    let c = rf.coords();
    assert!(close(c[0], 1.0, 1e-9) && close(c[1], 2.0, 1e-9) && close(c[2], -3.0, 1e-9));
    // reflecting twice returns the point
    let back = rf.reflect([0.0, 0.0, 1.0]).coords();
    assert!(close(back[0], 1.0, 1e-9) && close(back[2], 3.0, 1e-9));
}

#[test]
fn to_matrix() {
    let m = motor([0.0, 0.0, 1.0], std::f64::consts::PI / 2.0, [1.0, 2.0, 3.0]);
    let mtx = m.to_matrix();
    assert!(close(mtx[3], 1.0, 1e-9) && close(mtx[7], 2.0, 1e-9) && close(mtx[11], 3.0, 1e-9));
    // row-major [R|t]: R e1 = first column (mtx[0], mtx[4], mtx[8]) = (0,1,0)
    assert!(close(mtx[0], 0.0, 1e-9) && close(mtx[4], 1.0, 1e-9) && close(mtx[8], 0.0, 1e-9));
    // matrix action matches the motor action on a point
    let p = m.apply(point(1.0, 0.0, 0.0)).coords();
    assert!(close(
        mtx[0] * 1.0 + mtx[1] * 0.0 + mtx[2] * 0.0 + mtx[3],
        p[0],
        1e-9
    ));
    assert!(close(
        mtx[4] * 1.0 + mtx[5] * 0.0 + mtx[6] * 0.0 + mtx[7],
        p[1],
        1e-9
    ));
    assert!(close(
        mtx[8] * 1.0 + mtx[9] * 0.0 + mtx[10] * 0.0 + mtx[11],
        p[2],
        1e-9
    ));
}
