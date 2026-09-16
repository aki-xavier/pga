// primitives.rs — the geometry layer for PGA: points, planes, lines and the
// distance / angle / projection toolbox. Port of primitives.v.
//
//   plane:  n + d e0            (unit normal n, signed distance d)
//   line:   bivector direction + moment (join of two points; Euclidean lines
//           square to -|d|^2, ideal lines are null)
//   point:  trivector e123 + x e032 + y e013 + z e021
//
// Native incidence: point on plane  <=> point ^ plane = 0
//                   point on line   <=> line ^ point = 0
//                   plane ^ plane   = line;  plane ^ line = point.

use crate::multivector::{e0, e1, e123, e2, e3, mv_vector, Multivector};

/// point returns the PGA point at (x, y, z).
pub fn point(x: f64, y: f64, z: f64) -> Multivector {
    let mut p = e123();
    p = p.add(e0().op(e3()).op(e2()).scale(x));
    p = p.add(e0().op(e1()).op(e3()).scale(y));
    p = p.add(e0().op(e2()).op(e1()).scale(z));
    p
}

/// plane returns the unit plane n . x + d = 0 (normal normalized internally).
pub fn plane(normal: [f64; 3], d: f64) -> Multivector {
    let nl = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    if nl < 1e-12 {
        panic!("plane normal is zero");
    }
    mv_vector(normal[0] / nl, normal[1] / nl, normal[2] / nl, d)
}

/// line_from_points returns the join of two points: the line through them
/// (the ideal line if the points coincide in direction).
pub fn line_from_points(p1: Multivector, p2: Multivector) -> Multivector {
    p1.join(p2)
}

/// line returns the line through `p` in direction `dir` (direction need not be
/// unit; the vanish point of `dir` is joined with p).
pub fn line(p: Multivector, dir: [f64; 3]) -> Multivector {
    let u = mv_vector(dir[0], dir[1], dir[2], 0.0);
    let vanishing = e0().op(u.lc(e123()));
    p.join(vanishing)
}

impl Multivector {
    /// coords extracts the euclidean coordinates (x, y, z) of a trivector point
    /// (normalising by the e123 weight, so scaled or reflected points work).
    pub fn coords(&self) -> [f64; 3] {
        let w = self.values[7];
        if w.abs() < 1e-9 {
            panic!("multivector has no e123 component; not a finite point");
        }
        [
            -self.values[14] / w,
            self.values[13] / w,
            -self.values[11] / w,
        ]
    }

    /// ideal_direction returns the direction of a line: the direction part of its
    /// vanishing point (line ^ e0), which is the ideal point of the line.
    pub fn ideal_direction(&self) -> [f64; 3] {
        let v = self.meet(e0());
        [-v.values[14], v.values[13], -v.values[11]]
    }

    /// point_on_line returns either a finite point of `line` (intersecting it
    /// with the first axis plane that is not parallel to it) or the ideal point.
    pub fn point_on_line(&self) -> Multivector {
        for normal in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            let q = self.meet(plane(normal, 0.0));
            if !q.is_zero() && q.values[7].abs() > 1e-6 {
                return q;
            }
        }
        self.meet(e0())
    }

    /// reflect mirrors self across the plane with the given normal: X' = p X p
    /// with p the unit plane vector.
    pub fn reflect(&self, normal: [f64; 3]) -> Multivector {
        let p = plane(normal, 0.0);
        p.gp(*self).gp(p)
    }
}

// --- distances ---------------------------------------------------------------

/// point_dist returns the euclidean distance between two points.
pub fn point_dist(a: Multivector, b: Multivector) -> f64 {
    let c1 = a.coords();
    let c2 = b.coords();
    let dx = c1[0] - c2[0];
    let dy = c1[1] - c2[1];
    let dz = c1[2] - c2[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// point_plane_dist returns the signed distance of point p to plane pi
/// (positive on the side of the normal).
pub fn point_plane_dist(p: Multivector, pi: Multivector) -> f64 {
    p.meet(pi).values[15]
}

/// point_line_dist returns the distance from point p to the line l.
pub fn point_line_dist(p: Multivector, l: Multivector) -> f64 {
    let q = l.point_on_line();
    let d = l.ideal_direction();
    let an = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if an < 1e-12 {
        return point_dist(p, q);
    }
    let c = p.coords();
    // vector from q to p
    let qc = q.coords();
    let rx = c[0] - qc[0];
    let ry = c[1] - qc[1];
    let rz = c[2] - qc[2];
    // distance = |r x d| / |d|
    let crx = ry * d[2] - rz * d[1];
    let cry = rz * d[0] - rx * d[2];
    let crz = rx * d[1] - ry * d[0];
    (crx * crx + cry * cry + crz * crz).sqrt() / an
}

/// line_line_dist returns the distance between two lines (0 for intersecting
/// lines).
pub fn line_line_dist(a: Multivector, b: Multivector) -> f64 {
    let da = a.ideal_direction();
    let db = b.ideal_direction();
    let na = (da[0] * da[0] + da[1] * da[1] + da[2] * da[2]).sqrt();
    let nb = (db[0] * db[0] + db[1] * db[1] + db[2] * db[2]).sqrt();
    if na < 1e-12 || nb < 1e-12 {
        return 1e300;
    }
    // common normal direction
    let nx = da[1] * db[2] - da[2] * db[1];
    let ny = da[2] * db[0] - da[0] * db[2];
    let nz = da[0] * db[1] - da[1] * db[0];
    let nn = (nx * nx + ny * ny + nz * nz).sqrt();
    if nn < 1e-12 {
        // parallel lines: distance from a point of b to line a
        return point_line_dist(b.point_on_line(), a);
    }
    let pa = a.point_on_line();
    let pb = b.point_on_line();
    let ca = pa.coords();
    let cb = pb.coords();
    // (b - a) . n / |n|
    let dx = cb[0] - ca[0];
    let dy = cb[1] - ca[1];
    let dz = cb[2] - ca[2];
    (dx * nx + dy * ny + dz * nz).abs() / nn
}

// --- angles ------------------------------------------------------------------

/// plane_angle returns the angle (radians, [0, pi/2]) between two planes.
pub fn plane_angle(a: Multivector, b: Multivector) -> f64 {
    let mut dot = a.ip(b).values[0].abs();
    if dot > 1.0 {
        dot = 1.0;
    }
    dot.acos()
}

/// line_angle returns the angle (radians, [0, pi/2]) between two lines.
pub fn line_angle(a: Multivector, b: Multivector) -> f64 {
    let da = a.ideal_direction();
    let db = b.ideal_direction();
    let na = (da[0] * da[0] + da[1] * da[1] + da[2] * da[2]).sqrt();
    let nb = (db[0] * db[0] + db[1] * db[1] + db[2] * db[2]).sqrt();
    if na < 1e-12 || nb < 1e-12 {
        return 0.0;
    }
    let dot = ((da[0] * db[0] + da[1] * db[1] + da[2] * db[2]) / (na * nb)).clamp(-1.0, 1.0);
    dot.abs().acos()
}

// --- projections -------------------------------------------------------------

/// project_point_onto_plane projects p onto the plane pi along its normal.
pub fn project_point_onto_plane(p: Multivector, pi: Multivector) -> Multivector {
    let pn = pi.vector_part();
    let c = p.coords();
    let dist = pn[0] * c[0] + pn[1] * c[1] + pn[2] * c[2] + pi.values[8];
    point(
        c[0] - dist * pn[0],
        c[1] - dist * pn[1],
        c[2] - dist * pn[2],
    )
}

/// project_point_onto_line projects p onto the line l.
pub fn project_point_onto_line(p: Multivector, l: Multivector) -> Multivector {
    let q = l.point_on_line();
    let d = l.ideal_direction();
    let an = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if an < 1e-12 {
        return q;
    }
    let c = p.coords();
    let qc = q.coords();
    let lam = ((c[0] - qc[0]) * d[0] + (c[1] - qc[1]) * d[1] + (c[2] - qc[2]) * d[2]) / (an * an);
    point(qc[0] + lam * d[0], qc[1] + lam * d[1], qc[2] + lam * d[2])
}
