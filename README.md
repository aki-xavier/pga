# pga — Projective Geometric Algebra Cl(3,0,1)

Rust port of the `pga` V module, vendored into this repository so `simu` depends
on nothing outside it. The original is MIT-licensed (see `LICENSE`, carried over
unchanged); the algebra, the basis layout, the dual signs and the motion
conventions are its author's, and its test suite is ported here as the oracle.

P(R*_{3,0,1}) = Cl(3,0,1) — the plane-based (dual) PGA: the algebra of Euclidean
3D point/line/plane geometry, with a degenerate fourth generator `e0`
(`e0^2 = 0`) that makes the translation part of a rigid motion multiplicative.

```text
plane = grade-1 vector  n + d e0       (n unit normal, distance d)
line  = grade-2 blade                   (squares negative for Euclidean lines)
point = grade-3 blade                   (null)
meet  = outer product, join = regressive product
motor = even-grade versor, M X M~ applies it
```

The full derivation (basis bitmask indexing, the XOR/sign product, the Hodge
complement table, the screw bivector `exp`/`log` closed forms) is in the module
docs of `src/lib.rs` and `src/multivector.rs`, as it is in the original.

## Port notes

- Identifiers keep the V spelling (`mv_scalar`, `num_components`, …) so the crate
  diffs against the original line by line; `src/lib.rs` says why and carries the
  `allow` for it.
- `popcount`, `gp_blade` and `generator_metric` are public here because the
  ported tests use them from outside the crate (V keeps them module-private).
- V's approximate-equality method `eq` is `approx_eq`: `eq` is the name of a std
  trait method, and this crate derives exact `PartialEq` alongside it.
- The error cases panic, as in V: the inverse of a null element, `exp` of
  anything but scalar+bivector, `log` of the −1 motor, `coords` of a
  non-point, a zero plane normal, a zero rotor axis.
- `blade_grade` is kept although the original never calls it (unused there too).

## Tests

The original suite, ported with its tolerances:

- `tests/multivector.rs` cross-checks the whole 16×16 basis-blade product table
  against an independent reference that folds blades back into generators and
  bubble-sorts them — so a wrong sign or slot is caught against a different
  derivation, not against the table itself;
- `tests/primitives.rs` covers coordinates, distances, angles, projections, the
  join/meet constructions and `to_matrix` on real primitive data.
