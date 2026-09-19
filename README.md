# pga — Projective Geometric Algebra Cl(3,0,1)

MIT-licensed (see `LICENSE`).

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

The derivation (basis bitmask indexing, the XOR/sign product, the Hodge
complement table, the screw bivector `exp`/`log` closed forms) is in the module
docs of `src/lib.rs` and `src/multivector.rs`.

## Tests

The test suite:

- `tests/multivector.rs` cross-checks the whole 16×16 basis-blade product table
  against an independent reference that folds blades back into generators and
  bubble-sorts them — so a wrong sign or slot is caught against a different
  derivation, not against the table itself;
- `tests/primitives.rs` covers coordinates, distances, angles, projections, the
  join/meet constructions and `to_matrix` on real primitive data.
