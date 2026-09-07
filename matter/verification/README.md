# Verification for matter-sector claims

Empty by design: there is no matter-sector claim yet.

A numerical claim merges only if it ships with the code that produces it and a
check here that fails when the claim fails. Follow the shape of
`../../v2.9/verification/verify_v291.py`: a `check(name, condition, detail)`
helper, one block per claim, a `RESULT: n/m PASS` line and a non-zero exit on any
failure.

Two rules specific to this sector, from `../CONSTRAINTS.md`:

- a check that a computed **ratio** matches a measured one is a real check;
  a check that a mass matches after an index was chosen to make it match is not —
  it re-derives its own input (§4);
- state the maturity level in the check name. `Reparametrized` and `Derived` must
  not be reported the same way.
