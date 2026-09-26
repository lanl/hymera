# Bit-exactness policy

Normative rules for changing the C MHD solver. Every change declares a class
before work starts. The class fixes what must be proved and how.

Scope: `src/mhd/*.c` and `*.h`. The C++/Parthenon side is out of scope except
where the C ABI forces a change.

## The build guarantees this policy depends on

From `build/compile_commands.json`, the solver compiles with:

```
-O2 -g -DNDEBUG
```

and nothing else of numerical significance. Specifically verified absent: `-flto`,
`-ffast-math`, `-Ofast`, `INTERPROCEDURAL_OPTIMIZATION`. Also verified, via `nm`
over all eight `mhd_core` objects: no `.data`, no `.bss`, no file-scope statics,
no static functions anywhere.

Three consequences, and the whole policy rests on them:

1. **Deletion is provable without running anything.** With no LTO there is no
   legal cross-function inlining path between translation units, and with no
   statics an `extern` function's removal cannot alter any surviving function's
   generated code. So for a pure deletion the proof is: compile before and after,
   and confirm every surviving function's instruction stream is unchanged.
2. **Hoisting is exact, not approximate.** No `-ffast-math` means the compiler
   does not reassociate, and on this target doubles are not held at extended
   precision. Recomputing `f(x)` twice and caching it once give bit-identical
   results when the arguments are identical.
3. **All mutable state is reachable from one `User*`.** Nothing hides in globals,
   so a change's blast radius is bounded by what it touches in that struct.

**If any of these flags changes, this policy is void.** A test should assert they
stay off.

## Classes

### `EXACT` — 0 ULP, no exceptions

Qualifying transformations:

- deleting a dead function, variable, or branch
- adding `PetscCall`, `PetscFunctionBeginUser`, `PetscFunctionReturn`
- renaming anything
- extracting code into a helper **preserving operation order and association**
- replacing a magic literal with a named constant of the identical value
- hoisting a value recomputed from identical arguments (see consequence 2)
- comment and formatting changes

Proof, strongest first:

1. **Object-code identity.** Normalize `objdump -d` output (strip addresses and
   relocation offsets, keep mnemonics, operands, symbol names) and diff per
   function. Complete proof; runs in seconds; needs no fixture, no grid, no MPI.
   Applies to deletion and renaming. Does **not** apply to `PetscCall` insertion,
   which legitimately adds branches.
2. **Fingerprint identity at 0 ULP.** For `PetscCall` insertion, helper
   extraction, and predicate substitution. Compare the full run fingerprint; every
   value must match as a string.
3. **Symbol diff.** `nm --defined-only` and `nm -u`, sorted, before and after.
   Mandatory alongside 1 and 2 — it catches deleting more than intended.

One expected complication: wrapping a call whose return value is currently
discarded will abort if that call is *already* failing silently in normal
operation. That is a **bug discovery, not a regression**. Report it and stop; do
not suppress it. Expect this at least once in `mhd.c`'s construction path, which
is entirely unchecked.

### `REASSOC` — 1e-13 relative, justified per site

Changing sum order or association; replacing `a*b + a*c` with `a*(b+c)`;
replacing repeated division with multiplication by a reciprocal.

Proof: a **single residual evaluation**, not a multi-step run. Report the observed
`‖ΔX‖∞/‖X‖∞`. If it is exactly 0, reclassify the change as `EXACT`. If it exceeds
1e-13, stop and report.

Why a single evaluation: see amplification, below.

### `SOLVER` — 1e-8 relative, and the dangerous one

Anything that can shift the Newton or Krylov path.

`TSConvergedReason` and the step count must match **exactly**. Iteration counts
may vary and are reported but not asserted. State values compare at 1e-8
relative.

## Solver-tolerance amplification

`-ts_adapt_type none` fixes the step sequence, so `tⁿ` is reproducible. But `Xⁿ`
is whatever the inexact Newton–Krylov solve converged to, and convergence is only
to `snes_rtol`/`ksp_rtol`. Therefore:

- A 1e-16 perturbation of the residual perturbs the Newton step, which can change
  a KSP iteration count by one, which moves `X¹` by `O(ksp_rtol·‖X‖)` —
  potentially 1e-8. Eight orders of amplification, in one step.
- Over many steps it can grow further, and near a marginal mode in resistive MHD
  it can grow a great deal further.

Practical rules:

1. **Never validate a `REASSOC` change with a multi-step run alone.** The signal
   is buried in solver noise. Use a single residual evaluation, where a 1e-16
   change stays 1e-16.
2. **Measure the noise floor before trusting any multi-step tolerance.** Perturb
   one input by 1 ULP and record the resulting `‖ΔX‖` at the last step. That
   number is the smallest difference a multi-step run can resolve. If it exceeds
   the tolerance you intended to assert, then the multi-step run can only check
   structural properties — convergence reason, step count, no NaN, bounded div-B
   ratio — and the numerical burden falls on the single-evaluation check.
3. **Tightening tolerances is a legitimate test tool.** Regression runs may set
   `-snes_rtol`/`-ksp_rtol` far below production to shrink the floor. This became
   possible only after the command-line precedence fix; previously the hardcoded
   options silently overrode anything passed on the command line.

## Universal assertions

Every run, every class:

- exit code 0, no NaN, no Inf
- no `DIVERGED` in the log
- `max‖div B‖/max‖B‖` within 2x of baseline

The last one is the sharpest detector available. A mimetic scheme satisfies
discrete `div·curl = 0` exactly, so this quantity sits at machine zero
(measured: 1.34388e-13). A broken discretization moves it by orders of magnitude
even when the solution norms look plausible.

## Prohibited

Each of these is a reasonable-sounding idea that would cause disproportionate
damage. Do not attempt them as part of a refactoring pass.

1. **Multi-instance or re-entrancy work.** `PETSC_COMM_WORLD` is hardcoded in
   ~14 places in `mhd.c`, and `PetscInitialize`/`PetscFinalize` are baked into
   `mhd_PetscInit`/`mhd_destroy`. Structurally possible later given the absence
   of globals; no requirement and no coverage now.
2. **Splitting or reordering the `User` struct.** ~90 fields, and
   zero-initialization via `calloc` is load-bearing: it is the only initializer
   for roughly 40 fields nobody assigns. Permitted: add comments, add banner
   comments without reordering, add fields.
3. **Changing the `stag_vec_io` file format** before the dedicated hardening
   step. Its rank-independence is what makes state comparison possible at all,
   and because records are matched purely by loop order, any change invalidates
   every stored baseline.
4. **Removing the unscoped macros in `mfd_config.h`** (`LEFT`, `RIGHT`, `UP`,
   `DOWN`, `BACK`, `FRONT`, `ELEMENT`, ...) as part of another change. Thousands
   of sites across all six `.c` files. It needs one dedicated machine-applied
   whole-repo change with an object-code identity proof, never concurrent with
   anything else.
5. **Re-enabling the `if(0)` inside `ComputeCurrent`**
   (`monitor_functions.c`, near line 2511). It changes reported physics. Owner
   decision.
6. **Removing `target_link_libraries(infrastructure PUBLIC mhd_core)`** from
   `src/CMakeLists.txt`. Load-bearing for static link ordering; without it,
   executables whose driver does not directly reference the MHD symbols fail to
   link.

## Working rules for scoped changes

- **One file, one owner, one change at a time.** Never two concurrent edits in
  the same `.c` file even in disjoint line ranges: line numbers shift and every
  later range goes stale. Parallelize across files, serialize within a file.
- **Anchor on symbol names, never line numbers.** Locate a function by grepping
  for its definition and taking its extent to the matching closing brace.
  Regenerate `docs/mhd/reference/function-index.md` after each merge.
- **Verify after every change:** clean build with no new warnings; the class's
  proof from above; symbol diff; and a diffstat whose size matches expectation.
