# CuttingPlane.jl — V-Polyhedral Disjunctive Cuts in Julia

A compact, readable implementation of **V-polyhedral disjunctive cuts (VPCs)**
for mixed-integer linear programs, with an optional **monoidal strengthening**
step. Built on [JuMP](https://jump.dev) and the open-source
[HiGHS](https://highs.dev) solver.

Given a MIP, the code solves its LP relaxation, splits on a fractional integer
variable, describes each side of the split by a vertex and the rays of the LP
cone at that vertex, and solves a small *point-ray LP* (PRLP) to obtain a cut
`α'x ≥ β` that removes the LP optimum while keeping every integer-feasible point.
It then reports how much of the integrality gap the cut closes.

---

## Contents

- [Quick start](#quick-start)
- [Method](#method)
- [Repository layout](#repository-layout)
- [Usage](#usage)
- [Example output](#example-output)
- [API reference](#api-reference)
- [Assumptions and limitations](#assumptions-and-limitations)
- [Testing](#testing)
- [References](#references)

---

## Quick start

Requires **Julia ≥ 1.9** (developed with 1.11).

```bash
git clone <this-repo-url> Cutting_Plane_Julia
cd Cutting_Plane_Julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # install dependencies once
julia --project=. run.jl                               # run the built-in example
```

---

## Method

The full pipeline lives in [`run.jl`](run.jl); each step maps to one function.

| # | Step | Function(s) |
|---|------|-------------|
| 1 | Solve the MIP (for reference) and its LP relaxation; get the LP optimum `ω`. | `optimize_model!`, `relax_model` |
| 2 | Pick an integer variable `x_k` with fractional value `x̄_k`. | `fractional_integer_variables` |
| 3 | Rewrite the relaxation in standard form `max c'x, Ax ≤ b, x ≥ 0`. | `create_standard_form_model` |
| 4 | Build the split disjunction `x_k ≤ ⌊x̄_k⌋ ∨ x_k ≥ ⌈x̄_k⌉` (two subproblems). | `create_disjunctive_subproblems` |
| 5 | For each side: add slacks (normal form), solve, build the optimal simplex tableau, read off the optimal vertex `p_t` and the cone rays `R_t`. | `process_subproblem` → `create_normal_form_model`, `generate_optimal_simplex_tableau`, `compute_rays` |
| 6 | Solve the PRLP for the cut `α'x ≥ β`. | `solve_prlp` |
| 7 | Add the cut to the relaxation, re-solve, report the percentage of the gap closed. | `add_cut`, `gap_closed` |
| 8 | *(Optional)* Strengthen the cut coefficients with monoidal strengthening. | `prepare_monoidal_data`, `monoidal_strengthening` |

### Point-ray LP

Every point in the union of the two sides lies in `conv{p_t} + cone{R_t}`, so
any `(α, β)` satisfying

```
α'p_t ≥ β          for each side's vertex p_t
α'r   ≥ 0          for each ray r of each side
```

is a valid cut. Among those, the PRLP picks the one most violated at `ω`.
Two normalizations are available (`normalization` keyword of `solve_prlp`):

| Option | Formulation | Notes |
|--------|-------------|-------|
| `:beta` *(default)* | fix `β = 1`, minimize `α'ω` | Original formulation. Unbounded if a valid cut with `β = 0` separates `ω`; infeasible if every valid cut has `β ≤ 0`. |
| `:l1` | `β` free, `‖α‖₁ ≤ 1`, minimize `α'ω − β` | Always bounded; the most violated cut under ℓ1 scaling. |
| `:auto` | try `:beta`, fall back to `:l1` | Used by `run.jl`. |

### Gap closed

```
gap closed (%) = 100 · (z_cut − z_LP) / (z_IP − z_LP)
```

where `z_LP` is the LP bound, `z_IP` the MIP optimum and `z_cut` the LP bound
after adding the cut.

---

## Repository layout

```
Cutting_Plane_Julia/
├── Project.toml          # package metadata and dependencies
├── Manifest.toml         # exact dependency versions (reproducibility)
├── run.jl                # end-to-end pipeline script (start here)
├── examples/
│   └── toy_model.jl      # 2-variable example MIP
├── src/
│   ├── CuttingPlane.jl   # module definition and exports
│   ├── utils.jl          # read / relax / solve / report helpers
│   ├── forms.jl          # standard-form and normal-form reformulations
│   ├── tableau.jl        # simplex tableaux, basis recovery, cone rays
│   ├── vpc.jl            # disjunctions, PRLP, cut generation
│   └── monoidal.jl       # monoidal strengthening
└── test/
    └── runtests.jl       # validity checks on the toy example
```

---

## Usage

### Command line

```bash
julia --project=. run.jl                       # built-in toy example
julia --project=. run.jl path/to/model.mps     # split on the first fractional integer variable
julia --project=. run.jl path/to/model.mps 3   # split on variable 3 (1-based, in file order)
```

Any format JuMP can read works (`.mps`, `.mps.gz`, `.lp`, `.mof.json`, …).

### As a library

```julia
using CuttingPlane, JuMP

model = read_model("model.mps")
relaxed = relax_model(model);  optimize_model!(relaxed)
ω, _ = get_variable_values(relaxed)

integer_vars = lp_matrix_data(model).integers
k = first(fractional_integer_variables(relaxed, integer_vars))

std = create_standard_form_model(relaxed)
sides, D0 = create_disjunctive_subproblems(std, k, integer_vars)
results = [process_subproblem(s) for s in sides]

_, α, β = solve_prlp(ω, [r.point for r in results], [r.rays for r in results];
                     normalization = :auto)
cut_model = add_cut(relaxed, α, β)
optimize_model!(cut_model)
```

Each call to `process_subproblem` also returns the initial and optimal simplex
tableaux as `DataFrame`s (`first_tableau`, `optimal_tableau`), which is handy
for inspecting or teaching the method.

### Choosing a solver

Every solving function takes an `optimizer` keyword (default `HiGHS.Optimizer`),
so any JuMP-compatible LP solver that reports basis status can be used:

```julia
using Gurobi
optimize_model!(model; optimizer = Gurobi.Optimizer)
```

---

## Example output

`julia --project=. run.jl` on the toy problem

```
min  x₁ + 2x₂
s.t. -2x₁ - 2x₂ ≤ -0.5
      2x₁ -  x₂ ≤  1.5
      x₁ ≤ 1,  x₂ ≤ 1,  x ∈ ℤ²₊
```

prints (abridged):

```
╔═══════════════╗
║ LP relaxation ║
╚═══════════════╝
  x[1] = 0.25
  x[2] = 0.0
LP bound: 0.25    MIP optimum: 2.0
Splitting on variable 1 (x[1] = 0.25)
PRLP with β = 1 has no optimal solution; using the ℓ1-normalised PRLP instead.

╔══════════════════╗
║ V-polyhedral cut ║
╚══════════════════╝
Added cut: cut : -0.2 x[1] + 0.8 x[2] ≥ 0.2
Bound with cut: 0.5    gap closed: 14.29%
```

---

## API reference

| Function | Description |
|----------|-------------|
| `read_model(path)` | Read a model from file. |
| `relax_model(model)` | Copy of the model with integrality removed. |
| `optimize_model!(model; optimizer, verbose)` | Attach a solver and solve. |
| `get_variable_values(model)` | `(values, names)` of all variables; errors if not optimal. |
| `retrieve_objective_values(model)` | `(objective_value, termination_status)`. |
| `gap_closed(z_cut, z_int, z_lp)` | Percentage of integrality gap closed. |
| `create_standard_form_model(model)` | `max c'x, Ax ≤ b, x ≥ 0` reformulation. |
| `create_normal_form_model(model)` | Adds one slack per row: `[A I][x; s] = b`. |
| `generate_first_simplex_tableau(m)` | Initial tableau (`DataFrame`). |
| `generate_optimal_simplex_tableau(m)` | Optimal tableau + basis info (named tuple). |
| `compute_rays(tableau, structural, nonbasic)` | Cone rays at the optimal vertex (one column per ray). |
| `find_fractional_elements(values)` | Fractional entries and their indices. |
| `fractional_integer_variables(relaxed, ints)` | Integer variables that are fractional at the LP optimum. |
| `create_disjunctive_subproblems(model, k, ints)` | The two sides of the split on `x_k`, plus their RHS `D0`. |
| `process_subproblem(side)` | Tableaux, rays, vertex and normal-form model for one side. |
| `solve_prlp(ω, points, rays; normalization)` | Point-ray LP → `(model, α, β)`. |
| `add_cut(model, α, β)` | Copy of the model with `α'x ≥ β` added. |
| `prepare_monoidal_data(sides, results, α)` | Active-constraint multipliers per side. |
| `monoidal_strengthening(sides, k, D0, U, α)` | Strengthened coefficients `α̃`. |

All functions have docstrings; in the REPL type `?solve_prlp`, etc.

---

## Assumptions and limitations

- **Non-negative variables.** The standard-form step assumes `x ≥ 0`; finite upper bounds and positive lower bounds are turned into rows.
- **One split, one cut.** A single cut is generated per run, from a single variable split. There is no multi-round cutting loop yet.
- **Dense linear algebra.** Tableaux are built with dense matrices and an explicit basis inverse, which is fine for teaching and small or medium instances but not for large MIPLIB models.
- **Basis from the solver.** The tableau is rebuilt from the basis status the solver reports, so the solver must support `VariableBasisStatus` and `ConstraintBasisStatus` (HiGHS, Gurobi, CPLEX, Clp, … do).
- **Infeasible split side.** If one side of the split is infeasible, the split row itself is the strongest cut, and `run.jl` stops with a message.
- **Monoidal strengthening** is experimental. It uses the right-hand sides `-D0` of the disjunctive rows and the multipliers of the active constraints, and it may leave the cut unchanged.

---

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests generate a cut on the toy model and check that it:

- separates the LP optimum;
- is satisfied by every integer-feasible point (checked by enumeration);
- yields a bound between the LP bound and the MIP optimum.

---

## References

- E. Balas, *Disjunctive Programming*. Springer, 2018.
- E. Balas and R. G. Jeroslow, “Strengthening cuts for mixed integer programs,” *European Journal of Operational Research* 4 (1980) 224–234.
- E. Balas and A. M. Kazachkov, “V-polyhedral disjunctive cuts” (preprint on arXiv).
- M. Lubin et al., “JuMP 1.0: Recent improvements to a modeling language for mathematical optimization,” *Mathematical Programming Computation*, 2023.
