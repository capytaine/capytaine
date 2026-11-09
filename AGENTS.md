# Instructions for agents

Help the user follow the AI policy of this repository described in `CONTRIBUTING.md`.

## Dev environment

See also `docs/developer_manual/installation.rst`, `docs/developer_manual/testing.rst`, `Justfile`.

```bash
sudo apt install python3 python3-pip python3-venv python-is-python3 python3-dev git gfortran just
python -m venv /tmp/capy_dev_venv    # NOT inside the source tree (meson issue #396)
source /tmp/capy_dev_venv/bin/activate
pip install uv                       # needed by the just recipe below
just editable_install                # see justfile and pyproject.toml for details
```

Python edits are picked up immediately; the Fortran core is rebuilt on import after Fortran changes.

## Testing

The full test suite can be run in the current environment with `just _test`.

To run only the pytest suite (about 1 min), e.g. to select out a single test:

```bash
OMP_NUM_THREADS=1 MPLBACKEND=pdf CAPYTAINE_PROGRESS_BAR=False CAPYTAINE_CACHE_DIR=<tmpdir> python -m pytest pytest/
```

`OMP_NUM_THREADS=1` halves the run time (Fortran multithreading is overhead on small test meshes).
Some tests are skipped when optional dependencies are missing; a few are `xfail`.
If necessary, optional dependencies (`vtk`, `gmsh`, ...) can be skipped.
In particular, if `gmsh` fails at import with `OSError: libGLU.so.1`, either `apt install libglu1-mesa` or uninstall `gmsh` altogether.

## CI reference environments

CI runs `just _test` in isolated environments pinned by the lockfiles in `pytest/envs/`:

```bash
# Assuming `uv` is installed, for instance by sourcing /tmp/capy_dev_venv
just test_in_py38_reference_env > /tmp/ref38.log 2>&1; echo "exit $?" >> /tmp/ref38.log
just test_in_py314_reference_env > /tmp/ref314.log 2>&1; echo "exit $?" >> /tmp/ref314.log
```

- Need network; build the **current working tree** (uncommitted changes included).
- Take several minutes: run both in parallel in the background (the tests are single-threaded,
  `OMP_NUM_THREADS=1`) and wait for both to end.
- After adding a dependency to `pyproject.toml`, regenerate the lockfiles with
  `just create_test_env_file 3.8 <date>` and `just create_test_env_file 3.14 <date>`
  (dates as in `pytest/envs/` file names; not older than the required dependency releases).

## Style

### Tests

- Avoid pytest fixtures: use plain functions (possibly `lru_cache`d) or context managers,
  and normal Python scoping. Pytest builtin fixtures are acceptable.
- Shared helpers are in `pytest/capytaine_test_helpers.py`, used as
  `import capytaine_test_helpers as helpers` (module import avoids name clashes, e.g. `solver`).
  Reuse them:
  - `helpers.solver(method=..., **gf_kwargs)` instead of `cpt.BEMSolver()` (reuses a cached
    `helpers.green_function(**gf_kwargs)`; building `Delhommeau()` reloads its tabulation from disk).
  - `helpers.small_sphere_mesh()` / `helpers.small_sphere_body(dofs, lid=..., name=...)` when the geometry does not matter.
  - `helpers.broken_bem_solver()`, `helpers.single_panel()`, `helpers.symmetric_meshes_of_single_panel()`.
- Only cache immutable objects (meshes, Green functions). Bodies and solvers are mutable: build them anew.

### Changelog

List changes in `docs/changelog.rst`, under one of `Major changes`, `Minor changes`,
`Bug fixes`, `Internals`. Ask the user if unsure.
Before adding an item, read the unreleased section: if the change extends or renames something
that is itself unreleased, update that item instead of adding a new one (no new section either).

## Commit

If you make a commit, follow conventional commits and add a trailer:
`Assisted-by: <harness>:<model>`, where `<harness>` is the current agent
harness, and `<model>` is the AI model.
Do not add AI co-author trailers, only humans are co-authors.
Avoid emojis in commit messages.

The repository has a pre-commit hook (`.pre-commit-config.yaml`; whitespace, end of file, Python
syntax, `codespell`, no commit on `master`). It is shared by all worktrees, and fails with
"`pre-commit` not found" unless `pre-commit` is on `PATH`:

```bash
source /tmp/capy_dev_venv/bin/activate   # or any venv: pre-commit does not need capytaine
pip install pre-commit
pre-commit install                       # only if .git/hooks/pre-commit is missing
```

- The first run fetches the hook repositories (network needed, then cached in `~/.cache/pre-commit`).
- Do not bypass the hooks with `--no-verify`.
- When a hook fixes files (trailing whitespace, final newline), the commit is aborted: review the
  changes, `git add` them and commit again.
