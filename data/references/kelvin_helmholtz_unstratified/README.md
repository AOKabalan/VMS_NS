# Unstratified Kelvin–Helmholtz reference data

Retrieved 2026-11-03 from
[KHdata/KelvinHelmholtz](https://gitlab.gwdg.de/KHdata/KelvinHelmholtz),
commit `42eba2cf80ee5e9d59903274fd72c7f2f5d93cbe` (2018-07-27).
The [benchmark specification](https://www.ngsolve.org/kh-benchmark/) accompanies
Schroeder, John, Lederer, Lehrenfeld, Lube and Schöberl,
[On reference solutions and the sensitivity of the 2D Kelvin–Helmholtz instability problem](https://arxiv.org/abs/1803.06893).

This is a constant-density, incompressible Navier–Stokes benchmark without
thermal buoyancy or moisture. It supplies a hydrodynamic KH reference; it does
not validate the temperature/buoyancy coupling of our stratified case. Its
initial profiles, domain and scaling also differ from that application.

## Local files

- `upstream-42eba2cf.tar.gz`: complete, unchanged upstream snapshot, 24,438,580 bytes.
- `upstream/`: unpacked copy, preserving all upstream paths and file bytes.
- `manifest.json`: source/revision, archive SHA-256, per-file sizes, SHA-256 and
  Git blob identifiers, plus numeric-table and integrity checks.

The snapshot contains 495 files and 92,312,986 unpacked bytes: 61 time-series
tables, 401 spectrum tables, five VTK snapshots, 26 READMEs, and a plotting
example in TEX/PDF. The archive is versioned; the duplicate unpacked tree is
locally available and ignored by Git. To recreate it from the archive:

```bash
mkdir -p data/references/kelvin_helmholtz_unstratified/upstream
tar -xzf data/references/kelvin_helmholtz_unstratified/upstream-42eba2cf.tar.gz \
  --strip-components=1 -C data/references/kelvin_helmholtz_unstratified/upstream
```

Archive SHA-256:
`ddc00d92f51c1e4d80cc815cf01e1e9f9e7441d17bdb28bad575d9ccdb01666b`.
No separate LICENSE file is present in this revision; upstream attribution and
documentation are preserved in the snapshot.

## Data organization

- `upstream/reference_results/`: proposed reference runs at Re = 100, 1000 and
  10000, with multiple mesh resolutions. The finest Re=10000 directory also
  contains polynomial-order comparisons, spectra and VTK fields.
- `upstream/perturbation_results/`: sensitivity experiments involving numerical
  integration, linear solves, compiler arithmetic and mesh choice. These are
  distinct from the proposed reference solutions.
- `upstream/how_to_plot/`: upstream plotting example and copies of its input tables.

One upstream inconsistency is retained: the perturbation README says Re=100000,
while its filenames and the benchmark reference describe Re=10000. No upstream
file has been corrected or reclassified locally.

## Table conventions

Despite the `.csv` extension, the tables are **tab-separated**, with a `#`
commented header. All 61 have 401 rows. They can be read with
`numpy.loadtxt(path, comments="#")`, using its default whitespace separator.

`Cons2DKH-*.csv` columns are physical time, scaled time, kinetic energy,
enstrophy, palinstrophy, and the L2 divergence norm. `VortThick-*.csv` columns
are scaled time and normalized vorticity thickness.

The upstream `how_to_plot/plot_example.tex` confirms:

- scaled time is `physical_time * u_infinity / delta_0 = 28 * physical_time`;
- kinetic energy is `0.5 * integral |u|^2`;
- enstrophy is `0.5 * integral |curl u|^2`;
- palinstrophy is `0.5 * integral |grad_h curl u|^2`;
- the stored vorticity thickness is `delta(t) / delta_0`.

Here `u_infinity=1`, `delta_0=1/28`, and `Re=1/(28*nu)`. Use the downloaded
table/plot conventions when comparing quantities: the benchmark web page's
displayed kinetic-energy expression omits the factor `0.5`. An independent
integration of its specified initial condition gives `0.48221164950178574`,
matching the finest Re=1000 table's initial value `0.48221164950178835`.

## Download verification

The gzip integrity check passed and every extracted file matches its upstream
Git blob identifier. All 462 numeric tables parse with finite values; time-series
times are strictly increasing, and the two time columns obey the stated factor
28. All five VTK headers were checked. These checks verify the downloaded data;
no solver comparison has been run in this download task.
