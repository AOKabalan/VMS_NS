# Kelvin–Helmholtz shear-layer conventions (B7.1–B7.2)

This is an application/integration case, not an exact solution. B7.1 checks
initialization and a few coupled steps. Reference growth rates, resolution,
nonlinear billows and clouds remain separate studies.

Coordinates are `(x,z)`, with `z` increasing upwards. Length, velocity and thermal
scales are `delta = U0 = DeltaT = 1`. The code uses
`nu = 1/Re`, `kappa_T = 1/(Re Pr)`, `kappa_q = 1/(Re Sc)`,
`expansion_coefficient = Ri`, `gravity = (0,-1)`, `T_ref = 0`, so that
`b(theta) = (0,Ri theta)`. Domain length is `n_waves * 2*pi/wave_number`;
height defaults to 20. The value 0.4446 is an **unstratified inviscid** reference,
not a computed most unstable wavenumber of this bounded stratified problem.
The literature links and qualifications are in
[formulation §8](../../../docs/boussinesq_formulation.md#8-stratified-kelvinhelmholtz-application-with-a-moisture-proxy).

## Actual wall-adjusted profiles

Let `H = height/2`, `d = wall_transition_width`,
`c = H - wall_buffer - d`, `f = H - wall_buffer`.
Defaults give `c=7`, `f=9`, `d=2`. In `|z| <= c`, both horizontal base velocity
and temperature are exactly `tanh(z)`. In the positive transition set
`t=(z-c)/d`, `r=1-t`, `v=tanh(c)`, `a=sech(c)^2`,
`b=a*(1-2*d*v/5)`. Prescribe the gradient

```
G(t) = a*(1-t)^5 + 5*b*t*(1-t)^4
I(t) = a*(1-r^6)/6 + b*(1/6-r^5+5*r^6/6)
J(t) = a/6*(t-(1-r^7)/7)
       + b*(t/6-(1-r^6)/6+5*(1-r^7)/42)
```

The profile is `v+d*I(t)` in the transition and the constant
`v+d*(a+b)/6` outside it. Extend oddly to negative `z`.
Here `I'=G`, `J'=I`. The restriction `0 < d <= 2.5` ensures `b >= 0`;
the gradient stays nonnegative. The value, first and second derivatives match
the core and the constant plateau: at the core `G=a`, `G'/d=-2*v*a`,
and at the plateau `G=G'=0`. Thus no reversed temperature gradient is introduced
by the wall adjustment. The core Richardson formula `Ri*cosh(z)^2` applies only
in the core. No Richardson quotient is assigned to the zero-gradient plateaus.

The compact envelope is 1 in the core, `r^3*(1+3*t+6*t^2)` in the transition,
and 0 in the plateau, extended evenly. It is C². Multiply `exp(-z²)` by this
envelope in the streamfunction, and multiply the initial moisture by the same
envelope. Velocity perturbations and their tangential normal derivative vanish
near each wall; temperature/moisture have zero normal derivatives there.
The smoke mesh has a complete cell row inside each plateau, preserving those
conditions after P1/P2 interpolation. Arbitrarily coarse user meshes may not;
the initialization report measures scalar wall fluxes.

## Hydrostatic pressure and perturbation

The physical pressure initializer is `Ri*P(z)`. In the core
`P=log(cosh(z))`, evaluated with a stable log/exponential expression. In the
positive transition use
`P=log(cosh(c))+d*v*t+d²*J(t)`; continue linearly on the plateau with slope
`v+d*(a+b)/6`. Extend evenly. Thus `P'=theta_0` analytically.
Interpolate the initializer in the pressure space, then remove its mean.
Its gradient generally differs from `Ri*theta_h` after interpolation; record
the L² defect and initial divergence rather than asserting discrete hydrostatic
equilibrium. Initial pressure is a solve guess; the first coupled step determines
pressure through the discrete equations. Velocity/scalar subscale histories
and projection fields start at zero. No force or source holds the profiles steady.

The streamfunction is

```
psi = epsilon * exp(-z²) * envelope(z)
      * (cos(k*x) + epsilon_noise * noise(x))
noise = sum(cos(2*pi*m*x/Lx + phase_m)) / sqrt(number_of_modes)
```

Modes are distinct positive integers (defaults 1,2,3). Phases come from NumPy
`Generator(PCG64(random_seed))`, uniform on `[0,2*pi)`, generated independently
with the same seed on every MPI rank. Metadata stores the actual phases.
For these distinct integer modes the noise has RMS `1/sqrt(2)`, equal to the
base cosine; `epsilon_noise` scales their relative RMS **before addition**.
A noise mode may coincide with the base cosine, so the combined perturbation
RMS need not be the root-sum-square of the two separate RMS values.
Interpolate `(base+psi_z,-psi_x)` and apply the normal velocity constraints.
The continuum seed is solenoidal; its FE interpolant generally is not.

## Moisture and outputs

Transport `q0 = rh_peak*exp(-(z/moisture_width)^2 + alpha*theta0 - gamma*z)`
times the compact envelope, with independent diffusivity `1/(Re Sc)`.
It does not affect momentum or temperature. It is a visualization proxy,
not a condensation model. B7.1 writes `u,p,temperature,q` and named scalar
contents; it does not compute condensate, cloud fraction or growth rates.

`thermal_metadata.json` records physical scales, parameters, mesh counts,
minimum edge, cutoff and noise details. `initial_state.json` records physical
pressure mean, pressure–buoyancy mismatch, vertical energy, divergence,
normal flow and scalar diffusive wall flux norms. `advective_cfl` reports the
MPI-global maximum of `k*dt*|u_h|/MinCellEdgeLength` sampled at the common
quadrature points (`k` is velocity degree). It is a sampled endpoint maximum,
not an analytic supremum or an implicit-solver stability bound. Smoke configs
write it at every step.

## B7.2 mesh and startup preparation

The optional nonnegative `clustering` parameter is zero by default. For positive
`lambda`, map the uniform vertex heights by
`z=H*sinh(lambda*z_uniform/H)/sinh(lambda)` and retain affine triangles.
An exponential/expm1 form avoids overflow in this expression; collapsed rows
are rejected. Mapping precedes all spaces/geometry-dependent assembly. The
last cell at each wall must fit inside `wall_buffer`; increase nz if clustering
makes the wall cells too large. Uniform meshes retain their previous behavior,
including measured compatibility diagnostics for very coarse meshes.

Metadata stores the actual vertical vertices and global ranges of absolute cell
determinant, vertical-map determinant relative to the uniform mesh, shortest
edge, and `h_max²/(2 area)`, the longest-edge/shortest-altitude aspect ratio.
Triangle reference orientations in Firedrake can be either sign; orientation
preservation follows from ordered mapped rows, not the sign of a raw reference
determinant. The default uniform and tested lambda=2 meshes pass these checks.

Set `[output].horizontal_profiles_every` to a positive step interval for
`horizontal_profiles.csv` and quadrature metadata in `horizontal_profiles.json`.
The optional `Problem.horizontal_profile_geometry` describes an x-periodic
structured triangle channel starting at x=0; only cases supplying that geometry
can request the output. Means are taken on every mesh row and midrow, including
t=0 and the final/early-stop state. Two Gauss points per horizontal half-cell
integrate CG1/CG2 traces exactly; higher degree uses enough Gauss points for its
trace polynomial. All queried points are gathered/broadcast collectively by
Firedrake's `PointEvaluator`, so the means are global rather than rank-local.
The line quadrature is independent postprocessing; it does not alter the
common volume quadrature or OSS projections.
Profile CSV names `step`, `time`, `z`, `velocity_x` and `velocity_y` are reserved;
colliding scalar names are rejected for this output rather than overwriting a column.

The B7.2 study uses six matched pairs: uniform 16×24, 32×48 and 48×72 P1 meshes;
32×48 with half dt; 16×24 P2; and lambda=2 on 32×48 with moisture. Dynamic
nonlinear OSS, BDF2 with lower-order startup, and epsilon=0/1e-4 are compared to
t=0.5. No body force, changed seed amplitude or altered closure removes startup
motion. On these meshes zero-seed energy still dominates the seeded energy;
the report does not certify a KH fitting window or the proposed production mesh.
The analysis-only script reads completed outputs, matches numerical/physical
setups, and reports raw energy ratios and maximum mean-profile changes at the
recorded heights. It neither solves equations nor subtracts energies as if that
gave perturbation energy. See [the B7.2 report](../../../docs/verification/kelvin_helmholtz_b72.md).
