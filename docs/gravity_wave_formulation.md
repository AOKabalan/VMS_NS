# Stable internal gravity waves and a moist-layer visualization

This application uses the existing thermal Boussinesq equations and passive
scalar transport. It illustrates an air-motion mechanism associated with
undulatus patterns. It does not identify a cloud genus or simulate high-cloud
ice microphysics. WMO defines
[undulatus](https://cloudatlas.wmo.int/en/clouds-varieties-undulatus.html) by its
wavy appearance and separately defines KH breaking-wave curls as
[fluctus](https://cloudatlas.wmo.int/en/clouds-supplementary-features-fluctus.html).

## Background, boundaries and linear reference

Use the periodic channel [0,L]×[-H/2,H/2], L=n(2π/k), with free-slip rigid
horizontal walls. There is no imposed mean flow or velocity shear. The
unperturbed state is

$$
u_0=0,\qquad T_0=\Gamma z,\qquad p_0=\tfrac12\beta\Gamma z^2,
\qquad g=(0,-1),\quad T_{ref}=0,\quad N^2=\beta\Gamma>0.
$$

Hold T=Γz at the two walls. An insulated wall would conflict with this
nonzero background gradient; fixed wall values maintain the conductive
equilibrium ΔT_0=0. Velocity has zero normal component and natural zero
tangential derivative. Moisture has natural zero diffusive flux. Pressure
retains its constant null space.

Set m=π/H, K²=k²+m² and equal momentum/thermal diffusivity μ. The linearized
equations about the background give

$$
\partial_t u'+\nabla p'-\mu\Delta u'=\beta T'e_z,\quad
\nabla\cdot u'=0,\quad
\partial_t T'+\Gamma w-\mu\Delta T'=0,
\qquad \omega^2=N^2\frac{k^2}{k^2+m^2}.
$$

This is the nonrotating, nonhydrostatic Boussinesq dispersion relation in
[MIT's internal-gravity-wave notes, pp. 1–2](https://ocw.mit.edu/courses/12-810-dynamics-of-the-atmosphere-spring-2023/mit12_810_s23_intro_igw.pdf).
It is also verified directly against the differential equations by the case's
analytic test. Define d=μK² and φ=cos(kx)cos(mz). A standing-wave reference is

$$
\begin{aligned}
u'_x&=A\omega\frac{m}{k}e^{-dt}\cos(\omega t)\sin(kx)\sin(mz),\\
w&=A\omega e^{-dt}\cos(\omega t)\phi,\\
T'&=-\Gamma A e^{-dt}\sin(\omega t)\phi,\\
p'&=-A\omega^2\frac{m}{k^2}e^{-dt}\sin(\omega t)\cos(kx)\sin(mz).
\end{aligned}
$$

The code initializes this reference at t=0: zero temperature displacement and
maximum wave velocity. No time-dependent forcing holds the wave. The solver
retains nonlinear advection, so the linear fields are not registered as an exact
solution of that nonlinear problem. For the linear reference the minimum
vertical temperature gradient is at least Γ(1-mA), positive for the chosen
amplitudes. The application must check the actual FE gradient, rather than
assuming this inequality survives discretization.

P2 velocity/pressure/temperature allow exact representation of the quadratic
hydrostatic background pressure. Coarse P1 trials develop startup motion;
P2's two-step zero-kick control preserves rest to the specified tolerance.
This is a case/discretization choice, with no background-force subtraction or
change to the momentum equation.

## Moisture and the cloud proxy

Advect and diffuse q with its own positive diffusivity κ_q and zero source,
without feedback into buoyancy. Initially

$$
q_0=RH_*\exp[-(z/\sigma)^2-\gamma z]\cos^2(mz),\qquad
q_{sat}=e^{-\gamma z},\qquad q_c=\max(q-q_{sat},0).
$$

Thus initial relative humidity is at most RH_*<1 and q has zero normal
derivative at both walls. The case records the continuum initial tracer
maximum by locating the unique root of
`-2z/σ²-γ-2m tan(mz)=0`. The proxy removes no moisture and releases no latent
heat. Height-dependent saturation represents a qualitative lifting threshold;
it does not calculate adiabatic cooling or ice formation.

For an initially central, nondiffusing parcel, the upward displacement required
for saturation is -ln(RH_*)/γ. With RH_*=0.98 and γ=0.3 this is about 0.0673,
below the proposed A=0.2 kick. This motivates the initial data before running;
it is not a guarantee of cloud formation after mixing and diffusion. Cloud
should correlate with upward displacement, not necessarily instantaneous
upward velocity: velocity and displacement have a quarter-period offset.

The visualization evaluates the original P1/P2 polynomials at triangle
quadrature points. It reports nodal and quadrature tracer extrema, raw cloud
fraction and quadrature sensitivity. P2 display images use four subtriangles
per parent cell. Native VTK q_c is a nodal interpolation and may have negative
values between P2 nodes; the transported q is never clipped. All quantitative
cloud diagnostics use the raw positive-part expression.

## Bounded execution and checks

The application config is `configs/boussinesq/gravity_wave_cloud_pilot.toml`.
The smaller `gravity_wave_period_check.toml` checks the low-amplitude reference.
The pilot uses four wavelengths, P2 on 64×16 cells, A=0.2, Γ=β=k=1, H=2π,
RH_*=0.98, σ=1.2 and γ=0.3. Its period is 2π/√0.8≈7.02481. BDF2 uses
80 steps per reference period for three periods and saves every two steps.
Momentum and thermal diffusion are μ=0.01; moisture diffusion is 0.0002.
Dynamic linear OSS and skew-symmetric scalar convection use existing kernels.
These are dimensionless illustrative choices, separate from the KH parameters.

Required checks are hydrostatic-rest preservation, low-amplitude modal period
and amplitude against the damped reference, positive actual ∂zT, tracer bounds,
cloud correlation with lifted air and consistent saved-field integrals.
The small-amplitude run requires period error below 2% and both modal RMS
errors below 5% of the initial respective amplitudes, across three periods.
The current thermal-content-budget group rejects fixed-temperature walls even
when q is insulated, so this config does not request that group. Scalar content
is still recorded. No insulated-domain KE+potential identity is asserted across
the fixed-temperature walls.

Very-low-diffusion P2 trials at μ=0.0002, with conservative and skew-symmetric
scalar convection, follow the initial wave but subsequently develop short-scale
temperature oscillations. Their completed solves are rejected as stable-wave
verification. The more damped μ=0.01 setup is a distinct physical parameter
choice whose actual stability and period must still pass the checks above.
The report retains unsuccessful trials and separates them from accepted data.
