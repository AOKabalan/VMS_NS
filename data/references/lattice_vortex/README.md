# Lattice vortex reference

Kabalan, MSc thesis, §6.2, equations (6.3)–(6.5), printed pp.44–45 of
`Kabalan_Ali_Thesis.pdf`; read by PDF text extraction on 2026-10-22.
On (0,1)², u=(sin(2πx)sin(2πy), cos(2πx)cos(2πy)) exp(-8νπ²t),
p=-0.5(sin²(2πx)+cos²(2πy)) exp(-16νπ²t), ν=7e-5.
Analytic integration gives E(t)=0.25 exp(-16νπ²t).
The thesis uses h=0.02, linear elements, dt=0.01, T=10, Crank–Nicolson,
and strong Dirichlet data. Its unstructured mesh has 5678 triangles;
our requested 50×50 structured mesh has 5000 triangles.
Re-checked on 2026-10-22 against the thesis text: the formulas, ν, domain, h, element count,
δt, T and the scheme agree. The thesis used a fixed 3 Newton iterations per step; our runs
iterate every step to the solver tolerances.
