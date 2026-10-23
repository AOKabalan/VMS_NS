# Chapter 5 renumbering after adding the consistent-splitting section

- **Old chapter:** `720be759…`, 19 pages. This is the version `formulation_spec.md` cites, as of Phase 11b and 12a.
- **New chapter:** SHA-256 `c4995958c5e410591261f676aae661f5854def033feb534dc930e40668491784`, 23 pages.

## What changed

- **New §5.7, "Consistent splitting of the convective term"** (theory only, no new equations for the code):
  - **Lemma 5.3:** the standard coupling −(ũ, (a·∇)v_h) is exactly the conservative form's fine-scale part.
  - **Propositions 5.6–5.7:** momenta and energy when every form is split consistently.
  - **Table 5.1:** the defects of the consistent splitting.
  - **Remarks 5.8–5.9:** comparison with the fixed coupling, and practicality.
- **Moved:** the symmetric-stress lemma moved from the EMAC–VMS section into §5.7 (now Lemma 5.4). It gained the identity c_E(a; u, u) = ((u·∇)u, a) + ½(∇·a, |u|²), equation (5.24).
- **Unchanged and reframed:**
  - §5.3–5.6 are labelled as the "fixed coupling".
  - EMAC–VMS (now §5.8) is presented as the consistent splitting of EMAC, with the subscale kinetic pressure moved into p̃.
  - The §5.8.5 variant is stated to be exactly the consistent splitting.
  - No equation of the implemented formulation changed.

## Numbers that moved (old → new)

| Item | Old | New |
|---|---|---|
| Fixed-coupling defect tables | Tables 5.1, 5.2, 5.3 | Tables 5.2, 5.3, 5.4 |
| Consistent-splitting table | — | Table 5.1 |
| EMAC–VMS section | §5.7 | §5.8 |
| Momenta / variant / implementation subsections | §5.7.3, §5.7.6, §5.9.3 | §5.8.2, §5.8.5, §5.10.3 |
| Summary, implementation | §5.8, §5.9 | §5.9, §5.10 |
| Symmetric-stress lemma | Lemma 5.3, (5.22) | Lemma 5.4, (5.23) |
| S(w) | (5.21) | (5.22) |
| Resolved equation, G | (5.24), (5.25) | (5.28), (5.29) |
| Residual | (5.26) | (5.30) |
| T̃ | (5.23) | (5.27) |
| G as stress, expansion | (5.27), (5.28) | (5.31), (5.32) |
| Consistency equations | (5.29), (5.30) | (5.33), (5.34) |
| EMAC–VMS energy balance | (5.31) | (5.35) |
| Linear splitting | (5.32) | (5.36) |
| Variant T̃_full, N_E, N_E − N | (5.33)–(5.35) | (5.37)–(5.39) |
| DS, polarization, symmetric Picard | (5.36)–(5.38) | (5.40)–(5.42) |
| Picard, Newton (implementation) | (5.39), (5.40) | (5.43), (5.44) |
| Subgrid-stress element matrix | (5.41) | (5.45) |
| Consistency proposition | Prop. 5.6 | Prop. 5.8 |
| Momenta, energy, linear, iteration | Props. 5.7, 5.8, 5.9, 5.10 | Props. 5.9, 5.10, 5.11, 5.12 |
| Interpretation remark | Remark 5.8 | Remark 5.10 |
| Eliminated subscale | Remark 5.9 | Remark 5.11 |
| Time integration, stabilization | Remarks 5.10, 5.11 | Remarks 5.12, 5.13 |
| Boundary conditions and pressure | Remark 5.12 | Remark 5.14 |
| Related work | Remark 5.13 | Remark 5.15 |
| Picard remark (implementation) | Remark 5.14 | Remark 5.16 |

Unchanged:
- (5.1)–(5.20);
- Lemma 5.2;
- Remarks 5.1–5.7, including the projection-assumption remark 5.6;
- Propositions 5.1–5.5;
- Corollary 5.1;
- Definition 5.1.

## For the agent

Phase 12a's code is unaffected: no implemented equation changed. When the 12b spec phase starts:
1. Point the spec at the new hash.
2. Update its chapter-reference table using the map above.

The consistent splitting needs no new code:
- for the conservative form it coincides with the fixed coupling;
- for the non-conservative form it conserves nothing;
- for the skew-symmetric forms it needs ∇·ũ in the pointwise subscale residual;
- for EMAC it is the §5.8.5 variant.
