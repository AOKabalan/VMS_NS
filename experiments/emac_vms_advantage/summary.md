# EMAC–VMS advantage experiments

## E1 — translating Gresho vortex (energy and momentum invariants)

| Variant | Form | Status | E(T)/E(0) | ‖u_h − u‖(T) | centre error (T) | ∫𝒩 dt / E(0) | max |M_x − 1| |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| oss | Nonconservative | ok | 0.9838 | 0.1669 | 0.0149 | +2.00e-04 | 4.1e-13 |
| oss | Conservative | ok | 0.9845 | 0.1757 | 0.0096 | -2.53e-04 | 1.8e-15 |
| oss | Skew-symmetric 1 | ok | 0.9841 | 0.1717 | 0.0122 | +0.00e+00 | 4.4e-15 |
| oss | Skew-symmetric 2 | ok | 0.9841 | 0.1717 | 0.0122 | -2.98e-16 | 4.0e-15 |
| oss | EMAC | ok | 0.9920 | 0.1782 | 0.0145 | -8.98e-03 | 4.3e-03 |
| oss | EMAC–VMS | ok | 0.9832 | 0.1780 | 0.0112 | +4.20e-16 | 1.8e-15 |
| oss_weak | Nonconservative | ok | 0.9986 | 0.2595 | 0.0375 | +5.85e-06 | 3.7e-13 |
| oss_weak | Conservative | ok | 0.9987 | 0.2555 | 0.0380 | -7.17e-06 | 1.8e-15 |
| oss_weak | Skew-symmetric 1 | ok | 0.9987 | 0.2574 | 0.0377 | +0.00e+00 | 1.3e-15 |
| oss_weak | Skew-symmetric 2 | ok | 0.9987 | 0.2574 | 0.0377 | -9.62e-17 | 1.3e-15 |
| oss_weak | EMAC | ok | 0.9995 | 0.2500 | 0.0361 | -1.06e-03 | 5.2e-04 |
| oss_weak | EMAC–VMS | ok | 0.9985 | 0.2503 | 0.0370 | -8.06e-17 | 1.6e-15 |
| asgs | Nonconservative | ok | 0.9838 | 0.1686 | 0.0152 | +2.00e-04 | 5.3e-05 |
| asgs | Conservative | ok | 0.9845 | 0.1770 | 0.0099 | -2.53e-04 | 5.7e-05 |
| asgs | Skew-symmetric 1 | ok | 0.9842 | 0.1732 | 0.0125 | +0.00e+00 | 5.5e-05 |
| asgs | Skew-symmetric 2 | ok | 0.9842 | 0.1732 | 0.0125 | +2.22e-18 | 5.5e-05 |
| asgs | EMAC | ok | 0.9916 | 0.1750 | 0.0136 | -8.56e-03 | 4.1e-03 |
| asgs | EMAC–VMS | ok | 0.9832 | 0.1785 | 0.0114 | +2.37e-16 | 7.7e-05 |

## E3 — Gresho in a large box (angular-momentum invariant)

| Variant | Form | Status | A(T)/A(0) − 1 | ∫D_φ dt / A(0) | E(T)/E(0) | ‖u_h − u‖(T) | frame ‖u_h‖ (T) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| oss | Nonconservative | ok | +1.04e+00 | +2.90e-02 | 0.7246 | 0.5210 | 1.25e-01 |
| oss | Conservative | ok | -1.55e-02 | -2.97e-18 | 0.7586 | 0.2934 | 3.00e-02 |
| oss | Skew-symmetric 1 | ok | -9.53e-02 | +1.57e-02 | 0.7452 | 0.3433 | 3.71e-02 |
| oss | Skew-symmetric 2 | ok | -9.53e-02 | +1.57e-02 | 0.7452 | 0.3433 | 3.71e-02 |
| oss | EMAC | ok | +6.99e-02 | -3.72e-02 | 0.6638 | 0.3959 | 3.97e-02 |
| oss | EMAC–VMS | ok | -6.09e-02 | -1.47e-16 | 0.6013 | 0.3914 | 3.56e-02 |
| oss_weak | Nonconservative | ok | +3.95e-01 | +8.11e-03 | 0.9594 | 0.5593 | 1.86e-01 |
| oss_weak | Conservative | ok | -1.05e+00 | +5.59e-19 | 0.9698 | 0.5671 | 1.74e-01 |
| oss_weak | Skew-symmetric 1 | ok | -1.19e+00 | +3.66e-03 | 0.9589 | 0.5622 | 1.87e-01 |
| oss_weak | Skew-symmetric 2 | ok | -1.19e+00 | +3.66e-03 | 0.9589 | 0.5622 | 1.87e-01 |
| oss_weak | EMAC | ok | -4.90e-01 | -6.09e-03 | 0.9556 | 0.5610 | 1.79e-01 |
| oss_weak | EMAC–VMS | ok | -2.56e-01 | -2.98e-16 | 0.9368 | 0.5584 | 1.77e-01 |

## E2 — thin double shear layer (under-resolved robustness)

| Run | Status | E(T)/E(0) | Z(T)/Z(0) | ∫𝒩 dt / E(0) |
| --- | --- | ---: | ---: | ---: |
| E2_conservative_oss_32 | ok | 0.9245 | 0.4111 | -1.26e-02 |
| E2_conservative_oss_64 | ok | 0.9449 | 0.3514 | -4.62e-03 |
| E2_emac_oss_32 | ok | 0.8261 | 0.2972 | -6.89e-03 |
| E2_emac_oss_64 | ok | 0.9036 | 0.2925 | -4.37e-03 |
| E2_emac_vms_oss_32 | ok | 0.8196 | 0.2934 | -1.36e-17 |
| E2_emac_vms_oss_64 | ok | 0.8995 | 0.2898 | +1.03e-16 |
| E2_nonconservative_oss_32 | ok | 0.9158 | 0.3979 | +5.59e-03 |
| E2_nonconservative_oss_64 | ok | 0.9398 | 0.3461 | +2.63e-03 |
| E2_skew_symmetric_1_oss_32 | ok | 0.9187 | 0.4002 | +0.00e+00 |
| E2_skew_symmetric_1_oss_64 | ok | 0.9420 | 0.3456 | +0.00e+00 |
| E2_skew_symmetric_2_oss_32 | ok | 0.9187 | 0.4002 | -6.64e-19 |
| E2_skew_symmetric_2_oss_64 | ok | 0.9420 | 0.3456 | +4.86e-17 |
