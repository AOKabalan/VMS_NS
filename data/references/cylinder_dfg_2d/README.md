# DFG 2D-2 references

Retrieved 2026-10-22 from
https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html

`featflow.json`: the page’s finest listed level (6, k=1/200), primary scalar
reference. Mean=(max+min)/2 and amplitude=max-min, as defined on the page.
`dfg_1996.json`: contextual intervals transcribed from Schäfer & Turek (1996),
Table 4, printed p.13 of the supplied `dfg_benchmark_results.pdf`.

All five linked ZIPs are preserved unchanged and pass ZIP CRC checks. The four
older archives come from `~featflow/media/dfg_bench2_2d/`; the 20141124 archive
comes from `~featflow/media/dfg_bench2new_2d/`. The latter uses higher-order
line integration for forces. `featflow_curve.csv` extracts columns 2,4,5
(time, C_D, C_L) of `logs/results_cc2d_bench2_referencenew/bdforces_q2_lv6_dt4`,
following that archive’s plot_drag/plot_lift scripts (level 6, dt=1/1600).
Scalar page references and newer curves are kept distinct; the latter need
not reproduce the older page values exactly. Curves are compared after
aligning their last complete cycles at a lift minimum, without rescaling time.

SHA-256 of downloaded archives:
- `20141124_dfg_bench2_official.zip`: `08db59605b2f60efc512d2bbc79628ea5d1cc720c0986938a16920884a0215be`
- `draglift_q1t_cn_lv3-6_dt1-4.zip`: `9a9cc81374236b08fdbe5ffc80d2ba2dbad76e7d032d042f14639681b3570202`
- `draglift_q2_cn_lv3-6_dt1-4.zip`: `a193c85f19966f5fd0935dfa92f30d646ea33b06dbf9942b74dcf67365ce6f6e`
- `pressure_q1t_cn_lv3-6_dt1-4.zip`: `6e83193dd98d5a3b01b38d3976fc55fb017aff47fe3890808eff1497b2eba6c4`
- `pressure_q2_cn_lv3-6_dt1-4.zip`: `e187134f084e614db4e25cdbe174677f969b4bb45afd1a5ea7cd8ba3e340b43e`
