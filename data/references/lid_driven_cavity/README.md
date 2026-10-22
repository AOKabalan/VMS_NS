# Ghia centreline data

Ghia, Ghia & Shin (1982), JCP 48, 387–411, Tables I and II, printed pages
398–399 (PDF pages 12–13 of `Ghia1982.pdf`). Transcribed on 2026-10-22 from
the local source PDF. `pdftotext -layout` garbled signs; both table pages
were rendered with `pdftoppm` and every Re=100/1000 entry checked visually.
No values were reconstructed from figures. Coordinates retain the paper’s
rounded 17 stations, in descending order. `u.csv` is u at x=0.5;
`v.csv` is v at y=0.5. Endpoints are u(bottom)=0, u(lid)=1 and v=0 at both walls.
Ghia’s finite-difference accuracy limits comparisons near walls at Re=1000.

Independent re-check, 2026-10-22: both table pages rendered at 400 dpi (`pdftoppm -r 400`),
the y/x, Re = 100 and Re = 1000 columns cropped and read; all 68 values and all 34 stations
agree with `u.csv` and `v.csv`. The commonly quoted extremes agree with the PDF, where the paper
underscores them as local extrema: Re = 100, u min −0.21090 at y = 0.4531, v max 0.17527 at
x = 0.2344, v min −0.24533 at x = 0.8047; Re = 1000, u min −0.38289 at y = 0.1719, v max 0.37095
at x = 0.1563, v min −0.51550 at x = 0.9063.
