# Published wave and cloud demos

These are selected artifacts from the completed local simulation runs. The
MP4 files and PNG posters are copied unchanged; the GIF previews are resized
to 840 pixels wide at 8 frames/s for the repository README. Full raw results
remain under the ignored `results/` directory.

| Demo | Video | Preview | Poster | Details |
| --- | --- | --- | --- | --- |
| Stable gravity waves | [MP4](gravity-wave-clouds.mp4) | [GIF](gravity-wave-clouds.gif) | [PNG](gravity-wave-clouds.png) | [Report](../verification/gravity_wave_cloud_pilot.md) |
| Kelvin–Helmholtz billows | [MP4](kelvin-helmholtz-clouds.mp4) | [GIF](kelvin-helmholtz-clouds.gif) | [PNG](kelvin-helmholtz-clouds.png) | [Report](../verification/kelvin_helmholtz_cloud_pilot.md) |

The top panels show a diagnostic cloud indicator over temperature contours;
the lower panels show air rotation. Both are idealized 2D cases without ice
microphysics or latent heating. The KH cloud signal is affected by numerical
moisture overshoot, as its report documents.

[manifest.json](manifest.json) records source paths, configurations, video
metadata and SHA-256 checksums. Diagnostic plots are published beside their
verification reports. The simulation/post-processing reproduction commands
are in those reports.

To recreate a GIF preview from its MP4, run from this directory:

```bash
ffmpeg -i gravity-wave-clouds.mp4 \
  -filter_complex 'fps=8,scale=840:-1:flags=lanczos,split[a][b];[a]palettegen=stats_mode=full[p];[b][p]paletteuse=dither=bayer:bayer_scale=3' \
  -loop 0 gravity-wave-clouds.gif
```
