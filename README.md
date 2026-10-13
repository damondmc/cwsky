# cwsky

Interactive emulator for continuous gravitational wave (CW) signals from
spinning neutron stars.

## Static web app (GitHub Pages)

`docs/` is a self-contained JavaScript version of the emulator. All physics runs
in the browser, so sliders update live. Preview it locally with

```
python3 -m http.server 8001 --directory docs
```

and open http://127.0.0.1:8001. To publish, enable GitHub Pages for this repo
with source "Deploy from a branch", branch `main`, folder `/docs`. The site then
appears at https://damondmc.github.io/cwsky/.

`docs/data/` holds the inputs the page needs: a daily Earth ephemeris table,
detector geometry exported from LAL, and coastlines. Regenerate the first two
with

```
uv run python -m cwsky.export_data
```

## Python Dash app

The original server-side version, same physics, useful as the reference:

```
uv sync
uv run app.py
```

Open http://127.0.0.1:8050.

- Left panel: CW parameters (sliders and boxes are linked).
- Spectrogram: drag to zoom, scroll to zoom, double-click to reset. The
  emulated signal is re-rendered at the resolution of the zoomed window, so
  zooming in shows the daily Doppler wiggle and zooming out shows the annual
  pattern and the spindown slope.
- Time slider: moves the dashed marker and drives the two 3D views.
- Earth view: drag to rotate, scroll to zoom. Red arrow is the incoming CW
  propagation direction, orange points to the Sun. The globe spins with
  sidereal time.

## Physics

`cwsky/emulator.py` computes everything analytically, no SFTs are generated:

- Doppler: f(t) = [F0 + F1 dt + F2 dt^2 / 2] (1 + v.n / c), with v the
  detector's barycentric velocity (Earth orbit from the JPL DE432s ephemeris
  via astropy, plus the site's GCRS rotation velocity).
- Antenna pattern F+, Fx from lal.ComputeDetAMResponse at the Greenwich mean
  sidereal time, giving the daily amplitude modulation.
- Power per SFT bin relative to the noise: A^2 Tsft / (2 sqrtSX^2), spread
  with a sinc^2 of width 1/Tsft around the track, plus exponential noise.
