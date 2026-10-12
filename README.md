# cwsky

Interactive emulator for continuous gravitational wave (CW) signals from
spinning neutron stars.

## Run

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
