# cwsky: CW Signal Simulator

Interactive simulator for continuous gravitational wave (CW) signals from
non-axisymmetric spinning neutron stars. Set the frequency, spindown, sky
position and detector, and watch the Doppler-modulated signal appear in a
frequency-time spectrogram. Zoom from the annual pattern down to the daily
wiggle, and follow Earth and the incoming wave direction in a linked 3D view.

**Play online:** https://damondmc.github.io/cwsky/

## Run locally

The web app in `docs/` is self-contained JavaScript. All physics runs in the
browser, so sliders update live. Serve the folder and open http://127.0.0.1:8001:

```
python3 -m http.server 8001 --directory docs
```

Controls:

- Left panel: CW parameters (sliders and boxes are linked).
- Spectrogram: drag to zoom, scroll to zoom, double-click to reset. The signal
  is re-rendered at the resolution of the zoomed window, so zooming in shows the
  daily Doppler wiggle and zooming out shows the annual pattern and the spindown
  slope.
- Time slider: moves the dashed marker and drives the two 3D views.
- Earth view: drag to rotate, scroll to zoom. Red arrow is the incoming CW
  propagation direction, orange points to the Sun. The globe spins with
  sidereal time.

## Physics

The observed frequency is

    f(t) = [f0 + f1 dt + f2 dt^2 / 2] (1 + (v/c).n)

with v the detector's barycentric velocity (Earth's orbit from the JPL DE432s
ephemeris plus the site's rotation about Earth's axis) and n the unit vector to
the source. The antenna pattern F+, Fx follows the same algorithm as
`lal.ComputeDetAMResponse` at the Greenwich mean sidereal time, giving the
daily amplitude modulation. Power per SFT bin relative to the noise is
A^2 Tsft / (2 Sn), spread with a sinc^2 of width 1/Tsft around the track, plus
exponential noise. No SFTs are generated.

## Python reference implementation

`cwsky/emulator.py` is the same model in Python on top of astropy and LAL, and
`app.py` is a Plotly Dash version of the interface. Both are kept as the
reference the JavaScript port is checked against.

```
uv sync
uv run app.py
```

`docs/data/` holds the inputs the web app needs: a daily Earth ephemeris table
and detector geometry exported from LAL, plus coastlines. Regenerate the first
two with

```
uv run python -m cwsky.export_data
```
