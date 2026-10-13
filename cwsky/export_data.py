"""Export the data files the static JavaScript app needs into docs/data/.

  ephemeris.json : Earth barycentric position (AU) and velocity (m/s), daily samples
  detectors.json : site location (m, Earth-fixed) and response tensor per detector
"""

import json
import os
import numpy as np
import lal
from astropy.time import Time

from . import emulator as em

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "data")


def export_ephemeris(start="2010-01-01", end="2040-01-01", step_s=86400):
    t0 = Time(start, scale="utc").gps
    t1 = Time(end, scale="utc").gps
    n = int((t1 - t0) // step_s) + 1
    t = t0 + step_s * np.arange(n)
    pos, vel = em.earth_posvel(t)
    data = {
        "t0": float(t0), "dt": step_s, "n": n,
        "pos_au": [round(float(v), 7) for v in pos.ravel()],
        "vel_ms": [round(float(v), 3) for v in vel.ravel()],
    }
    with open(os.path.join(OUT, "ephemeris.json"), "w") as f:
        json.dump(data, f, separators=(",", ":"))
    return n


def export_detectors():
    out = {}
    for name in em.DETECTOR_INDEX:
        d = em.detector(name)
        out[name] = {
            "name": d.frDetector.name,
            "lat_deg": float(np.degrees(d.frDetector.vertexLatitudeRadians)),
            "lon_deg": float(np.degrees(d.frDetector.vertexLongitudeRadians)),
            "location_m": [float(x) for x in d.location],
            "response": [[float(x) for x in row] for row in np.asarray(d.response)],
        }
    with open(os.path.join(OUT, "detectors.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    n = export_ephemeris()
    export_detectors()
    print(f"wrote {n} ephemeris samples and {len(em.DETECTOR_INDEX)} detectors to {OUT}")
