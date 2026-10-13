"""Fast CW signal emulator: Doppler track, antenna pattern, synthetic spectrogram.

No SFT generation. Everything is computed analytically from the ephemeris and
detector geometry, so it is fast enough to re-render interactively.
"""

import numpy as np
import lal
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    solar_system_ephemeris,
    get_body_barycentric_posvel,
)
from astropy import units as u
from astropy.utils import iers

# Never touch the network: UT1-UTC is below 1 s, irrelevant at this level.
iers.conf.auto_download = False
iers.conf.iers_degraded_accuracy = "ignore"

solar_system_ephemeris.set("de432s")

C = 299_792_458.0  # m/s

DETECTOR_INDEX = {
    "H1": lal.LALDetectorIndexLHODIFF,
    "L1": lal.LALDetectorIndexLLODIFF,
    "V1": lal.LALDetectorIndexVIRGODIFF,
    "K1": lal.LALDetectorIndexKAGRADIFF,
}


def detector(name: str) -> lal.Detector:
    return lal.CachedDetectors[DETECTOR_INDEX[name]]


def source_unit_vector(alpha: float, delta: float) -> np.ndarray:
    """Unit vector from the SSB toward the source (equatorial ICRS frame)."""
    return np.array([
        np.cos(delta) * np.cos(alpha),
        np.cos(delta) * np.sin(alpha),
        np.sin(delta),
    ])


def earth_posvel(t_gps: np.ndarray) -> tuple:
    """Earth barycentric position (AU) and velocity (m/s), shapes (N, 3)."""
    times = Time(np.atleast_1d(t_gps), format="gps", scale="utc")
    pos, vel = get_body_barycentric_posvel("earth", times.tdb)
    p = pos.xyz.to(u.AU).value.T
    v = vel.xyz.to(u.m / u.s).value.T
    return p, v


def site_posvel(t_gps: np.ndarray, det: lal.Detector) -> tuple:
    """Detector position (m) and velocity (m/s) relative to Earth centre, GCRS, shapes (N, 3)."""
    times = Time(np.atleast_1d(t_gps), format="gps", scale="utc")
    loc = EarthLocation.from_geocentric(*det.location, unit=u.m)
    pos, vel = loc.get_gcrs_posvel(times)
    return pos.xyz.to(u.m).value.T, vel.xyz.to(u.m / u.s).value.T


def detector_velocity(t_gps: np.ndarray, det: lal.Detector, include_spin: bool = True) -> np.ndarray:
    """Barycentric velocity of the detector (m/s), shape (N, 3)."""
    _, v_orb = earth_posvel(t_gps)
    if not include_spin:
        return v_orb
    _, v_spin = site_posvel(t_gps, det)
    return v_orb + v_spin


def doppler_track(
    t_gps: np.ndarray,
    F0: float,
    F1: float,
    F2: float,
    alpha: float,
    delta: float,
    tref: float,
    det: lal.Detector,
    include_spin: bool = True,
) -> np.ndarray:
    """Observed frequency f(t) = f_intrinsic(t) * (1 + v.n / c)."""
    n = source_unit_vector(alpha, delta)
    beta = detector_velocity(t_gps, det, include_spin) @ n / C
    dt = t_gps - tref
    f_int = F0 + F1 * dt + 0.5 * F2 * dt**2
    return f_int * (1.0 + beta)


def gmst_rad(t_gps: np.ndarray) -> np.ndarray:
    times = Time(np.atleast_1d(t_gps), format="gps", scale="utc")
    return times.sidereal_time("mean", "greenwich").rad


def antenna_pattern(t_gps: np.ndarray, alpha: float, delta: float, psi: float, det: lal.Detector) -> tuple:
    """F+(t), Fx(t) for the detector, shapes (N,)."""
    g = gmst_rad(t_gps)
    fp = np.empty(len(g))
    fc = np.empty(len(g))
    for i, gi in enumerate(g):
        fp[i], fc[i] = lal.ComputeDetAMResponse(det.response, alpha, delta, psi, gi)
    return fp, fc


def signal_amplitude_sq(fp: np.ndarray, fc: np.ndarray, h0: float, cosi: float) -> np.ndarray:
    """Squared strain amplitude seen by the detector: (A+ F+)^2 + (Ax Fx)^2."""
    a_plus = h0 * (1.0 + cosi**2) / 2.0
    a_cross = h0 * cosi
    return (a_plus * fp) ** 2 + (a_cross * fc) ** 2


def emulate_spectrogram(
    params: dict,
    t_range: tuple,
    f_range: tuple,
    n_t_max: int = 500,
    n_f_max: int = 400,
    noise: bool = True,
    seed: int = 0,
) -> dict:
    """Render an emulated power spectrogram over a time/frequency window.

    params keys: F0, F1, F2, Alpha, Delta, h0, cosi, psi, sqrtSX, Tsft,
                 tstart, tref, detector, include_spin.
    t_range : (t0, t1) GPS seconds of the window.
    f_range : (f0, f1) Hz of the window.

    Time samples are SFT midpoints inside the window (decimated to n_t_max).
    Frequency samples are SFT bins (decimated to n_f_max). When decimated,
    the signal line width is widened to the pixel size so it stays visible.

    Returns dict with t (GPS), f (Hz), power (Nt, Nf) in units of mean noise
    power, track (Nt,), amp2 (Nt,) relative signal power per bin.
    """
    Tsft = params["Tsft"]
    tstart = params["tstart"]
    det = detector(params["detector"])

    # SFT midpoints covering the window
    k0 = int(np.floor((t_range[0] - tstart) / Tsft))
    k1 = int(np.ceil((t_range[1] - tstart) / Tsft))
    k = np.arange(max(k0, 0), max(k1, k0 + 1) + 1)
    if len(k) > n_t_max:
        k = np.unique(np.linspace(k[0], k[-1], n_t_max).astype(int))
    t = tstart + (k + 0.5) * Tsft

    # Frequency bins covering the window
    dF = 1.0 / Tsft
    j0 = int(np.floor(f_range[0] / dF))
    j1 = int(np.ceil(f_range[1] / dF))
    j = np.arange(j0, max(j1, j0 + 1) + 1)
    if len(j) > n_f_max:
        f = np.linspace(j[0] * dF, j[-1] * dF, n_f_max)
    else:
        f = j * dF
    df_pix = f[1] - f[0] if len(f) > 1 else dF
    width = max(dF, df_pix)

    track = doppler_track(
        t, params["F0"], params["F1"], params["F2"],
        params["Alpha"], params["Delta"], params["tref"], det,
        params.get("include_spin", True),
    )
    fp, fc = antenna_pattern(t, params["Alpha"], params["Delta"], params["psi"], det)
    a2 = signal_amplitude_sq(fp, fc, params["h0"], params["cosi"])
    # Signal power per SFT bin relative to mean noise power
    amp2 = a2 * Tsft / (2.0 * params["sqrtSX"] ** 2)

    x = (f[None, :] - track[:, None]) / width
    power = amp2[:, None] * np.sinc(x) ** 2
    if noise:
        rng = np.random.default_rng(seed)
        power = power + rng.exponential(1.0, size=power.shape)

    return {"t": t, "f": f, "power": power, "track": track, "amp2": amp2, "fp": fp, "fc": fc}
