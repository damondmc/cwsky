#!/usr/bin/env python3
"""Interactive CW signal simulator: parameter panel, zoomable f-t spectrogram, 3D Earth view."""

import json
import os

import numpy as np
import plotly.graph_objects as go
from astropy.time import Time
from dash import Dash, dcc, html, Input, Output, State, Patch, ctx, no_update

from cwsky import emulator as em

PI = np.pi
BG = "#d6dae0"
PANEL = "#c6ccd4"
SCENE_BG = "#dde1e7"
FG = "#1a1a1a"
MUTED = "#4a4f57"
INPUT_STYLE = {"background": "#eef0f3", "color": FG, "border": "1px solid #9aa0a8", "fontSize": 12}

with open(os.path.join(os.path.dirname(__file__), "cwsky", "data", "coastlines_110m.json")) as _f:
    COASTLINES = json.load(_f)  # list of segments, each a list of [lon, lat] in degrees

# name, label (markdown with $latex$), slider(min, max, step), default, log10 slider?
PARAMS = [
    ("F0", r"$f_0$ (Hz)", (20.0, 1500.0, 0.1), 100.0, False),
    ("F1", r"$\dot f$ (Hz/s)", (-1e-8, 0.0, 1e-11), -1e-9, False),
    ("F2", r"$\ddot f$ (Hz/s$^2$)", (-1e-17, 1e-17, 1e-19), 0.0, False),
    ("Alpha", r"$\alpha$ right ascension (rad)", (0.0, 2 * PI, 0.01), 1.5, False),
    ("Delta", r"$\delta$ declination (rad)", (-PI / 2, PI / 2, 0.01), 0.5, False),
    ("h0", r"$h_0$ (log slider)", (-25.0, -22.0, 0.05), 1e-23, True),
    ("sqrtSX", r"$\sqrt{S_n}$ (Hz$^{-1/2}$, log slider)", (-24.0, -22.0, 0.05), 1e-23, True),
    ("cosi", r"$\cos\iota$", (-1.0, 1.0, 0.01), 0.0, False),
    ("psi", r"$\psi$ polarisation (rad)", (-PI / 4, PI / 4, 0.01), 0.0, False),
    ("Tsft", r"$T_\mathrm{SFT}$ (s)", (300, 86400, 300), 1800, False),
    ("duration", r"$T_\mathrm{obs}$ (days)", (1, 730, 1), 500, False),
]
DEFAULT_TSTART = 1_400_000_000


def label(text, **style):
    return dcc.Markdown(text, mathjax=True,
                        style={"fontSize": 12, "marginTop": 8, "marginBottom": 2, **style})


def control(name, text, rng, default, log):
    sval = np.log10(default) if log else default
    return html.Div([
        label(text),
        html.Div([
            html.Div(dcc.Slider(
                id=f"sl-{name}", min=rng[0], max=rng[1], step=rng[2], value=sval,
                marks=None, tooltip={"placement": "bottom", "always_visible": False},
                updatemode="mouseup",
            ), style={"flex": 1}),
            dcc.Input(id=f"in-{name}", type="number", value=default, debounce=True,
                      style={"width": 90, **INPUT_STYLE}),
        ], style={"display": "flex", "alignItems": "center", "gap": 6}),
    ])


app = Dash(__name__)
app.title = "CW Signal Simulator"

app.layout = html.Div([
    dcc.Store(id="view"),
    dcc.Store(id="params"),
    html.Div([
        html.H3("CW parameters", style={"margin": "4px 0 8px"}),
        *[control(*p) for p in PARAMS],
        label("Detector"),
        dcc.Dropdown(id="detector", options=list(em.DETECTOR_INDEX), value="H1",
                     clearable=False),
        label("Start time (GPS s)"),
        dcc.Input(id="tstart", type="number", value=DEFAULT_TSTART, debounce=True,
                  style={"width": "100%", **INPUT_STYLE}),
        html.Div(id="tstart-date", style={"fontSize": 11, "color": MUTED}),
        dcc.Checklist(
            id="opts",
            options=[
                {"label": " Earth rotation (daily Doppler)", "value": "spin"},
                {"label": " Show track line", "value": "track"},
                {"label": " Noise", "value": "noise"},
            ],
            value=["spin", "track", "noise"],
            style={"fontSize": 12, "marginTop": 10},
        ),
        dcc.Markdown(
            r"""
Drag on the spectrogram to zoom, double-click to reset. Drag the Earth view to rotate,
scroll to zoom.

$f(t) = \left[f_0 + \dot f\,\Delta t + \tfrac12 \ddot f\,\Delta t^2\right]
\left(1 + \frac{\vec v}{c}\cdot\hat n\right)$
""",
            mathjax=True, style={"fontSize": 11, "color": MUTED, "marginTop": 12},
        ),
    ], style={"width": 330, "padding": 12, "background": PANEL, "overflowY": "auto"}),
    html.Div([
        dcc.Graph(id="spec", style={"height": "46vh"}, config={"scrollZoom": True}, mathjax=True),
        html.Div([
            html.Div("Time (days from start)", style={"fontSize": 12, "width": 150}),
            html.Div(dcc.Slider(id="tslider", min=0, max=365, step=0.01, value=0, marks=None,
                                tooltip={"placement": "bottom"}, updatemode="drag"),
                     style={"flex": 1}),
        ], style={"display": "flex", "alignItems": "center", "padding": "0 12px"}),
        html.Div([
            dcc.Graph(id="earth", style={"flex": 3, "height": "52vh"}, mathjax=True),
            dcc.Graph(id="orbit", style={"flex": 2, "height": "52vh"}, mathjax=True),
        ], style={"display": "flex"}),
    ], style={"flex": 1, "minWidth": 0, "overflowY": "auto"}),
], style={"display": "flex", "height": "100vh", "background": BG, "color": FG,
          "fontFamily": "sans-serif"})


# ---------------------------------------------------------------------------
# Slider <-> input sync
# ---------------------------------------------------------------------------

def _make_sync(name, log):
    @app.callback(
        Output(f"sl-{name}", "value"), Output(f"in-{name}", "value"),
        Input(f"sl-{name}", "value"), Input(f"in-{name}", "value"),
        prevent_initial_call=True,
    )
    def sync(sval, ival):
        if ctx.triggered_id == f"sl-{name}":
            return no_update, float(10 ** sval) if log else sval
        if ival is None:
            return no_update, no_update
        if log:
            return (np.log10(ival) if ival > 0 else no_update), no_update
        return ival, no_update


for _name, _, _, _, _log in PARAMS:
    _make_sync(_name, _log)


@app.callback(Output("tstart-date", "children"), Input("tstart", "value"))
def show_date(tstart):
    if tstart is None:
        return ""
    return Time(tstart, format="gps").utc.iso[:19] + " UTC"


# ---------------------------------------------------------------------------
# Spectrogram
# ---------------------------------------------------------------------------

def full_view(p):
    """Full time range and a frequency window that contains the whole track."""
    t = np.linspace(p["tstart"], p["tstart"] + p["duration"], 400)
    det = em.detector(p["detector"])
    tr = em.doppler_track(t, p["F0"], p["F1"], p["F2"], p["Alpha"], p["Delta"],
                          p["tref"], det, p["include_spin"])
    dF = 1.0 / p["Tsft"]
    lo, hi = tr.min(), tr.max()
    pad = max(0.15 * (hi - lo), 20 * dF)
    return [0.0, p["duration"] / 86400.0], [float(lo - pad), float(hi + pad)]


def spectrogram_figure(p, view, show_track, noise, t_mark):
    t0 = p["tstart"] + view["t"][0] * 86400
    t1 = p["tstart"] + view["t"][1] * 86400
    r = em.emulate_spectrogram(p, (t0, t1), tuple(view["f"]), noise=noise)
    t_days = (r["t"] - p["tstart"]) / 86400.0
    zmax = max(np.percentile(r["power"], 99.5), 1.0)

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        x=t_days, y=r["f"], z=r["power"].T, colorscale="Viridis", zmin=0, zmax=zmax,
        colorbar={"title": "power / noise"},
        hovertemplate="t = %{x:.3f} d<br>f = %{y:.6f} Hz<br>P = %{z:.1f}<extra></extra>",
    ))
    if show_track:
        fig.add_trace(go.Scatter(
            x=t_days, y=r["track"], mode="lines", name="Doppler track",
            line={"color": "#ff4444", "width": 1}, hoverinfo="skip",
        ))
    fig.update_layout(
        shapes=[{"type": "line", "x0": t_mark, "x1": t_mark, "y0": 0, "y1": 1,
                 "yref": "paper", "line": {"color": "white", "dash": "dash", "width": 1.5}}],
        xaxis={"title": "Time from start (days)", "range": view["t"]},
        yaxis={"title": "Frequency (Hz)", "range": view["f"],
               "tickformat": ".6f", "exponentformat": "none"},
        margin={"l": 90, "r": 10, "t": 30, "b": 40},
        paper_bgcolor=BG, plot_bgcolor="#440154", font={"color": FG},
        showlegend=False, dragmode="zoom",
        title={"text": f"{p['detector']}   Tsft = {p['Tsft']:.0f} s   "
                       f"bin = {1e3 / p['Tsft']:.3f} mHz   "
                       f"pixels {len(r['t'])} x {len(r['f'])}",
               "font": {"size": 13}},
    )
    return fig


@app.callback(
    Output("spec", "figure"), Output("view", "data"), Output("params", "data"),
    Output("tslider", "max"), Output("tslider", "value"),
    *[Input(f"in-{n}", "value") for n, *_ in PARAMS],
    Input("detector", "value"), Input("tstart", "value"), Input("opts", "value"),
    Input("spec", "relayoutData"),
    State("view", "data"), State("tslider", "value"),
)
def render_spectrogram(*args):
    n = len(PARAMS)
    vals = args[:n]
    detector, tstart, opts, relayout, view, t_mark = args[n:]
    if any(v is None for v in vals) or tstart is None:
        return (no_update,) * 5

    p = {name: float(v) for (name, *_), v in zip(PARAMS, vals)}
    if p["Tsft"] <= 0 or p["duration"] <= 0 or p["sqrtSX"] <= 0 or p["h0"] < 0 or p["F0"] <= 0:
        return (no_update,) * 5
    p["duration"] *= 86400.0
    p.update(tstart=float(tstart), tref=float(tstart), detector=detector,
             include_spin="spin" in opts)
    dur_days = p["duration"] / 86400.0
    t_mark = min(t_mark or 0.0, dur_days)

    if ctx.triggered_id == "spec":
        if not relayout or view is None:
            return (no_update,) * 5
        new = {"t": list(view["t"]), "f": list(view["f"])}
        full = None
        changed = False
        if "xaxis.range[0]" in relayout:
            new["t"] = [relayout["xaxis.range[0]"], relayout["xaxis.range[1]"]]
            changed = True
        elif "xaxis.range" in relayout:
            new["t"] = list(relayout["xaxis.range"])
            changed = True
        if "yaxis.range[0]" in relayout:
            new["f"] = [relayout["yaxis.range[0]"], relayout["yaxis.range[1]"]]
            changed = True
        elif "yaxis.range" in relayout:
            new["f"] = list(relayout["yaxis.range"])
            changed = True
        if relayout.get("xaxis.autorange") or relayout.get("yaxis.autorange"):
            full = full_view(p)
            if relayout.get("xaxis.autorange"):
                new["t"] = full[0]
            if relayout.get("yaxis.autorange"):
                new["f"] = full[1]
            changed = True
        if not changed:
            return (no_update,) * 5
        new["t"][0] = max(new["t"][0], 0.0)
        new["t"][1] = min(new["t"][1], dur_days)
        view = new
    else:
        t_rng, f_rng = full_view(p)
        view = {"t": t_rng, "f": f_rng}

    fig = spectrogram_figure(p, view, "track" in opts, "noise" in opts, t_mark)
    return fig, view, p, dur_days, t_mark


# ---------------------------------------------------------------------------
# 3D views
# ---------------------------------------------------------------------------

def _sphere(radius, n=40):
    u_ = np.linspace(0, 2 * PI, n)
    v_ = np.linspace(0, PI, n // 2)
    x = radius * np.outer(np.cos(u_), np.sin(v_))
    y = radius * np.outer(np.sin(u_), np.sin(v_))
    z = radius * np.outer(np.ones(n), np.cos(v_))
    return x, y, z


def _solid_sphere(center, radius, color, name):
    x, y, z = _sphere(radius, n=24)
    return go.Surface(
        x=x + center[0], y=y + center[1], z=z + center[2],
        colorscale=[[0, color], [1, color]], showscale=False, name=name, showlegend=True,
        hoverinfo="name", lighting={"ambient": 0.75},
    )


def _arrow(tail, head, color, name, width=5):
    d = np.asarray(head) - np.asarray(tail)
    return [
        go.Scatter3d(x=[tail[0], head[0]], y=[tail[1], head[1]], z=[tail[2], head[2]],
                     mode="lines", line={"color": color, "width": width}, name=name,
                     hoverinfo="name"),
        go.Cone(x=[head[0]], y=[head[1]], z=[head[2]], u=[d[0]], v=[d[1]], w=[d[2]],
                sizemode="absolute", sizeref=0.25 * np.linalg.norm(d), anchor="tip",
                colorscale=[[0, color], [1, color]], showscale=False, hoverinfo="skip"),
    ]


def _scene(uirev, rng):
    ax = {"visible": False, "range": [-rng, rng]}
    return {
        "uirevision": uirev,
        "scene": {"xaxis": ax, "yaxis": ax, "zaxis": ax, "aspectmode": "cube",
                  "bgcolor": SCENE_BG, "uirevision": uirev},
        "paper_bgcolor": BG, "font": {"color": FG},
        "margin": {"l": 0, "r": 0, "t": 30, "b": 0},
        "legend": {"x": 0, "y": 1, "font": {"size": 10}, "bgcolor": "rgba(255,255,255,0.6)"},
    }


def _lonlat_to_xyz(lon_deg, lat_deg, gmst, r=1.003):
    """Geographic lon/lat to GCRS unit-sphere coordinates at sidereal angle gmst."""
    lam = np.radians(lon_deg) + gmst
    phi = np.radians(lat_deg)
    return r * np.cos(phi) * np.cos(lam), r * np.cos(phi) * np.sin(lam), r * np.sin(phi)


def _coastline_trace(gmst):
    xs, ys, zs = [], [], []
    for seg in COASTLINES:
        a = np.asarray(seg)
        x, y, z = _lonlat_to_xyz(a[:, 0], a[:, 1], gmst)
        xs += list(x) + [None]
        ys += list(y) + [None]
        zs += list(z) + [None]
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", line={"color": "#e8f0d8", "width": 2},
                        hoverinfo="skip", showlegend=False)


def earth_figure(p, t_gps):
    n = em.source_unit_vector(p["Alpha"], p["Delta"])
    pos_e, _ = em.earth_posvel(np.array([t_gps]))
    pos_e = pos_e[0]
    gmst = em.gmst_rad(np.array([t_gps]))[0]
    sites = {}
    for name in em.DETECTOR_INDEX:
        pos_s, _ = em.site_posvel(np.array([t_gps]), em.detector(name))
        sites[name] = 1.02 * pos_s[0] / np.linalg.norm(pos_s[0])
    sun_dir = -pos_e / np.linalg.norm(pos_e)

    traces = []
    x, y, z = _sphere(1.0)
    traces.append(go.Surface(
        x=x, y=y, z=z, surfacecolor=z, colorscale=[[0, "#0b3d91"], [1, "#2d7dd2"]],
        showscale=False, opacity=1.0, hoverinfo="skip", lighting={"ambient": 0.7},
    ))
    traces.append(_coastline_trace(gmst))
    # Equator and Greenwich meridian so the spin with sidereal time is visible
    lat_line = np.linspace(-90, 90, 60)
    gx, gy, gz = _lonlat_to_xyz(np.zeros_like(lat_line), lat_line, gmst)
    traces.append(go.Scatter3d(
        x=gx, y=gy, z=gz, mode="lines",
        line={"color": "#ffcc00", "width": 3}, hoverinfo="skip", name="Greenwich meridian",
    ))
    lon_line = np.linspace(0, 360, 90)
    ex, ey, ez = _lonlat_to_xyz(lon_line, np.zeros_like(lon_line), gmst)
    traces.append(go.Scatter3d(
        x=ex, y=ey, z=ez, mode="lines",
        line={"color": "#ffffff", "width": 2}, hoverinfo="skip", name="Equator",
    ))
    others = [k for k in sites if k != p["detector"]]
    traces.append(go.Scatter3d(
        x=[sites[k][0] for k in others], y=[sites[k][1] for k in others],
        z=[sites[k][2] for k in others], mode="markers+text",
        marker={"size": 5, "color": "#000000", "line": {"color": "#fff", "width": 1}},
        text=others, textposition="top center", textfont={"color": "#000000", "size": 11},
        name="Other detectors", hoverinfo="text",
    ))
    sel = sites[p["detector"]]
    traces.append(go.Scatter3d(
        x=[sel[0]], y=[sel[1]], z=[sel[2]], mode="markers+text",
        marker={"size": 8, "color": "#ff2d2d", "line": {"color": "#fff", "width": 1}},
        text=[p["detector"]], textposition="top center",
        textfont={"color": "#ff2d2d", "size": 13}, name=f"{p['detector']} (selected)",
    ))
    traces += _arrow(1.9 * n, 1.15 * n, "#ff4444", "CW propagation (from source)")
    traces += _arrow([0, 0, 0], 1.5 * sun_dir, "#ff9900", "To Sun", width=4)
    traces.append(go.Scatter3d(
        x=[0, 1.4 * n[0]], y=[0, 1.4 * n[1]], z=[0, 1.4 * n[2]], mode="lines",
        line={"color": "#ff4444", "width": 2, "dash": "dot"}, showlegend=False, hoverinfo="skip",
    ))

    fig = go.Figure(traces)
    date = Time(t_gps, format="gps").utc.iso[:19]
    fig.update_layout(**_scene("earth", 1.5),
                      title={"text": f"Earth view (GCRS axes)   GPS {t_gps:.0f}   {date} UTC",
                             "font": {"size": 13}})
    return fig


def orbit_figure(p, t_gps, ring):
    n = em.source_unit_vector(p["Alpha"], p["Delta"])
    pos_e, _ = em.earth_posvel(np.array([t_gps]))
    pos_e = pos_e[0]
    t_path = np.linspace(p["tstart"], max(t_gps, p["tstart"] + 1), 120)
    path, _ = em.earth_posvel(t_path)

    traces = [
        go.Scatter3d(x=ring[:, 0], y=ring[:, 1], z=ring[:, 2], mode="lines",
                     line={"color": "#99a", "width": 2}, name="Orbit (1 yr)", hoverinfo="skip"),
        go.Scatter3d(x=path[:, 0], y=path[:, 1], z=path[:, 2], mode="lines",
                     line={"color": "#1aa34a", "width": 4}, name="Observed so far", hoverinfo="skip"),
        _solid_sphere([0, 0, 0], 0.24, "#ff9900", "Sun"),
        _solid_sphere(pos_e, 0.06, "#2d7dd2", "Earth"),
    ]
    traces += _arrow(1.6 * n, 1.15 * n, "#ff4444", "CW propagation (from source)")
    fig = go.Figure(traces)
    layout = _scene("orbit", 1.6)
    layout["scene"].update(zaxis={"visible": False, "range": [-1.2, 1.2]}, aspectmode="manual",
                           aspectratio={"x": 1, "y": 1, "z": 0.75},
                           camera={"eye": {"x": 1.0, "y": 1.0, "z": 0.6}})
    fig.update_layout(**layout, title={"text": "Orbit view (SSB, AU)", "font": {"size": 13}})
    return fig


_ring_cache = {}


def orbit_ring(tstart):
    key = round(tstart)
    if key not in _ring_cache:
        t = tstart + np.linspace(0, 365.25 * 86400, 200)
        _ring_cache[key], _ = em.earth_posvel(t)
    return _ring_cache[key]


@app.callback(
    Output("earth", "figure"), Output("orbit", "figure"),
    Output("spec", "figure", allow_duplicate=True),
    Input("tslider", "value"), Input("params", "data"),
    prevent_initial_call=True,
)
def render_3d(t_days, p):
    if p is None or t_days is None:
        return (no_update,) * 3
    t_gps = p["tstart"] + t_days * 86400.0
    patch = Patch()
    patch["layout"]["shapes"][0]["x0"] = t_days
    patch["layout"]["shapes"][0]["x1"] = t_days
    return earth_figure(p, t_gps), orbit_figure(p, t_gps, orbit_ring(p["tstart"])), patch


if __name__ == "__main__":
    app.run(debug=True, port=8050)
