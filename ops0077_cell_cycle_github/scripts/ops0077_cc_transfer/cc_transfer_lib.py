"""Shared code to freeze the NTC cell-cycle phase score and transfer it to other cells (e.g. gene knockouts).

The NTC score (ops0077_cell_cycle_score.py) is a closed principal curve through four axes per cell:
  phos = log(pRb) - log(Rb), p21, pRb (nuclear intensities, 4i) and relative DNA content (DAPI),
each background-subtracted / log-transformed / robust-z-scored PER WELL.

For transfer, EVERYTHING that was estimated from data is frozen from the NTC cells: per-well background, floor,
median and MAD, the curve nodes and the phase anchoring. New cells are normalised with those NTC constants, never
with their own statistics (renormalising a knockout population on itself would remove the very shift being measured),
and are projected onto the frozen curve without refitting.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

MARKERS = ["pRb", "Rb", "p21"]
AXES = ["phos", "p21", "pRb", "dna"]
WEIGHTS = np.array([1.0, 1.0, 1.0, 1.5])
DNA_COL, P21_COL = 3, 1
MIN_NUC_AREA_PX = 50
Z_CLIP = 4.0
SECTORS = ["G1", "S", "G2/M", "M-to-G1 exit"]
GRID_N = 720


# ---------------------------------------------------------------- circular helpers
def circ_diff(a, b):
    """Signed circular difference a - b in turns, in [-0.5, 0.5)."""
    return (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5


def to_phase(sincos):
    """Angle of a (sin, cos) pair, in turns [0, 1)."""
    return (np.arctan2(sincos[:, 0], sincos[:, 1]) / (2 * np.pi)) % 1.0


def mean_vector(phase):
    z = np.exp(2j * np.pi * np.asarray(phase)[~np.isnan(phase)]).mean()
    return z


# ---------------------------------------------------------------- closed principal curve
def project_closed(Xw, nodes):
    """Project onto the closed polyline through nodes; returns raw phase in [0,1) and distance."""
    P = np.vstack([nodes, nodes[:1]])
    seg = np.diff(P, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    cum = np.r_[0, np.cumsum(seglen)]
    best_d = np.full(len(Xw), np.inf)
    best_t = np.zeros(len(Xw))
    for i in range(len(seg)):
        s = seg[i]
        u = np.clip(((Xw - P[i]) @ s) / (s @ s + 1e-12), 0, 1)
        d = np.linalg.norm(Xw - (P[i] + u[:, None] * s), axis=1)
        m = d < best_d
        best_d[m], best_t[m] = d[m], cum[i] + u[m] * seglen[i]
    return (best_t / cum[-1]) % 1.0, best_d


def node_positions(nodes):
    P = np.vstack([nodes, nodes[:1]])
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    return cum[:-1] / cum[-1]


def fit_closed_curve(Xw, init_phase, n_nodes=16, n_iter=30, smooth=3):
    phase = init_phase.copy()
    for _ in range(n_iter):
        rank = np.argsort(np.argsort(phase, kind="stable"), kind="stable")
        idx = np.minimum(rank * n_nodes // len(phase), n_nodes - 1)  # equal-count bins along the phase
        nodes = np.array([Xw[idx == k].mean(0) for k in range(n_nodes)])
        for _ in range(smooth):  # periodic smoothing
            nodes = 0.25 * np.roll(nodes, 1, 0) + 0.5 * nodes + 0.25 * np.roll(nodes, -1, 0)
        phase, dist = project_closed(Xw, nodes)
    return nodes, phase, dist


def anchor_params(nodes, w):
    """Phase 0 = node with the lowest DNA content; forward direction passes the p21 minimum (S) before the DNA maximum."""
    pos = node_positions(nodes)
    dna_n, p21_n = nodes[:, DNA_COL] / w[DNA_COL], nodes[:, P21_COL] / w[P21_COL]
    p0, pmax, pp21 = pos[np.argmin(dna_n)], pos[np.argmax(dna_n)], pos[np.argmin(p21_n)]
    flip = not (((pp21 - p0) % 1.0) < ((pmax - p0) % 1.0))
    return float(p0), bool(flip)


def anchor_phase(raw_phase, p0, flip):
    out = (np.asarray(raw_phase) - p0) % 1.0
    return (-out) % 1.0 if flip else out


def init_phase_from_pca(Xw, rng):
    Xc = Xw - Xw.mean(0)
    _, _, vt = np.linalg.svd(Xc[rng.choice(len(Xc), min(len(Xc), 20000), replace=False)], full_matrices=False)
    p = Xc @ vt[:2].T
    return (np.arctan2(p[:, 1], p[:, 0]) / (2 * np.pi)) % 1.0


def fit_curve_on_ntc(X, w, rng, fit_n=20000):
    """Same fit as ops0077_cell_cycle_score.py (same RNG call order). Returns nodes, anchored phase, distance, p0, flip."""
    Xw = X * w
    init = init_phase_from_pca(Xw, rng)
    sub = rng.choice(len(X), min(fit_n, len(X)), replace=False)
    nodes, _, _ = fit_closed_curve(Xw[sub], init[sub])
    raw, dist = project_closed(Xw, nodes)
    p0, flip = anchor_params(nodes, w)
    return nodes, anchor_phase(raw, p0, flip), dist, p0, flip


# ---------------------------------------------------------------- 4i nuclear intensities and frozen per-well normalisation
def nuclear_table(get_feature, wells, index):
    """Nuclear mean intensity per marker from the CellProfiler cell / cytoplasm integrated and mean intensities."""
    df = pd.DataFrame(index=index)
    df["well"] = np.asarray(wells)
    for m in MARKERS:
        p = f"cp_single_object_{m}"
        cell_int, cell_mean = get_feature(f"{p}_cell_Intensity_IntegratedIntensity"), get_feature(f"{p}_cell_Intensity_MeanIntensity")
        cyt_int, cyt_mean = get_feature(f"{p}_cytoplasm_Intensity_IntegratedIntensity"), get_feature(f"{p}_cytoplasm_Intensity_MeanIntensity")
        with np.errstate(divide="ignore", invalid="ignore"):
            nuc_area = cell_int / cell_mean - cyt_int / cyt_mean  # px
            df[f"{m}_nuc"] = (cell_int - cyt_int) / nuc_area
        if m == "pRb":
            df["nuc_area_px"] = nuc_area
    bad = ~(df["nuc_area_px"] > MIN_NUC_AREA_PX)
    df.loc[bad, [f"{m}_nuc" for m in MARKERS] + ["nuc_area_px"]] = np.nan
    return df


def _pct(x, q):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, q)) if len(x) else np.nan


def _mad(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x - np.median(x))) * 1.4826) if len(x) else np.nan


def fit_frozen_constants(df, dna):
    """Estimate per-well constants on the NTC cells. df has {m}_nuc columns and 'well'; dna is the reliable DNA axis."""
    const = {"bg": {}, "floor": {}, "median": {}, "mad": {}}
    wells = sorted(df["well"].dropna().unique()) + ["__global__"]
    sel = lambda w: np.ones(len(df), bool) if w == "__global__" else (df["well"] == w).values
    logi = {}
    for m in MARKERS:
        col = df[f"{m}_nuc"].values
        const["bg"][m], const["floor"][m] = {}, {}
        logi[m] = np.full(len(df), np.nan)
        for w in wells:
            s = sel(w)
            bg = _pct(col[s], 0.5)
            sig = col[s] - bg
            floor = _pct(sig[sig > 0], 1)
            const["bg"][m][w], const["floor"][m][w] = bg, floor
            if w != "__global__":
                logi[m][s] = np.log(np.clip(sig, floor, None))
    axes = {"phos": logi["pRb"] - logi["Rb"], "p21": logi["p21"], "pRb": logi["pRb"], "dna": np.asarray(dna, float)}
    for a, v in axes.items():
        const["median"][a], const["mad"][a] = {}, {}
        for w in wells:
            x = v[sel(w)] if w != "__global__" else v
            const["median"][a][w] = float(np.nanmedian(x))
            const["mad"][a][w] = _mad(x)
    return const


def frozen_axes(df, dna, const):
    """Apply the frozen NTC constants. Wells unseen in the NTC fall back to the pooled NTC constants (flagged)."""
    wells = df["well"].astype(str).values
    known = np.array([w in const["bg"]["pRb"] for w in wells])
    key = np.where(known, wells, "__global__")
    pick = lambda d, name: np.array([d[name][k] for k in key], float)
    logi = {}
    for m in MARKERS:
        sig = df[f"{m}_nuc"].values - pick({m: const["bg"][m]}, m)
        floor = pick({m: const["floor"][m]}, m)
        with np.errstate(invalid="ignore"):
            logi[m] = np.log(np.maximum(sig, floor))
        logi[m][np.isnan(df[f"{m}_nuc"].values)] = np.nan
    raw = {"phos": logi["pRb"] - logi["Rb"], "p21": logi["p21"], "pRb": logi["pRb"], "dna": np.asarray(dna, float)}
    z = pd.DataFrame(index=df.index)
    for a, v in raw.items():
        med, mad = pick({a: const["median"][a]}, a), pick({a: const["mad"][a]}, a)
        z[a] = np.clip((v - med) / mad, -Z_CLIP, Z_CLIP)
    return z, ~known


# ---------------------------------------------------------------- sectors defined from the NTC curve only
def curve_points(nodes, p0, flip, n=GRID_N):
    """Points along the anchored closed curve at n equally spaced anchored phases (weighted space)."""
    pos = node_positions(nodes)
    g = (np.arange(n) + 0.5) / n
    raw = (p0 - g) % 1.0 if flip else (p0 + g) % 1.0
    P = np.vstack([nodes, nodes[:1]])
    xp = np.r_[pos, 1.0]
    return g, np.column_stack([np.interp(raw, xp, P[:, j]) for j in range(P.shape[1])])


def define_sectors(nodes, p0, flip, w):
    """Label the loop by the curve's DNA content (from the NTC fit only): G1 low, S rising, G2/M high, then the exit descent."""
    g, pts = curve_points(nodes, p0, flip)
    dna = pts[:, DNA_COL] / w[DNA_COL]
    f = (dna - dna.min()) / (dna.max() - dna.min())
    peak = int(np.argmax(dna))
    labels = np.empty(len(g), object)
    rising = np.arange(len(g)) <= peak
    labels[f < 0.25] = "G1"
    labels[(f >= 0.25) & (f < 0.75) & rising] = "S"
    labels[f >= 0.75] = "G2/M"
    labels[(f >= 0.25) & (f < 0.75) & ~rising] = "M-to-G1 exit"
    return labels.tolist(), f.tolist()


def sector_of(phase, ref):
    labels = np.array(ref["sector_grid"]["labels"], object)
    out = np.full(len(phase), "", object)
    ok = ~np.isnan(phase)
    out[ok] = labels[np.minimum((phase[ok] * len(labels)).astype(int), len(labels) - 1)]
    return out


def sector_intervals(ref):
    """Human-readable phase intervals of each sector."""
    labels = ref["sector_grid"]["labels"]
    n, out, start = len(labels), [], 0
    for i in range(1, n + 1):
        if i == n or labels[i] != labels[start]:
            out.append((labels[start], start / n, i / n))
            start = i
    return out


# ---------------------------------------------------------------- reference I/O and projection
def save_reference(ref, path):
    json.dump(ref, open(path, "w"), indent=1)


def load_reference(path):
    ref = json.load(open(path))
    ref["curve_nodes"] = np.array(ref["curve_nodes"], float)
    ref["weights"] = np.array(ref["weights"], float)
    return ref


def project_to_reference(z, ref, axes=AXES):
    """Project cells with all four axes finite onto the frozen NTC curve. Returns anchored phase and distance (NaN if incomplete)."""
    Z = z[axes].values
    ok = np.isfinite(Z).all(axis=1)
    phase, dist = np.full(len(z), np.nan), np.full(len(z), np.nan)
    if ok.any():
        raw, d = project_closed(Z[ok] * ref["weights"], ref["curve_nodes"])
        phase[ok] = anchor_phase(raw, ref["anchor"]["p0"], ref["anchor"]["flip"])
        dist[ok] = d
    return phase, dist
