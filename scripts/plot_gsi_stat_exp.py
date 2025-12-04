#!/usr/bin/env python3
"""
Plot O-F statistics from GSI stat files for multiple experiments.

Modernized for Python 3.11.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from matplotlib import rcParams, ticker, pyplot as plt
from matplotlib import gridspec as gspec

from pyGSI.gsi_stat import GSIstat


# ---- Configurable constants ----
IT_QUERY = "it == 1"
PLOT_TYPE = "mean"
OB_TYPES = [120, 220]
OB_VARS = ["t", "uv", "q"]
LEVELS = np.array([1000, 900, 800, 600, 400, 300, 250, 200, 150, 100, 50, 0])


@dataclass(frozen=True)
class DateRange:
    """Represents a start/end date range and produces 6-hourly cycle strings."""
    start: datetime
    end: datetime

    def cycles_6h(self) -> Tuple[List[str], List[str]]:
        """Return lists of expected gsistat filenames and cycle timestamps."""
        statfiles, cycles = [], []
        cur = self.start
        while cur <= self.end:
            cyc = cur.strftime("%Y%m%d%H")
            statfiles.append(f"gsistat.gdas.{cyc}")
            cycles.append(cyc)
            cur += timedelta(hours=6)
        return statfiles, cycles


def _validate_inputs(gsistat_dirs: List[Path], labels: List[str]) -> None:
    """Ensure each experiment directory exists and labels align with them."""
    if len(gsistat_dirs) != len(labels):
        raise ValueError(
            f"Number of --gsistats ({len(gsistat_dirs)}) must match number of --label ({len(labels)})."
        )
    missing = [p for p in gsistat_dirs if not p.is_dir()]
    if missing:
        raise FileNotFoundError(f"The following directories do not exist: {', '.join(map(str, missing))}")


def _set_matplotlib_defaults() -> None:
    """Apply consistent plot aesthetics for all figures."""
    rcParams["figure.subplot.left"] = 0.1
    rcParams["figure.subplot.top"] = 0.85
    rcParams["legend.fontsize"] = 12
    rcParams["axes.grid"] = True


def gen_figure(
    datadict: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    datatypestr: str,
    stattype: str,
    labels: List[str],
    sdate: datetime,
    edate: datetime,
    save: bool,
    plotdir: Path,
) -> None:
    """
    Generate a 1x3 figure of t / uv / q vs pressure (log scale).

    Args:
        datadict: Nested dictionary of experiment, stat type, variable arrays.
        datatypestr: Label for plot title (e.g., "RMSE", "Bias").
        stattype: Which statistic aggregation to plot ("mean", "aggr", or "sum").
        labels: List of experiment identifiers for the legend.
        sdate: Start datetime for annotation.
        edate: End datetime for annotation.
        save: If True, saves plots as PDF/PNG instead of showing interactively.
        plotdir: Directory path where figures are saved if `save` is True.
    """
    _set_matplotlib_defaults()
    colors = ["k", "r", "g", "b", "m", "c", "y"]
    fig = plt.figure(figsize=(10, 8))
    plt.subplots_adjust(hspace=0.3)
    gs = gspec.GridSpec(1, 3)
    y_levels = LEVELS[:-1]

    for v, var in enumerate(OB_VARS):
        ax = plt.subplot(gs[v])
        xmin, xmax = np.inf, -np.inf

        for e, expid in enumerate(labels):
            profile = datadict[expid][stattype][var][:-1]
            c = colors[e % len(colors)]
            ax.plot(profile, y_levels, marker="o", color=c, mfc=c, mec=c, label=expid)

            # Track min/max for consistent x-limits
            valid = profile[:-1] if var == "q" else profile
            pmin, pmax = np.nanmin(valid), np.nanmax(valid)
            xmin, xmax = min(xmin, pmin), max(xmax, pmax)

        # Axis labels and titles
        if v == 0:
            ax.legend(loc=0, numpoints=1)
            ax.set_ylabel("pressure (hPa)")

        var_labels = {"uv": ("m/s", "Winds"), "t": ("K", "Temperature"), "q": ("%", "Relative Humidity")}
        var_unit, var_name = var_labels.get(var, ("", var))

        ax.set_title(var_name, fontsize=14)
        ax.set_yscale("log")
        ax.set_ylim(1020, 50)
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=np.arange(1, 10)))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        ax.set_xlabel("count" if stattype == "sum" else f"magnitude ({var_unit})")

        # Add horizontal padding to x-limits
        if np.isfinite(xmin) and np.isfinite(xmax) and xmax > xmin:
            pad = (xmax - xmin) * 0.1
            ax.set_xlim(xmin - pad, xmax + pad)

    # Annotate figure header
    sdatestr, edatestr = sdate.strftime("%Y%m%d%H"), edate.strftime("%Y%m%d%H")
    plt.figtext(0.5, 0.93, f"{datatypestr} O-F ({sdatestr}-{edatestr})", ha="center", fontsize=18)

    # Save or show
    fname = f"gsistat_uvtq_{datatypestr}"
    if save:
        plotdir.mkdir(parents=True, exist_ok=True)
        plt.savefig(plotdir / f"{fname}.pdf", bbox_inches="tight")
        plt.savefig(plotdir / f"{fname}.png", bbox_inches="tight", dpi=150)
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    """Main CLI entry point: read data, compute aggregates, and plot."""
    # 1. Parse command-line args
    parser = argparse.ArgumentParser(
        description="Plot comparison of GSI O-F statistics (RMSE, Bias, Count) for multiple experiments."
    )
    parser.add_argument("-d", "--gsistats", nargs="+", required=True,
                        help="Directories containing GSI stat files (one per experiment).")
    parser.add_argument("-l", "--label", nargs="+", required=False,
                        help="Labels for experiments (must match --gsistats order).")
    parser.add_argument("-f", "--save-figure", action="store_true", dest="save_figure",
                        help="Save figures as PNG/PDF instead of showing interactively.")
    parser.add_argument("-p", "--plotdir", default="./", help="Output directory for saved figures.")
    parser.add_argument("-s", "--start-date", required=True, metavar="YYYYMMDDHH",
                        help="Start date/time (e.g., 2025010100).", dest="start_date")
    parser.add_argument("-e", "--end-date", required=True, metavar="YYYYMMDDHH",
                        help="End date/time (e.g., 2025010700).", dest="end_date")
    args = parser.parse_args()

    # 2. Prepare paths and dates
    save_figure = args.save_figure
    plotdir = Path(args.plotdir).expanduser().resolve()
    sdate = datetime.strptime(args.start_date, "%Y%m%d%H")
    edate = datetime.strptime(args.end_date, "%Y%m%d%H")
    date_range = DateRange(sdate, edate)
    statfiles, cycles = date_range.cycles_6h()

    gsistat_dirs = [Path(p).expanduser().resolve() for p in args.gsistats]
    labels = args.label if args.label else [p.name for p in gsistat_dirs]
    _validate_inputs(gsistat_dirs, labels)

    # 3. Read gsistat data
    rmses: Dict[str, dict] = {}
    counts: Dict[str, dict] = {}
    biases: Dict[str, dict] = {}

    for exp, gsistats_dir in zip(labels, gsistat_dirs):
        rmses[exp], counts[exp], biases[exp] = {}, {}, {}
        for gsistat, cycle in zip(statfiles, cycles):
            rmses[exp][cycle], counts[exp][cycle], biases[exp][cycle] = {}, {}, {}
            inputfile = gsistats_dir / gsistat
            if not inputfile.is_file():
                raise FileNotFoundError(f"Unable to find {inputfile} for cycle {cycle}")

            # Read gsistat file using pyGSI
            gdas = GSIstat(str(inputfile), cycle)

            # Extract desired variables and QC subsets
            for var in OB_VARS:
                stat = gdas.extract(var).query(IT_QUERY)
                tmp = stat[stat.index.isin(OB_TYPES, level="typ")]
                tmp = tmp[tmp.index.isin(["asm"], level="use")]
                rmses[exp][cycle][var] = tmp[tmp.index.isin(["rms"], level="stat")]
                counts[exp][cycle][var] = tmp[tmp.index.isin(["count"], level="stat")]
                biases[exp][cycle][var] = tmp[tmp.index.isin(["bias"], level="stat")]

    # 4. Aggregate across cycles
    n_cycles = len(cycles)
    n_levels = len(LEVELS)
    for exp in labels:
        rmses[exp]["mean"], rmses[exp]["aggr"] = {}, {}
        biases[exp]["mean"], biases[exp]["aggr"] = {}, {}
        counts[exp]["sum"] = {}

        for var in OB_VARS:
            # Build matrices across all cycles
            rmse_mat = np.empty((n_cycles, n_levels), dtype=float)
            bias_mat = np.empty((n_cycles, n_levels), dtype=float)
            cnt_mat = np.empty((n_cycles, n_levels), dtype=float)

            for i, cyc in enumerate(cycles):
                rmse_mat[i, :] = rmses[exp][cyc][var].values[0]
                bias_mat[i, :] = biases[exp][cyc][var].values[0]
                cnt_mat[i, :] = counts[exp][cyc][var].values[0]

            # Simple means across time
            rmses[exp]["mean"][var] = np.nanmean(rmse_mat, axis=0)
            biases[exp]["mean"][var] = np.nanmean(bias_mat, axis=0)

            # Weighted aggregates across time (using counts as weights)
            sum_c = np.nansum(cnt_mat, axis=0)
            with np.errstate(invalid="ignore", divide="ignore"):
                w_rmse = np.sqrt(np.nansum(cnt_mat * (rmse_mat ** 2.0), axis=0) / sum_c)
                w_bias = np.nansum(cnt_mat * bias_mat, axis=0) / sum_c
            w_rmse[sum_c <= 0] = np.nan
            w_bias[sum_c <= 0] = np.nan

            rmses[exp]["aggr"][var] = w_rmse
            biases[exp]["aggr"][var] = w_bias
            counts[exp]["sum"][var] = np.nansum(cnt_mat, axis=0)

    # 5. Plot results
    gen_figure(rmses, "RMSE", PLOT_TYPE, labels, sdate, edate, save_figure, plotdir)
    gen_figure(biases, "Bias", PLOT_TYPE, labels, sdate, edate, save_figure, plotdir)
    gen_figure(counts, "Count", "sum", labels, sdate, edate, save_figure, plotdir)


if __name__ == "__main__":
    main()
