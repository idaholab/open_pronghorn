#!/usr/bin/env python3
"""
plot_residuals.py
=================
Parse residuals from a MOOSE/Otter log file and render them with the same
polished dark dashboard style used by plot_debug_baffle.py.

Run simulation with:
    ../otter-opt -i yourfile.i 2>&1 | tee yourfile.log

Run plotter with:
    python plot_residuals.py yourfile.log -o all_residuals.png
"""

import argparse
import os
import re
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
FLOAT_RE = r"[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?"

# --- Line-level patterns -----------------------------------------------
# Different Otter/MOOSE log versions format residual lines slightly differently:
#   v3-style: "Momentum equation: Component 1 3.48e-05"
#             bare turbulence lines: "TKE_system 0.000453414"
#   v4-style: "Momentum equation: Component 3 [32m0.0014143[39m Linear its: 9"
#             turbulence prefixed:  "Advected system: TKE_system [32m0.0038[39m Linear its: 10"
# The patterns below match both, after ANSI codes are stripped.

MOMENTUM_RE = re.compile(rf"^\s*Momentum equation:\s*Component\s*(\d+)\s+({FLOAT_RE})\b")
PRESSURE_RE = re.compile(rf"^\s*Pressure equation:\s+({FLOAT_RE})\b")
ADVECTED_RE = re.compile(rf"^\s*Advected system:\s*(\S+)\s+({FLOAT_RE})\b")
SOLID_RE = re.compile(rf"^\s*(?:Currently Executing\s+)?Solid energy equation:\s+({FLOAT_RE})\b")
# Bare turbulence-quantity lines with no "Advected system:" prefix, e.g.
# "TKE_system 0.000453414" / "TKED_system 9.12324e-05" (v3-style logs).
BARE_SYSTEM_RE = re.compile(rf"^\s*(\w+_system)\s+({FLOAT_RE})\s*$")

MOMENTUM_LABELS = {1: "u momentum", 2: "v momentum", 3: "w momentum"}

ITER_RE = re.compile(r"^\s*Iteration\s+(\d+)\b")

# Keys that hint a line is residual-like but wasn't matched by any pattern above;
# used only for the --report-unparsed diagnostic.
RESIDUAL_HINT_KEYS = [
    "Momentum equation:",
    "Pressure equation:",
    "Advected system:",
    "Solid energy equation:",
]


def _momentum_label(component):
    return MOMENTUM_LABELS.get(component, f"momentum component {component}")


def _prettify_system_name(name):
    """Turn 'TKE_system' -> 'TKE', 'enthalpy_system' -> 'enthalpy', etc."""
    if name.lower().endswith("_system"):
        return name[: -len("_system")]
    return name

# Palette copied from plot_debug_baffle.py for visual consistency.
RESIDUAL_COLORS = [
    "#E74C3C",
    "#3498DB",
    "#2ECC71",
    "#F39C12",
    "#9B59B6",
    "#1ABC9C",
    "#E67E22",
    "#34495E",
]

FIGURE_BG = "#0F1117"
AXES_BG = "#1A1D27"
GRID_COLOR = "#222233"
SPINE_COLOR = "#333344"
TEXT_COLOR = "white"
MUTED_TEXT_COLOR = "#AAAAAA"
TICK_COLOR = "#888888"


def strip_ansi(text):
    return ANSI_RE.sub("", text)


def parse_log(log_path):
    series = defaultdict(lambda: {"iters": [], "res": []})
    current_it = None
    unparsed = []

    def record(name, value, lineno, line):
        if current_it is None:
            unparsed.append((lineno, "Residual found before first iteration", line))
        else:
            series[name]["iters"].append(current_it)
            series[name]["res"].append(value)

    with open(log_path, "r", errors="replace") as f:
        for lineno, raw_line in enumerate(f, start=1):
            line = strip_ansi(raw_line).strip()
            if not line:
                continue

            m_it = ITER_RE.match(line)
            if m_it:
                current_it = int(m_it.group(1))
                continue

            m = MOMENTUM_RE.match(line)
            if m:
                component = int(m.group(1))
                record(_momentum_label(component), float(m.group(2)), lineno, line)
                continue

            m = PRESSURE_RE.match(line)
            if m:
                record("pressure", float(m.group(1)), lineno, line)
                continue

            m = ADVECTED_RE.match(line)
            if m:
                record(_prettify_system_name(m.group(1)), float(m.group(2)), lineno, line)
                continue

            m = SOLID_RE.match(line)
            if m:
                record("solid energy", float(m.group(1)), lineno, line)
                continue

            m = BARE_SYSTEM_RE.match(line)
            if m:
                record(_prettify_system_name(m.group(1)), float(m.group(2)), lineno, line)
                continue

            if any(key in line for key in RESIDUAL_HINT_KEYS):
                unparsed.append((lineno, "Unmatched residual-like line", line))

    return series, unparsed


def apply_plot_style():
    """Apply the global dark/monospace style from plot_debug_baffle.py."""
    matplotlib.rcParams.update(
        {
            "axes.prop_cycle": matplotlib.cycler(color=RESIDUAL_COLORS),
            "font.family": "monospace",
            "figure.facecolor": FIGURE_BG,
            "savefig.facecolor": FIGURE_BG,
            "axes.facecolor": AXES_BG,
            "axes.edgecolor": SPINE_COLOR,
            "axes.labelcolor": MUTED_TEXT_COLOR,
            "xtick.color": TICK_COLOR,
            "ytick.color": TICK_COLOR,
            "text.color": TEXT_COLOR,
        }
    )


def style_ax(ax, title, xlabel, ylabel):
    """Style one axes to match the baffle debug dashboard."""
    ax.set_facecolor(AXES_BG)
    ax.set_title(title, color=TEXT_COLOR, fontsize=12, pad=8, fontfamily="monospace")
    ax.set_xlabel(xlabel, color=MUTED_TEXT_COLOR, fontsize=10)
    ax.set_ylabel(ylabel, color=MUTED_TEXT_COLOR, fontsize=10)
    ax.tick_params(colors=TICK_COLOR)

    for spine in ax.spines.values():
        spine.set_edgecolor(SPINE_COLOR)

    ax.grid(True, color=GRID_COLOR, linewidth=0.5, linestyle="--")
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())


def plot_series(series, output_file="residuals.png", show=True):
    apply_plot_style()

    fig, ax = plt.subplots(figsize=(12, 7), facecolor=FIGURE_BG)
    fig.suptitle(
        "MOOSE/Otter Residual Convergence",
        fontsize=16,
        color=TEXT_COLOR,
        y=0.98,
        fontweight="bold",
        fontfamily="monospace",
    )
    style_ax(ax, "Residuals vs SIMPLE Iteration", "Iteration", "Residual (log scale)")

    plotted_any = False
    for name, data in series.items():
        if data["iters"] and data["res"]:
            ax.semilogy(
                data["iters"],
                data["res"],
                linewidth=2.0,
                marker=".",
                markersize=4,
                label=name,
            )
            plotted_any = True

    if not plotted_any:
        raise RuntimeError("No residual data found to plot.")

    legend = ax.legend(
        loc="best",
        fontsize=9,
        framealpha=0.3,
        facecolor=AXES_BG,
        labelcolor=TEXT_COLOR,
        edgecolor="#444455",
        title="Equation",
        title_fontsize=9,
    )
    legend.get_title().set_color(TEXT_COLOR)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(output_file, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())

    if show:
        plt.show()
    else:
        plt.close(fig)


def default_output_name(log_path):
    """Derive 'logname_residuals.png' from 'path/to/logname.log'."""
    stem = os.path.splitext(os.path.basename(log_path))[0]
    return f"{stem}_residuals.png"


def main():
    parser = argparse.ArgumentParser(description="Plot residuals from a MOOSE/Otter log file.")
    parser.add_argument("logfile", nargs="?", default="step4s.log", help="Path to log file")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output image filename (default: <logname>_residuals.png)",
    )
    parser.add_argument("--show", action="store_true", help="Display the plot window")
    parser.add_argument(
        "--report-unparsed",
        action="store_true",
        help="Print lines that look like residuals but were not parsed",
    )
    args = parser.parse_args()

    output_file = args.output or default_output_name(args.logfile)

    series, unparsed = parse_log(args.logfile)

    print("Parsed residual counts:")
    for name, data in series.items():
        print(f"  {name:>16}: {len(data['res'])}")

    if args.report_unparsed and unparsed:
        print("\nUnparsed residual-like lines:")
        for lineno, reason, line in unparsed:
            print(f"  line {lineno}: {reason}: {line}")

    plot_series(series, output_file, show=args.show)
    print(f"\nSaved plot to {output_file}")


if __name__ == "__main__":
    main()
