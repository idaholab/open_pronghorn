#!/usr/bin/env python3

"""Shared plotting utilities for the vortex-shedding benchmark quantities."""

from pathlib import Path

import matplotlib.pyplot as plt

BACKGROUND_COLOR = "#F7F8FA"
TEXT_COLOR = "#24292F"
MUTED_TEXT_COLOR = "#57606A"
TRACK_COLOR = "#D0D7DE"
RANGE_COLOR = "#0072B2"
PASS_COLOR = "#00875A"
FAIL_COLOR = "#C9472D"


def plot_benchmark_quantity(
    *,
    value,
    bounds,
    title,
    symbol,
    output_file,
    decimals=3,
):
    """Plot one computed value against its accepted literature interval."""

    lower_bound, upper_bound = bounds
    range_width = upper_bound - lower_bound
    if range_width <= 0:
        raise ValueError("The upper benchmark bound must exceed the lower bound")

    within_range = lower_bound <= value <= upper_bound
    status_color = PASS_COLOR if within_range else FAIL_COLOR
    status = "WITHIN LITERATURE RANGE" if within_range else "OUTSIDE LITERATURE RANGE"
    value_format = f"{{:.{decimals}f}}"

    plot_min = min(lower_bound, value) - 0.75 * range_width
    plot_max = max(upper_bound, value) + 0.75 * range_width

    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.titlecolor": TEXT_COLOR,
        }
    ):
        fig, ax = plt.subplots(figsize=(7.2, 3.1), facecolor=BACKGROUND_COLOR)
        ax.set_facecolor(BACKGROUND_COLOR)

        # A neutral track provides context if a computed value falls outside the
        # blue literature interval.
        ax.plot(
            [plot_min, plot_max],
            [0, 0],
            color=TRACK_COLOR,
            linewidth=4,
            solid_capstyle="round",
            zorder=1,
        )
        ax.plot(
            [lower_bound, upper_bound],
            [0, 0],
            color=RANGE_COLOR,
            linewidth=16,
            solid_capstyle="round",
            zorder=2,
        )
        ax.scatter(
            value,
            0,
            marker="D",
            s=150,
            color=status_color,
            edgecolor="white",
            linewidth=2,
            zorder=4,
        )

        ax.annotate(
            f"Current: {value_format.format(value)}",
            xy=(value, 0),
            xytext=(0, 27),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color=TEXT_COLOR,
            fontweight="bold",
        )
        for bound, alignment, label in (
            (lower_bound, "right", "Lower limit"),
            (upper_bound, "left", "Upper limit"),
        ):
            ax.annotate(
                f"{label}\n{value_format.format(bound)}",
                xy=(bound, 0),
                xytext=(((-7 if alignment == "right" else 7), -24)),
                textcoords="offset points",
                ha=alignment,
                va="top",
                color=MUTED_TEXT_COLOR,
                fontsize=9.5,
                linespacing=1.3,
            )

        ax.text(
            0,
            1.04,
            title,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            color=TEXT_COLOR,
            fontsize=17,
            fontweight="bold",
        )
        ax.text(
            0,
            0.95,
            rf"Schäfer–Turek benchmark  ·  ${symbol}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color=MUTED_TEXT_COLOR,
            fontsize=10,
        )
        ax.text(
            1,
            1.04,
            status,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            color=status_color,
            fontsize=9,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.42",
                "facecolor": "white",
                "edgecolor": status_color,
                "linewidth": 1.2,
            },
        )

        ax.set_xlim(plot_min, plot_max)
        ax.set_ylim(-0.62, 0.66)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        fig.tight_layout(pad=1.6)
        fig.savefig(
            Path(output_file),
            dpi=300,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
        )
        plt.close(fig)
