#! /usr/bin/env python
# -*- coding: utf-8 -*-


"""Code to handle pescoid wetting/dewetting visualization."""


from typing import List, Union

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import pandas as pd
from scipy.signal import find_peaks  # type: ignore
from scipy.signal import savgol_filter  # type: ignore
import seaborn as sns  # type: ignore


def _lineplot(xlabel_text: str, ylabel_text: str, title_text: str) -> None:
    """Plot a line plot with the provided labels and title."""
    plt.xlabel(xlabel_text)
    plt.ylabel(ylabel_text)
    plt.title(title_text)
    plt.margins(x=0)
    plt.show()
    plt.clf()


class PescoidVisualizer:
    """Class to visualize edge velocities from Z-stack data

    Initialize the class using file paths to the velocities and area data
    without preformatting.

    Attributes (public):
        velocities: DataFrame containing edge velocities
        area: DataFrame containing area data

    if smoothing is applied via smooth_area_curve, the following attributes are
    available:
        smoothed_area: Series containing a smoothed version of the area

    if phase_area_curve is applied, the following attributes are available:
        max_area_idx: int, index of the change point
        max_area: float, maximum area of the pescoid (from raw data)

    Methods
    ----------
    plot_area:
        Plots the area of the pescoid over time (frames)
    plot_single_line_velocity:
        Make a single line plot of edge velocities at a given index
    plot_multiple_line_velocities:
        Make a larger plot with subplots for each index in the provided list

    Examples:
    --------
    Initialize the class
    >>> dewetvizobj = PescoidVisualizer(
        velocities_file="velocities.mat", area_file="area.csv"
    )

    Apply smoothing via savgol filter
    >>> dewetvizobj.smooth_area_curve()

    Apply change point analysis to phase wetting/dewetting
    >>> dewetvizobj.phase_area_curve()

    Get the value and index where the area of pescoid is at its maximum
    >>> dewetvizobj.max_area
    >>> dewetvizobj.max_area_idx

    Plot the pescoid's raw area over time
    >>> dewetvizobj.plot_area()

    Plot the pescoid's smoothed area over time with a vertical line to indicate
    wetting/dewetting phase
    >>> dewetvizobj.plot_area(smoothed=True, phased=True)

    Plot the rate of the change in pescoid area over time
    >>> dewetvizobj.plot_rate_of_change(smoothed=True)
    >>> dewetvizobj.plot_rate_of_change(smoothed=False)

    Plot a single line of velocities according to index
    >>> dewetvizobj.plot_single_line_velocity(index=0)

    Multi-plot velocities for a list of indexes
    >>> dewetvizobj.plot_multiple_line_velocities(indexes=[0, 1, 2])
    """

    def __init__(
        self, area_file: str, velocities_file: Union[str, None] = None
    ) -> None:
        """Initialize the class"""
        # load data into dataframes
        self.area = pd.read_csv(area_file, header=0, skiprows=[1, 2, 3])[
            ["FRAME", "AREA"]
        ]
        if velocities_file:
            self.velocities: Union[pd.DataFrame, None] = pd.read_csv(
                velocities_file, delimiter="\t", header=None, usecols=list(range(270))
            )
        else:
            self.velocities = None

        # set plotting params
        self._scatter_size = 1
        self._set_matplotlib_publication_parameters()

    def _set_matplotlib_publication_parameters(self) -> None:
        """Set matplotlib parameters for publication quality plots."""
        font_size = 7
        plt.rcParams.update(
            {
                "font.family": "sans-serif",
                "font.sans-serif": "Arial",
                "font.size": font_size,
                "axes.titlesize": font_size,
                "axes.labelsize": font_size,
                "xtick.labelsize": font_size,
                "ytick.labelsize": font_size,
                "legend.fontsize": font_size,
            }
        )

    def _get_baseline_size(self) -> float:
        """Get the baseline size of the pescoid"""
        return self.area["AREA"].iloc[0]

    def _find_change_point(self, area_data: pd.Series) -> int:
        """Find the change point to indicate wetting/dewetting phase. Assumes a
        smooth enough curve that monotonicity is preserved."""
        peaks, _ = find_peaks(area_data)
        max_peak = peaks[np.argmax(area_data[peaks])]

        # Use the gradient to find when the slope starts to decrease
        gradient = np.gradient(area_data[:max_peak])
        decreasing_points = np.where(gradient < 0)[0]
        return decreasing_points[0] if decreasing_points.size > 0 else max_peak

    def _avg_rate_of_change(
        self,
        area_data: pd.Series,
        change_point_idx: int,
        phase: str,
    ) -> float:
        """Calculate the average rate of change of the data"""
        if phase == "wetting":
            return (area_data[change_point_idx] - area_data[0]) / change_point_idx
        elif phase == "dewetting":
            points_after = len(area_data) - change_point_idx
            return (area_data[-1] - area_data[change_point_idx]) / points_after
        else:
            raise ValueError("Invalid phase, must be 'wetting' or 'dewetting'")

    def _plot_phase_annotations(self, area_data: pd.Series):
        """Plot phase annotations: vertical line at change point and text with
        average slopes"""
        # plot vertical line at change point
        try:
            change_point = self.max_area_idx
        except AttributeError as e:
            raise AttributeError(
                "Change point not found, run phase_area_curve first"
            ) from e

        plt.axvline(x=change_point, color="red", linestyle="dashed")
        slope_wetting = self._avg_rate_of_change(area_data, change_point, "wetting")
        slope_dewetting = self._avg_rate_of_change(area_data, change_point, "dewetting")

        # set text annotation positions
        text_x_position = plt.xlim()[0] + (plt.xlim()[1] - plt.xlim()[0]) * 0.05
        text_offset = (plt.ylim()[1] - plt.ylim()[0]) * 0.03
        text_y_position = plt.ylim()[0] + text_offset

        plt.text(
            x=text_x_position,
            y=text_y_position,
            s=f"Slope dewetting: {slope_dewetting:.2f}",
            ha="left",
            va="bottom",
        )
        plt.text(
            x=text_x_position,
            y=text_y_position + text_offset,
            s=f"Slope wetting: {slope_wetting:.2f}",
            ha="left",
            va="bottom",
        )

    def plot_rate_of_change(self, smoothed: bool = False) -> None:
        """Plot the rate of change of the area curve"""
        title_text = "Pescoid area rate of change"
        try:
            if smoothed:
                rate_of_change = np.gradient(self.smoothed_area)
                title_text += " (smoothed)"
            else:
                rate_of_change = np.gradient(self.area["AREA"])
        except AttributeError:
            if smoothed:
                print("Smoothing the area curve before finding the rate of change")
                self.smooth_area_curve()
                rate_of_change = np.gradient(self.smoothed_area)
                title_text += " (smoothed)"
            else:
                print("Change point not found, run phase_area_curve first")
                return

        sns.lineplot(data=rate_of_change)

        try:
            plt.axvline(x=self.max_area_idx, color="red", linestyle="dashed")
        except AttributeError:
            print("Change point not found, run phase_area_curve first")

        _lineplot(
            xlabel_text="Time (frames)",
            ylabel_text="Rate of change (nm^2/s)",
            title_text=title_text,
        )

    def smooth_area_curve(self) -> None:
        """Convolutionally smooth the area curve via the Savitzky-Golay
        filter."""
        self.smoothed_area = savgol_filter(
            self.area["AREA"], window_length=15, polyorder=3
        )

    def phase_area_curve(self) -> None:
        """Add the change point to the class attributes"""
        # Apply smoothing if not already done
        if not hasattr(self, "smoothed_area"):
            print("Smoothing the area curve before finding the change point")
            self.smooth_area_curve()

        self.max_area_idx = self._find_change_point(self.smoothed_area)
        self.max_area = self.area["AREA"][self.max_area_idx]

    def plot_area(
        self,
        baseline: bool = True,
        smoothed: bool = False,
        phased: bool = False,
    ) -> None:
        """Plot the area of the pescoid over time

        Arguments:
            baseline: bool, whether to plot a horizontal line at the start as
            reference size
            smoothed: bool, whether to plot the smoothed data
            phased: bool, whether to indicate wetting/dewetting
            phase with a change point analysis
        """
        title_text = "Pescoid size during wetting/dewetting" + (
            " (smoothed)" if smoothed else ""
        )
        area_data = self.smoothed_area if smoothed else self.area["AREA"]

        if phased:
            if not smoothed:
                raise ValueError("Smoothing must be enabled to find change point")
            self._plot_phase_annotations(area_data)

        sns.lineplot(data=area_data)

        if baseline:
            plt.axhline(y=self._get_baseline_size(), color="gray", linestyle="dashed")

        _lineplot(
            xlabel_text="Time (frames)",
            ylabel_text="Area (nm^2)",
            title_text=title_text,
        )

    def plot_single_line_velocity(self, index: int) -> None:
        """Make a line plot of velocities at a given index"""
        if self.velocities is None:
            raise ValueError("No velocity data provided")
        data = self.velocities.iloc[index]
        fix, ax = plt.subplots()

        # plot positive and negative values in different colors
        positive_indices = data.index[data >= 0]
        negative_indices = data.index[data < 0]
        positive_values = data[data >= 0]
        negative_values = data[data < 0]

        ax.scatter(
            positive_indices,
            positive_values,
            color="blue",
            label="Positive",
            s=self._scatter_size,
        )
        ax.scatter(
            negative_indices,
            negative_values,
            color="red",
            label="Negative",
            s=self._scatter_size,
        )

        # plot dashed line through values
        ax.plot(data.index, data, color="gray", linestyle="dashed", linewidth=1)

        plt.axhline(y=0, color="black", linewidth=0.8)
        plt.margins(x=0)

        plt.xlabel("Time (frames)")
        plt.ylabel("Velocity (nm/s)")
        plt.title(f"Edge velocity at index {index}")
        plt.show()
        plt.clf()

    def plot_multiple_line_velocities(self, indexes: List[int]) -> None:
        """Create a larger plot with subplots for each index in the provided list.

        Arguments:
            indexes: List of indexes to plot
        """
        if self.velocities is None:
            raise ValueError("No velocity data provided")
        num_plots = len(indexes)
        fig, axs = plt.subplots(
            num_plots, 1, constrained_layout=True, sharex=True, sharey=True
        )

        for i, index in enumerate(indexes):
            data = self.velocities.iloc[index]
            positive = data[data >= 0]
            negative = data[data < 0]

            ax = axs[i] if num_plots > 1 else axs
            ax.scatter(
                positive.index,
                positive.values,
                color="blue",
                label="Positive",
                s=self._scatter_size,
            )
            ax.scatter(
                negative.index,
                negative.values,
                color="red",
                label="Negative",
                s=self._scatter_size,
            )

            # gray dashed line through values
            ax.plot(
                data.index, data.values, color="gray", linestyle="dashed", linewidth=1
            )

            # set labels
            ax.axhline(y=0, color="black", linewidth=0.8)
            ax.set_title(f"Index {index}")
            ax.margins(x=0)

        # Adjust layout, show legend, and show the plot
        fig.supxlabel("Time (frames)")
        fig.supylabel("Velocity (nm/s)")
        fig.suptitle("Edge velocities at multiple indexes")
        plt.show()
        plt.clf()
