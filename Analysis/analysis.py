import os
import pandas as pd
import matplotlib.pyplot as plt
import ast
import re
import numpy as np
import scipy.stats as stats
from PAHMC import Input
import sys
from scipy.stats import gaussian_kde
from quantiphy import Quantity
import matplotlib.ticker as ticker

figsize: set = (5, 4)
dpi: int = 150


class Settings:
    def __init__(self, figsize=(5, 4), dpi=150):
        self.figsize = figsize
        self.dpi = dpi


class RRKM(Settings):
    def __init__(self, folder: str):
        self.folder = folder

        super().__init__(figsize, dpi)
        self.rates = {"D_hop": [], "H_hop": [], "D_diss": [], "H_diss": []}
        self.pairs = {"Hops": [], "Dissociations": []}

        for file in os.listdir(folder):
            if file.startswith("D"):
                D_file, H_file = file, f"H{file[1:]}"
                if os.path.exists(os.path.join(folder, H_file)):
                    if "to" in file:
                        self.pairs["Hops"].append([D_file, H_file])
                    elif "diss" in file:
                        self.pairs["Dissociations"].append([D_file, H_file])
                if "to" in file:
                    self.rates["D_hop"].append(file)
                elif "diss" in file:
                    self.rates["D_diss"].append(file)
            elif file.startswith("H"):
                if "to" in file:
                    self.rates["H_hop"].append(file)
                elif "diss" in file:
                    self.rates["H_diss"].append(file)

    def plot(
        self,
        subset: str,
        filter_function=None,
        title: str = "",
        ylim: list = [10e-20, 10e13],
        xlim: list = [0.0, 5.2],
    ):
        """
        Generate a plot for a specified subset of data files.

        Parameters:
        subset (str): The key to identify the subset of data files to plot.
        filter_function (callable, optional): A function to filter the filenames in the subset. Defaults to None.
        title (str, optional): The title of the plot. Defaults to an empty string.
        ylim (list, optional): A list specifying the y-axis limits. Defaults to an empty list.
        xlim (list, optional): A list specifying the x-axis limits. Defaults to [0.0, 5.2].

        Example Usage:
        analysis.plot(subset="H_diss", filter_function=lambda x: 'a' not in x, title="Hydrogen Dissociation", ylim=[0, 1e5], xlim=[0.0, 5.2])
        """

        plt.figure(figsize=self.figsize, dpi=self.dpi)

        subset = self.rates.get(subset)

        if filter_function is not None:
            subset = list(filter(filter_function, subset))

        colors = self._get_colormap(subset)

        for i, filename in enumerate(subset):
            filepath = os.path.join(self.folder, filename)
            data = np.loadtxt(filepath, skiprows=2)
            label = self._parse_label(filename)

            color, line_style = self._get_color_and_style(filename, colors)
            plt.plot(
                data[:, 0], data[:, 1], label=label, color=color, linestyle=line_style
            )

        if len(ylim) > 0:
            plt.ylim(ylim)

        if len(xlim) > 0:
            plt.xlim(xlim)
            plt.xticks(np.arange(xlim[0], xlim[1] + 0.4, 0.4))

        plt.title(title)
        self._plot_plot()

    def plot_pairs(
        self, subset: str, title: str = "", ylim: list = [], xlim: list = [0.0, 5.2]
    ):
        """
        Generate plots for pairs of data files.

        Parameters:
        subset (str): The key to identify the subset of pairs of data files to plot.
        title (str, optional): The title of the plot. Defaults to an empty string.
        ylim (list, optional): A list specifying the y-axis limits. Defaults to an empty list.
        xlim (list, optional): A list specifying the x-axis limits. Defaults to [0.0, 5.2].

        Example Usage:
        analysis.plot_pairs(subset="H_pairs", title="Hydrogen Pairs", ylim=[0, 1e5], xlim=[0.0, 5.2])
        """

        subset = self.pairs.get(subset)

        for pair in subset:
            plt.figure(figsize=self.figsize, dpi=self.dpi)
            for i, filename in enumerate(pair):
                filepath = os.path.join(self.folder, filename)
                data = np.loadtxt(filepath, skiprows=2)
                label = f"{self._parse_label(filename)}"

                color = "dimgray" if "H" in filename else "darkgray"

                plt.plot(data[:, 0], data[:, 1], label=label, color=color)

            if len(ylim) > 0:
                plt.ylim(ylim)

            if len(xlim) > 0:
                plt.xlim(xlim)
                plt.xticks(np.arange(xlim[0], xlim[1] + 0.4, 0.4))

            plt.title(title)
            self._plot_plot()

    def _plot_plot(self, sort_legend_by="pairs"):
        """
        Finalize the plot with labels, scales, and legend.

        Parameters:
        sort_legend_by (str): Criteria to sort the legend. Options are "alphabetical" or "pairs".
        """
        plt.xticks
        plt.yscale("log")
        plt.xlabel("Energy (eV)")
        plt.ylabel("Reaction rate (s$^{-1}$)")

        handles, labels = plt.gca().get_legend_handles_labels()
        if sort_legend_by == "alphabetical":
            labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        elif sort_legend_by == "pairs":
            labels, handles = zip(
                *sorted(zip(labels, handles), key=lambda t: self._extract_pair(t[0]))
            )

        plt.legend(handles, labels, loc="lower right")
        plt.show()
        plt.clf()

    def _get_colormap(self, subset):
        colors = {}
        color_cycle = plt.cm.Greys(np.linspace(0.3, 1.0, len(subset)))
        color_cycle = plt.cm.tab10.colors
        unique_numbers = sorted(
            set(
                int(num)
                for filename in subset
                for num in self._extract_pair(filename)
                if num.isdigit()
            )
        )

        for filename in subset:
            pair = self._extract_pair(filename)
            if "diss" in filename:
                number = (
                    int(pair[0])
                    if len(pair) == 1 and pair[0].isdigit()
                    else int(pair[1]) if len(pair) > 1 and pair[1].isdigit() else None
                )
                if number is not None:
                    index = unique_numbers.index(number)
                    colors[pair] = color_cycle[index % len(color_cycle)]
            elif pair not in colors:
                colors[pair] = color_cycle[len(colors) % len(color_cycle)]
        return colors

    def _get_color_and_style(self, filename, colors):
        pair = self._extract_pair(filename)
        color = colors.get(pair, "black")
        if "diss" in filename or len(pair) == 1:
            line_style = "-"
        else:
            numbers = [x for x in self._parse_label(filename)[2:].split(" to ")]
            line_style = "-" if numbers[0] > numbers[1] else "--"
        return color, line_style

    def _extract_pair(self, filename):
        base_name = os.path.splitext(filename)[0]
        pair = base_name[1:]
        pair = pair.replace("diss", "")  # Remove "diss" if present
        return tuple(sorted(pair.split("to")))

    def _parse_label(self, label):
        base_name = os.path.splitext(label)[0]
        if "to" in base_name:
            parts = base_name.split("to")
            return f"{base_name[0]} {parts[0][1:]} to {parts[1]}"
        elif "diss" in base_name:
            parts = base_name.split("diss")
            return f"{base_name[0]} {parts[0][1:]}"
        return label


class MonteCarlo(Settings):
    def __init__(
        self,
        folder: str,
        simulations: list,
        numbering: list,
        name: str,
        input: str,
        perdeuterated: bool = False,
        symmetries: dict = {},
    ):
        super().__init__(figsize, dpi)
        file_types = {
            "data": "_data.csv",
            "hops": "_key_hops.csv",
            "log": ".log",
            "output": ".out",
            "end_structures": "_end_structures.out",
        }
        self.files = {
            label: [os.path.join(folder, sim + file_type) for sim in simulations]
            for label, file_type in file_types.items()
        }
        self.perdeuterated = perdeuterated
        self.name = name
        self.numbering = self.process_numbering(numbering)
        self.data = self.process_data()
        self.symmetries = symmetries
        self.dissociation_counts, self.dissociation_positions = self.process_output()
        try:
            self.hops = self.process_hops()
        except FileNotFoundError:
            pass
        try:
            self.end_structures = self.process_endstruct()
        except FileNotFoundError:
            self.end_structures = None
        self.input = input
        self.position_times = self.hop_test()
        self.ratios = self.RRKM_rates()

    def hop_test(self):
        yaml = Input.Input_reader(self.input)
        from_list = list(set([x.split("to")[0] for x in self.hops.keys()]))
        results = pd.DataFrame()

        for origin in from_list:
            if isinstance(origin, list):
                symmetrical_origin = origin
            else:
                symmetrical_origin = [origin]
            symmetrical_origin.extend(self.symmetries.get(origin, []))
            filtered_hops = pd.DataFrame()
            for prefix in symmetrical_origin:
                filtered_columns = [
                    col for col in self.hops.columns if col.split("to")[0] == prefix
                ]
                filtered_hops = pd.concat(
                    [filtered_hops, self.hops[filtered_columns]], axis=1
                )

            filtered_hops = filtered_hops.reset_index(drop=True)
            summed_results = pd.Series(0, index=filtered_hops.index)

            for r_key in filtered_hops.columns:
                column_data = filtered_hops[r_key]

                E = yaml.energy
                rate_idx = np.argmin((np.abs(yaml.reactionrates[r_key][0, :] - E)))
                r_rate = yaml.reactionrates[r_key][1, rate_idx]
                dt = 1 / np.sum(r_rate)

                summed_results += column_data * dt

            results[origin] = summed_results

        summed_results = pd.DataFrame(index=results.index)

        for num in set(
            col[1:] for col in results.columns if col.startswith(("H", "D"))
        ):
            summed_results[f"{num}"] = results.get(f"H{num}", 0) + results.get(
                f"D{num}", 0
            )

        # check if summes_results contains symmetries and sum them together
        for origin in summed_results.columns:
            for key, value in self.symmetries.items():
                if isinstance(value, int):
                    value = [value]
                values = [str(x) for x in value]
                if origin in values:
                    summed_results[str(key)] += summed_results[str(origin)]
                    summed_results.drop(
                        columns=[origin], inplace=True
                    )  # Drop the column

        summed_results = summed_results.div(summed_results.sum(axis=1), axis=0)

        return summed_results

    def RRKM_rates(self):
        yaml = Input.Input_reader(self.input)

        rates_avg = {}
        rates_H = {}
        rates_D = {}

        for pos in self.symmetries.keys():
            rate_avg = 0
            rate_H = 0
            rate_D = 0
            for mol in ["H", "D"]:
                r_key = f"{mol}{pos}diss"
                E = yaml.energy
                rate_idx = np.argmin((np.abs(yaml.reactionrates[r_key][0, :] - E)))
                r_rate = yaml.reactionrates[r_key][1, rate_idx]
                if self.perdeuterated:
                    if mol == "D":
                        multiplier = 10 / 11
                    elif mol == "H":
                        multiplier = 1 / 11
                else:
                    if mol == "D":
                        multiplier = 1 / 11
                    elif mol == "H":
                        multiplier = 10 / 11
                rate_avg += r_rate * multiplier
                if mol == "H":
                    rate_H += r_rate
                elif mol == "D":
                    rate_D += r_rate
            rates_avg[pos] = rate_avg
            rates_H[pos] = rate_H
            rates_D[pos] = rate_D

        min_value_avg = min(rates_avg.values())
        min_value_H = min(rates_H.values())
        min_value_D = min(rates_D.values())

        ratios_avg = {key: value / min_value_avg for key, value in rates_avg.items()}
        ratios_H = {key: value / min_value_H for key, value in rates_H.items()}
        ratios_D = {key: value / min_value_D for key, value in rates_D.items()}

        ratios_avg = pd.DataFrame([ratios_avg], columns=ratios_avg.keys())
        ratios_H = pd.DataFrame([ratios_H], columns=ratios_H.keys())
        ratios_D = pd.DataFrame([ratios_D], columns=ratios_D.keys())

        ratios_avg = ratios_avg.div(ratios_avg.sum(axis=1), axis=0)
        ratios_H = ratios_H.div(ratios_H.sum(axis=1), axis=0)
        ratios_D = ratios_D.div(ratios_D.sum(axis=1), axis=0)

        ratios = pd.concat([ratios_avg, ratios_H, ratios_D], keys=["Avg", "H", "D"])

        return ratios

    def process_data(self):
        data_frames = []

        for file in self.files["data"]:
            df = pd.read_csv(file, na_values=["None"])
            df.columns = df.columns.str.strip()
            df.set_index(df.columns[0], inplace=True)
            df.loc[df["Diss pos"].isna(), "Diss time"] = np.nan
            df["# H hops"] = df["# hops"] - df["# D hops"]
            if not self.perdeuterated:
                df["# H hops"] /= 10
                df["HH time"] /= 9
            else:
                df["# D hops"] /= 10
                df["DD time"] /= 9
            # df.dropna(subset=['Diss pos'], inplace=True) # TODO: Discuss wheter this data needs to be used or not
            data_frames.append(df)

        return pd.concat(data_frames)

    def process_output(self):
        dissociation_counts = []
        dissociation_positions = []

        for file in self.files["output"]:
            with open(file, "r") as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("Dissociations for Energy"):
                        counts = lines[i + 1].strip().split(", ")
                        counts_dict = {}
                        for item in counts:
                            parts = item.split(": ")
                            if len(parts) == 2:
                                k, v = parts
                                counts_dict[k] = int(v)
                        dissociation_counts.append(counts_dict)
                    elif line.startswith("Dissociation positions for Energy"):
                        positions_dict = {}
                        for pos_line in lines[i + 1 : i + 11]:
                            parts = pos_line.strip().split(": ")
                            if len(parts) == 2:
                                k, v = parts
                                positions_dict[int(k)] = int(v)
                        dissociation_positions.append(positions_dict)

        df = pd.DataFrame(dissociation_positions)

        dissociation_positions_merged = pd.DataFrame()

        for key, value in self.symmetries.items():
            if isinstance(value, list):
                for v in value:
                    if v in df.columns and key in df.columns:
                        dissociation_positions_merged[key] = pd.concat(
                            [df[key], df[v]], ignore_index=True
                        )
            else:
                if value in df.columns and key in df.columns:
                    dissociation_positions_merged[key] = pd.concat(
                        [df[key], df[value]], ignore_index=True
                    )

        return pd.DataFrame(dissociation_counts), pd.DataFrame(
            dissociation_positions_merged
        )

    def process_hops(self):
        data_frames = []

        for file in self.files["hops"]:
            df = pd.read_csv(file)
            df.set_index(df.columns[0], inplace=True)
            data_frames.append(df)

        return pd.concat(data_frames)

    def process_endstruct(self):
        data_frames = []

        for file in self.files["end_structures"]:
            df = pd.read_csv(file, sep="\t", header=None, names=["ID", "Structures"])
            df.set_index("ID", inplace=True)
            df["Structures"] = df["Structures"].apply(self.safe_literal_eval)
            df["End position"] = None
            atom = "D" if not self.perdeuterated else "H"

            df["End position"] = df["Structures"].apply(
                lambda structures: self.find_end_position(structures, atom)
            )

            data_frames.append(df)

        df = pd.concat(data_frames)
        df.replace(-1, pd.NA, inplace=True)
        df.dropna(subset=["End position"], inplace=True)

        return df

    def safe_literal_eval(self, value):
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value.split(",")

    def find_end_position(self, structures: list, atom: str):
        for struct_list, num_list in zip(structures, self.numbering):
            for i, s in enumerate(struct_list):
                if s == atom:
                    return num_list[i]
        return -1

    def process_numbering(self, numbering: list):
        if "(" in numbering:
            groups = re.findall(r"\((.*?)\)", numbering)
            return [list(map(int, group.split(","))) for group in groups]
        else:
            return [list(map(int, group.split(","))) for group in numbering.split(";")]

    def all_plots(self):
        self.bar_hops("^D")
        self.bar_hops("^H")
        self.histogram_hops()
        self.histogram_hops_comparison()
        self.histogram_time()
        self.histogram_time_comparison()
        self.histogram_time_percentage()

    def plot_ratios(self):
        result = self.ratios  # TODO: Split ratios into H/D
        mean_positions = result.mean()
        std_positions = result.std()

        index = np.arange(len(mean_positions.index))

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index,
            mean_positions.values,
            yerr=std_positions.values,
            capsize=5,
            color="dimgray",
        )
        plt.xlabel("Dissociation position")
        plt.ylabel(
            f"Relative dissociation rate at {Input.Input_reader(self.input).energy} eV"
        )
        plt.xticks(index, mean_positions.index)
        plt.show()

    def histogram_position_times(self):
        results = self.position_times

        mean_positions = results.mean()
        std_positions = results.std()
        index = np.arange(len(mean_positions.index))

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index,
            mean_positions.values,
            yerr=std_positions.values,
            capsize=5,
            color="dimgray",
        )
        plt.xlabel("Position")
        plt.ylabel("Total time spent")
        plt.xticks(index, mean_positions.index)
        plt.show()

    def bar_hops(self, regex: str):
        df = self.hops.filter(regex=regex)
        df = df.div(df.sum(axis=1), axis=0) * 100

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        df.mean().plot(kind="bar", color="dimgray", yerr=df.std(), capsize=4)
        plt.ylabel("Occurences [%]")
        plt.xlabel("Hop")
        plt.xticks(rotation=60)
        plt.show()

    def histogram_time(self):
        data = self.data.dropna(subset=["Diss pos"])
        data = data["Diss time"]  # * 1000
        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.hist(data, bins=128, color="dimgray")
        plt.xlabel("Dissociation time [s]")
        plt.ylabel("Occurences")
        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        plt.show()

    def histogram_time_kde(self, filter="Diss time"):
        data = self.data.dropna(subset=["Diss pos"])
        data = data[filter]

        plt.figure(figsize=self.figsize, dpi=self.dpi)

        bin_width_ms = 1e-7

        bins = np.arange(0, data.max() + bin_width_ms, bin_width_ms)
        hist, bin_edges = np.histogram(data, bins=bins)

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.plot(bin_edges[:-1], hist, drawstyle="steps-post", color="darkgray")

        kde = gaussian_kde(data, bw_method="scott")
        x = np.linspace(0, data.max(), 1000)
        kde_values = kde(x)
        plt.plot(
            x, kde_values * len(data) * bin_width_ms, label="KDE fit", color="dimgray"
        )

        median_peak = np.median(data)

        half_max = np.max(kde_values) / 2
        fwhm_indices = np.where(kde_values >= half_max)[0]
        fwhm = x[fwhm_indices[-1]] - x[fwhm_indices[0]]

        plt.axvline(
            median_peak,
            color="r",
            linestyle="--",
            label=f"Median peak: {Quantity(median_peak, 's')}",
        )
        plt.axvline(
            x[fwhm_indices[0]],
            color="g",
            linestyle="--",
            label=f"FWHM: {Quantity(fwhm, 's')}",
        )
        plt.axvline(x[fwhm_indices[-1]], color="g", linestyle="--")

        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

        plt.xlabel("Time")
        plt.ylabel("Frequency")
        plt.title(f"Histogram of {filter}")
        plt.legend()
        plt.show()

    def histogram_hops(self):
        data = self.data.dropna(subset=["Diss pos"])
        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.hist(data["# hops"], bins=128, color="dimgray")
        plt.xlabel("Total number of scrambling hops")
        plt.ylabel("Occurences")
        plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        plt.show()

    def histogram_hops_comparison(self):
        data = self.data.dropna(subset=["Diss pos"])

        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)

        axs[0].hist(data["# D hops"], bins=128, label="D hops", color="dimgray")
        axs[1].hist(data["# H hops"], bins=128, label="H hops", color="darkgray")

        fig.supxlabel("Number of scrambling hops")
        fig.supylabel("Occurences")

        for ax in axs:
            ax.label_outer()

        plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

        axs[0].legend()
        axs[1].legend()
        plt.show()

    def histogram_time_percentage(self):
        if not self.perdeuterated:
            hh_time = (self.data["HH time"] * 9 / self.data["Diss time"]).dropna() * 100
            hd_time = (self.data["HD time"] / self.data["Diss time"]).dropna() * 100
        else:
            hd_time = (self.data["HD time"] / self.data["Diss time"]).dropna() * 100
            dd_time = (self.data["DD time"] * 9 / self.data["Diss time"]).dropna() * 100

        fig, (ax1, ax2) = plt.subplots(
            1, 2, sharey=True, figsize=self.figsize, dpi=self.dpi
        )
        fig.subplots_adjust(hspace=0)

        if not self.perdeuterated:
            ax1.hist(hd_time, bins=128, label="HD", color="dimgray")
            ax2.hist(hh_time, bins=128, label="HH", color="darkgray")
        else:
            ax1.hist(hd_time, bins=128, label="HD", color="dimgray")
            ax2.hist(dd_time, bins=128, label="DD", color="darkgray")

        ax1.set_xlim(5, 15)
        ax2.set_xlim(85, 95)

        ax1.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.tick_params(labelleft=False)
        ax2.yaxis.tick_right()

        fig.supxlabel("Percentage time spent [%]")
        fig.supylabel("Occurrences")

        d = 0.02
        kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False)

        ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        kwargs.update(transform=ax2.transAxes)
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)

        fig.legend(loc="upper right")
        plt.show()

    def histogram_time_comparison(self):
        data = self.data.dropna(subset=["Diss pos"])
        if not self.perdeuterated:
            hh_time = data["HH time"].dropna()
            hd_time = data["HD time"].dropna()
        else:
            hd_time = data["HD time"].dropna()
            dd_time = data["DD time"].dropna()

        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)

        if not self.perdeuterated:
            axs[0].hist(hh_time, bins=128, label="HH time", color="dimgray")
            axs[1].hist(hd_time, bins=128, label="HD time", color="darkgray")
        else:
            axs[0].hist(hd_time, bins=128, label="HD time", color="dimgray")
            axs[1].hist(dd_time, bins=128, label="DD time", color="darkgray")

        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

        fig.supxlabel("Time [s]")
        fig.supylabel("Occurrences")

        for ax in axs:
            ax.label_outer()

        axs[0].legend()
        axs[1].legend()
        plt.show()

    def histogram_dissociation_positions(self):
        mean_positions = self.dissociation_positions.mean()
        std_positions = self.dissociation_positions.std()

        index = np.arange(len(mean_positions.index))

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index,
            mean_positions.values,
            yerr=std_positions.values,
            capsize=5,
            color="dimgray",
        )
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence")
        plt.xticks(index, mean_positions.index)
        plt.show()

    def histogram_diss_analysis(self, mol: list = ["H", "D"]):
        """Compute data"""
        data = self.data.dropna(subset=["Diss pos"])
        data = data[data["Diss atom"].isin(mol)]

        counts = data.groupby("Diss pos").size().reset_index(name="Occurences")
        counts.set_index("Diss pos", inplace=True)

        for key, value in self.symmetries.items():
            if isinstance(value, list):
                counts.loc[key] = counts.loc[value].sum()
                counts.drop(value, inplace=True)
            else:
                counts.loc[key] = counts.loc[value]
                counts.drop(value, inplace=True)

        index = np.arange(len(counts.index))

        ratios = self.ratios.loc["Avg"] if len(mol) == 2 else self.ratios.loc[mol[0]]

        """ Plot histogram of Occurences / Dissociation position"""
        counts = counts.reindex(self.symmetries.keys())

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, counts["Occurences"].values, color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurences")
        plt.xticks(index, counts.index)
        plt.show()

        plt.clf()

        """ Plot Occurence divided by time spent on position """
        result = pd.DataFrame()
        for pos in counts.index:
            result.loc[pos, "Occurences"] = (
                counts.loc[pos, "Occurences"] / self.position_times.mean().loc[str(pos)]
            )

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, result["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by time spent")
        plt.xticks(index, result.index)
        plt.show()

        plt.clf()

        """ Plot Occurence divided by reaction rate """

        occurence_rate = counts.div(ratios.iloc[0], axis=0)

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, occurence_rate["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by reaction rate")
        plt.xticks(index, occurence_rate.index)
        plt.show()

        plt.clf()

        """ Plot Occurence divided by reaction rate and time spent on position """
        result = result.div(ratios.iloc[0], axis=0)

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, result["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by reaction rate and time spent")
        plt.xticks(index, result.index)
        plt.show()

        plt.clf()

    def mass_spectra(self, parent_ion):
        no_diss = (self.data["Diss atom"].isna()).sum() / len(self.data) * 100
        h_diss = (self.data["Diss atom"] == "H").sum() / len(self.data) * 100
        d_diss = (self.data["Diss atom"] == "D").sum() / len(self.data) * 100

        masses = [parent_ion, parent_ion - 1, parent_ion - 2]

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        if not self.perdeuterated:
            plt.stem(
                masses,
                [no_diss, h_diss, d_diss],
                basefmt=" ",
                markerfmt="",
                linefmt="black",
            )
        else:
            plt.stem(
                masses,
                [no_diss, h_diss, d_diss],
                basefmt=" ",
                markerfmt="",
                linefmt="black",
            )

        plt.stem(parent_ion + 0.1, 100, basefmt=" ", markerfmt="", linefmt="dimgray")
        plt.xlabel("m/z (amu/e)")
        plt.ylabel("Intensity [%]")
        plt.xlim(min(masses) - 2, max(masses) + 2)
        plt.ylim(0, 105)
        plt.show()

    def data_tables(self):
        print(f"Data analysis for {self.name}")
        data = self.dissociation_counts
        data["Dissociated"] = (data["H"] + data["D"]) / data.sum(axis=1)
        data["%H"] = data["H"] / (data["H"] + data["D"])
        data["%D"] = data["D"] / (data["H"] + data["D"])
        print(data)

        print(
            f"Percentage dissociated = {data["Dissociated"].mean():.1%} ± {data["Dissociated"].std():.1%}"
        )
        combined_std = (
            data["%H"].std() ** 2 + data["%D"].std() ** 2
        ) ** 0.5  # TODO: check this
        print(
            f"H/D loss ratio = {data["%H"].mean():.1%}/{data["%D"].mean():.1%} ± {combined_std:.1%}"
        )

        # print(f"H/D loss ratio = {H_dissociated:.1%}/{D_dissociated:.1%} ± {std_H_dissociated:.1%}/{std_D_dissociated:.1%}")

        data = self.data.dropna(subset=["Diss pos"])
        print(f"Median dissociation time = {Quantity(data['Diss time'].median(), 's')}")
        print(f"Median Scrambling hops = {Quantity(self.data['# hops'].median())}\n")

        h_median, d_median = (
            self.data["# H hops"].median(),
            self.data["# D hops"].median(),
        )
        percentile_diff = (h_median - d_median) / d_median
        print(f"Median H scrambling hops = {Quantity(h_median)}")
        print(f"Median D scrambling hops = {Quantity(d_median)}")
        print(f"Percentile difference = {percentile_diff:.2%}\n")

        HH_or_DD_time = (
            self.data["DD time"].median()
            if self.perdeuterated
            else self.data["HH time"].median()
        )
        HD_time = self.data["HD time"].median()
        percentile_diff = (HH_or_DD_time - HD_time) / HD_time
        print(
            f"Median {'DD' if self.perdeuterated else 'HH'} time = {Quantity(HH_or_DD_time, 's')}"
        )
        print(f"Median HD time: {Quantity(HD_time, 's')}")
        print(f"Percentile difference: {percentile_diff:.1%}\n")

        data = self.data
        data["%HH time"] = data["HH time"] * 9 / data["Diss time"]
        data["%DD time"] = data["DD time"] * 9 / data["Diss time"]
        data["%HD time"] = data["HD time"] / data["Diss time"]
        data["% Time"] = data["%HH time"] + data["%DD time"] + data["%HD time"]
        print(f"Debugging: {data["% Time"].mean():.2%}, {data["% Time"].std():.2%}")

        if not self.perdeuterated:
            print(
                f"HH mean relative time: {data["%HH time"].mean():.1%} ± {data["%HH time"].std():.1%}"
            )
        else:
            print(
                f"DD mean relative time: {data["%DD time"].mean():.1%} ± {data["%DD time"].std():.1%}"
            )

        print(
            f"HD mean relative time: {data["%HD time"].mean():.1%} ± {data["%HD time"].std():.1%}"
        )
        print("")


class DissociationAnalysis(Settings):
    def __init__(self, alpha, beta):
        super().__init__(figsize, dpi)
        self.alpha, self.beta = alpha, beta

    def histogram(self):

        alpha_positions = self.alpha.dissociation_positions
        beta_positions = self.beta.dissociation_positions

        bar_width = 0.4  # Width of the bars
        index = np.arange(len(alpha_positions.mean().index))  # The label locations

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index - bar_width / 2,
            alpha_positions.mean().values,
            bar_width,
            yerr=alpha_positions.std().values,
            label=self.alpha.name,
            capsize=5,
            color="dimgray",
        )
        plt.bar(
            index + bar_width / 2,
            beta_positions.mean().values,
            bar_width,
            yerr=beta_positions.std().values,
            label=self.beta.name,
            capsize=5,
            color="darkgray",
        )

        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrences")
        plt.xticks(index, alpha_positions.mean().index)
        plt.legend()
        plt.show()
