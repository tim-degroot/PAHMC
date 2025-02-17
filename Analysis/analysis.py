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
from matplotlib.ticker import PercentFormatter

figsize: set = (5, 4)
dpi: int = 150


class Settings:
    def __init__(self, figsize: set = (5, 4), dpi: int = 150):
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
        steps: float = 0.4,
        sort_legend_by="colors",
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
        subset.sort()

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
            plt.xticks(np.arange(xlim[0], xlim[1] + steps, steps))

        plt.title(title)
        self._plot_plot(sort_legend_by=sort_legend_by)

    def quantify(self, param: str,filter_function=None, energy: float = 4.66):
        subset = self.rates.get(param)
        subset.sort()

        if filter_function is not None:
            subset = list(filter(filter_function, subset))

        result = pd.DataFrame(columns=[]).set_index(pd.Index([]))

        for i, filename in enumerate(subset):
            filepath = os.path.join(self.folder, filename)
            header = open(filepath, 'r').readline().strip()
            barrier = float(re.search(r"E0=([0-9.]+)", header).group(1))
            # print(barriers)
            data = np.loadtxt(filepath, skiprows=2)
            label_tuple = self._parse_label_tuple(filename)
            label = self._parse_label(filename)

            index = np.abs(data[:, 0] - energy).argmin()
            new_row = [label_tuple["From"], label_tuple["To"]]
            i = len(result)
            result.loc[i, ["From", "To"]] = new_row
            result.loc[i, f"{label_tuple["Atom"]} Rate"] = data[index, 1]
            result.loc[i, f"{label_tuple["Atom"]} Barrier"] = barrier
            # print(new_row, new_row.append(data[index, 1]) )

        return result


    def plot_pairs(
        self, subset: str, title: str = "", ylim: list = [], xlim: list = [0.0, 5.2], steps: float = 0.4,
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
                plt.xticks(np.arange(xlim[0], xlim[1] + steps, steps))

            plt.title(title)
            self._plot_plot()

    def _plot_plot(self, sort_legend_by="colors"):
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
        elif sort_legend_by == "colors":
            labels, handles = zip(
                *sorted(zip(labels, handles), key=lambda t: t[1].get_color())
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
    
    def _parse_label_tuple(self, label):
        base_name = os.path.splitext(label)[0]
        if "to" in base_name:
            parts = base_name.split("to")
            return {"From": parts[0][1:], "To": parts[1], "Atom": base_name[0]}
        elif "diss" in base_name:
            parts = base_name.split("diss")
            return {"From": parts[0][1:], "To": "diss", "Atom": base_name[0]}
            return f"{base_name[0]} {parts[0][1:]}"
        # f"{base_name[0]} {parts[0][1:]} to {parts[1]}"


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
        self.dissociation_counts, self.dissociation_raw, self.dissociation_positions = (
            self.process_output()
        )
        try:
            self.hops = self.process_hops()
        except FileNotFoundError:
            pass
        try:
            self.end_structures = self.process_endstruct()
        except FileNotFoundError:
            self.end_structures = None
        self.input = input
        self.input_data = Input.Input_reader(self.input)
        self.ratios = self.RRKM_rates()

    def position_times(self, mol: list = ["H", "D"], time=True):
        yaml = Input.Input_reader(self.input)
        full_list = list(set([x.split("to")[0] for x in self.hops.keys()]))
        from_list = [y for y in full_list if y[0] in mol]
        time_results = pd.DataFrame()
        occurrence_results = pd.DataFrame()

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
            summed_time_results = pd.Series(0, index=filtered_hops.index)
            summed_occurrence_results = pd.Series(0, index=filtered_hops.index)

            for r_key in filtered_hops.columns:
                column_data = filtered_hops[r_key]

                E = yaml.energy
                rate_idx = np.argmin((np.abs(yaml.reactionrates[r_key][0, :] - E)))
                r_rate = yaml.reactionrates[r_key][1, rate_idx]
                dt = 1 / np.sum(r_rate)

                summed_time_results += column_data * dt
                summed_occurrence_results += column_data

            time_results[origin] = summed_time_results
            occurrence_results[f"{origin}"] = summed_occurrence_results

            summed_time_results = pd.DataFrame(index=time_results.index)
            summed_occurrence_results = pd.DataFrame(index=occurrence_results.index)


        for num in set(
            col[1:] for col in time_results.columns if col.startswith(("H", "D"))
        ):
            summed_time_results[f"{num}"] = time_results.get(f"H{num}", 0) + time_results.get(
                f"D{num}", 0
            )
            summed_occurrence_results[f"{num}"] = occurrence_results.get(
                f"H{num}", 0
            ) + occurrence_results.get(f"D{num}", 0)

        # check if summes_results contains symmetries and sum them together
        for origin in summed_time_results.columns:
            for key, value in self.symmetries.items():
                if isinstance(value, int):
                    value = [value]
                values = [str(x) for x in value]
                if origin in values:
                    summed_time_results[str(key)] += summed_time_results[str(origin)]
                    summed_time_results.drop(
                        columns=[origin], inplace=True
                    )  # Drop the column

        for origin in summed_occurrence_results.columns:
            for key, value in self.symmetries.items():
                if isinstance(value, int):
                    value = [value]
                values = [str(x) for x in value]
                if origin in values:
                    summed_occurrence_results[str(key)] += summed_occurrence_results[str(origin)]
                    summed_occurrence_results.drop(
                        columns=[origin], inplace=True
                    )  # Drop the column

        factors = {}
        for key, value in self.symmetries.items():
            if isinstance(value, int):
                value = [value]
            factors[key] = 1+len(value)
        minimum = min(factors.values())
        factors = {str(k): int(v/minimum) for k, v in factors.items()}
        
        # print(summed_time_results.keys(), factors.keys())

        for key, factor in factors.items():
            summed_time_results[str(key)] /= factor
            summed_occurrence_results[str(key)] /= factor


        # summed_time_results = summed_time_results.div(summed_time_results.sum(axis=1), axis=0) # TODO: does this normalize the entire dataset or per row
        # summed_occurrence_results = summed_occurrence_results.div(summed_occurrence_results.sum(axis=1), axis=0)

        return summed_time_results if time else summed_occurrence_results

    def RRKM_rates(self):
        yaml = Input.Input_reader(self.input)

        rates = pd.DataFrame()

        for pos in self.symmetries.keys():
            for mol in ["H", "D"]:
                r_key = f"{mol}{pos}diss"
                E = yaml.energy
                rate_idx = np.argmin((np.abs(yaml.reactionrates[r_key][0, :] - E)))
                r_rate = yaml.reactionrates[r_key][1, rate_idx]
                rates.loc[mol, pos] = r_rate

        if self.perdeuterated:
            rates.loc["Avg"] = 10 / 11 * rates.loc["D"] + 1 / 11 * rates.loc["H"]
        else:
            rates.loc["Avg"] = 1 / 11 * rates.loc["D"] + 10 / 11 * rates.loc["H"]

        rates.loc["H/D ratio"] = rates.loc["H"] / rates.loc["D"]

        maximum = max(rates.loc["D"].max(), rates.loc["H"].max())
        rates.loc["D"] /= maximum
        rates.loc["H"] /= maximum
        rates.loc["Avg"] /= rates.loc["Avg"].max()

        return rates
    
    def formatted_ratios(self):
        yaml = Input.Input_reader(self.input)

        rates = pd.DataFrame()

        for pos in self.symmetries.keys():
            for mol in ["H", "D"]:
                r_key = f"{mol}{pos}diss"
                E = yaml.energy
                rate_idx = np.argmin((np.abs(yaml.reactionrates[r_key][0, :] - E)))
                r_rate = yaml.reactionrates[r_key][1, rate_idx]
                rates.loc[mol, pos] = r_rate
                print(r_rate)

        if self.perdeuterated:
            rates.loc["Avg"] = 10 / 11 * rates.loc["D"] + 1 / 11 * rates.loc["H"]
        else:
            rates.loc["Avg"] = 1 / 11 * rates.loc["D"] + 10 / 11 * rates.loc["H"]

        rates.loc["H/D ratio"] = rates.loc["H"] / rates.loc["D"]

        maximum = max(rates.loc["D"].max(), rates.loc["H"].max())
        rates.loc["D_normalized"] = rates.loc["D"] / maximum
        rates.loc["H_normalized"] = rates.loc["H"] / maximum
        rates.loc["Avg"] /= rates.loc["Avg"].max()

        return rates


    def process_data(self):
        data_frames = []

        for file in self.files["data"]:
            df = pd.read_csv(file, na_values=["None"])
            df.columns = df.columns.str.strip()
            df["MC#"] = df["MC#"].apply(lambda x: f"{x}-{file.split('/')[-1][0]}")
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

        dissociation_positions_merged = pd.DataFrame(0, index=df.index, columns=self.symmetries.keys())

        for key in self.symmetries.keys():
            dissociation_positions_merged[key] = df[key]

        for key, value in self.symmetries.items():
            if isinstance(value, list):
                for v in value:
                    if v in df.columns:
                        dissociation_positions_merged[key] += df[v]
                        # dissociation_positions_merged[key] = pd.concat(
                        #     [df[key], df[v]], ignore_index=True
                        # )
            else:
                if value in df.columns:
                    dissociation_positions_merged[key] += df[value]
                    # dissociation_positions_merged[key] = pd.concat(
                    #     [df[key], df[value]], ignore_index=True
                    # )

        return (
            pd.DataFrame(dissociation_counts),
            df,
            dissociation_positions_merged,
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
            df["ID"] = df["ID"].apply(lambda x: f"{x}-{file.split('/')[-1][0]}")
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
        # df.dropna(subset=["End position"], inplace=True)

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

    def plot_ratios(self, yscale="linear"):
        """Plot the dissociation rate relative to the highest rate"""
        energy = Input.Input_reader(self.input).energy
        index = np.arange(len(self.ratios.loc["H"].index))

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.scatter(
            index,
            self.ratios.loc["H"],
            label="H",
            color="dimgray",
        )
        # plt.scatter(
        #     index,
        #     self.ratios.loc["D"],
        #     label="D",
        #     color="darkgray",
        # )
        plt.xlabel("Dissociation position")
        plt.ylabel(f"Relative dissociation rate at {energy} eV")
        plt.xticks(index, self.ratios.loc["H"].keys())
        plt.legend()
        plt.ylim(0,1.1)
        plt.yscale(yscale)
        plt.show()

        plt.clf()

    def histogram_position_times(self, mol: list = ["H", "D"]):
        results = self.position_times(mol)

        mean_positions = results.sum()
        std_positions = np.sqrt(results.sum())

        if "8a" in mean_positions.index and "10a" in mean_positions.index:
            mean_positions["10a"] += mean_positions["8a"]
            std_positions["10a"] = np.sqrt(std_positions["10a"]**2 + std_positions["8a"]**2)
            mean_positions = mean_positions.drop("8a")
            std_positions = std_positions.drop("8a")


        index = np.arange(len(mean_positions.index))

        sorted_index = sorted(mean_positions.index)
        mean_positions = mean_positions[sorted_index]
        std_positions = std_positions[sorted_index]
        index = np.arange(len(sorted_index))


        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index,
            mean_positions.values,
            # yerr=std_positions.values,
            capsize=5,
            color="dimgray",
        )
        plt.xlabel("Position")
        plt.ylabel("Total time spent [s]")
        plt.xticks(index, mean_positions.index)
        # ax = plt.gca()
        # ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        # bar_width = 0.4

        # position_times_H = self.position_times("H")
        # position_times_D = self.position_times("D")

        # aligned_H, aligned_D = position_times_H.align(position_times_D, fill_value=0)

        # index = np.arange(len(aligned_H.sum().index))

        # plt.figure(figsize=self.figsize, dpi=self.dpi)
        # plt.bar(
        #     index - bar_width / 2,
        #     aligned_H.sum(),
        #     bar_width,
        #     # yerr=aligned_H.std(),
        #     label="H",
        #     color="dimgray",
        # )
        # plt.bar(
        #     index + bar_width / 2,
        #     aligned_D.sum(),
        #     bar_width,
        #     # yerr=aligned_D.std(),
        #     label="D",
        #     color="darkgray",
        # )
        # plt.xlabel("Position")
        # plt.ylabel("Time spent [s]")
        # plt.xticks(index, aligned_H.keys())
        # plt.legend()
        # # ax = plt.gca()
        # # ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        # plt.show()

        # plt.clf()

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
        data = data["Diss time"]  * 1000
        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.hist(data, bins=128, color="dimgray")
        plt.xlabel("Dissociation time [ms]")
        plt.ylabel("Frequency")
        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        # plt.gca().ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
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
        plt.ylabel("Frequency")
        plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        plt.show()

    def histogram_hops_comparison(self, ylim=None):
        data = self.data.dropna(subset=["Diss pos"])

        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        gs = fig.add_gridspec(2, hspace=0)
        axs = gs.subplots(sharex=True, sharey=True)

        axs[0].hist(data["# D hops"], bins=128, label="D hops", color="dimgray")
        axs[1].hist(data["# H hops"], bins=128, label="H hops", color="darkgray")

        fig.supxlabel("Number of scrambling hops")
        fig.supylabel("Frequency")

        for ax in axs:
            ax.label_outer()
            if ylim is not None:
                ax.set_ylim(ylim)

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
        fig.supylabel("Frequency")
        fig.subplots_adjust(left=0.15)

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
        result = self.dissociation_positions.sum()
        # std_positions = self.dissociation_positions.std()
        error = np.sqrt(result)

        index = np.arange(len(result.index))

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index,
            result.values,
            yerr=error.values,
            capsize=5,
            color="dimgray",
        )
        plt.xlabel("Dissociation position")
        plt.ylabel("Dissociation events")
        plt.xticks(index, result.index)
        plt.show()

    def histogram_diss_analysis(self, mol: list = ["H", "D"]):
        """Compute data"""
        data = self.data.copy()
        data = data[data["Diss atom"].isin(mol)]

        counts = data.groupby("Diss pos").size().reset_index(name="Occurences")
        counts.set_index("Diss pos", inplace=True)

        for key, value in self.symmetries.items():
            if isinstance(value, list):
                counts.loc[key] += counts.loc[value].sum()
                counts.drop(value, inplace=True)
            else:
                counts.loc[key] += counts.loc[value]
                counts.drop(value, inplace=True)

        index = np.arange(len(counts.index))

        ratios = self.ratios.loc["Avg"] if len(mol) == 2 else self.ratios.loc[mol]

        output = pd.DataFrame()

        """ Plot histogram of Occurences / Dissociation position"""
        counts = counts.reindex(self.symmetries.keys())
        counts["Occurences"] /= counts["Occurences"].sum()
        output["Occurences"] = counts

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, counts["Occurences"].values, color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurences")
        plt.xticks(index, counts.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        plt.clf()

        """ Plot Occurence divided by time spent on position """
        result = pd.DataFrame()
        for pos in counts.index:
            result.loc[pos, "Occurences"] = (
                counts.loc[pos, "Occurences"]
                / self.position_times(mol).mean().loc[str(pos)]
            )

        result["Occurences"] = result["Occurences"] / result["Occurences"].sum()
        output["Occurences/time"] = result

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, result["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by time spent (s$^{-1}$)")
        plt.xticks(index, result.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        plt.clf()

        """ Plot Occurence divided by reaction rate """

        occurence_rate = counts.div(ratios, axis=0)
        occurence_rate["Occurences"] = (
            occurence_rate["Occurences"] / occurence_rate["Occurences"].sum()
        )
        output["Occurence/rate"] = occurence_rate

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, occurence_rate["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by reaction rate")
        plt.xticks(index, occurence_rate.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        plt.clf()

        """ Plot Occurence divided by reaction rate and time spent on position """
        result = result.div(ratios, axis=0)

        result["Occurences"] = result["Occurences"] / result["Occurences"].sum()
        output["Occurence/rate/time"] = result

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(index, result["Occurences"], color="dimgray")
        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrence divided by reaction rate and time spent")
        plt.xticks(index, result.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        plt.clf()

        return output

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

        plt.stem(parent_ion + 0.1, 100, basefmt=" ", markerfmt="", linefmt="darkgray")
        plt.xlabel("m/z (amu/e)")
        plt.ylabel("Intensity [%]")
        plt.xlim(min(masses) - 2, max(masses) + 2)
        plt.ylim(0, 105)
        plt.show()

    def histogram_diss_pair_analysis(self):
        """Compute data"""
        data = self.data.copy()
        data = data[data["Diss atom"].isin(["H", "D"])]

        counts = (
            data.groupby(["Diss pos", "Diss atom"])
            .size()
            .reset_index(name="Occurences")
        )
        counts["Occurences"] = counts["Occurences"].astype(float)
        counts = counts.pivot(
            index="Diss pos", columns="Diss atom", values="Occurences"
        ).fillna(0)

        for column in counts.columns:
            for key, value in self.symmetries.items():
                if isinstance(value, list):
                    counts.loc[key, column] += counts.loc[value, column].sum()
                    counts.loc[key, column] /= 1 + len(value)
                else:
                    counts.loc[key, column] += counts.loc[value, column]
                    counts.loc[key, column] /= 2

        symmetry_values = []
        for value in self.symmetries.values():
            if isinstance(value, list):
                symmetry_values.extend(value)
            else:
                symmetry_values.append(value)
        counts.drop(symmetry_values, inplace=True)

        index = np.arange(len(counts.index))
        bar_width = 0.4

        """ Plot histogram of Occurences / Dissociation position"""  # TODO: Add error bars
        counts = counts.reindex(self.symmetries.keys())
        counts["Error_H"] = np.sqrt(counts["H"])
        counts["Error_D"] = np.sqrt(counts["D"])

        sum_H = counts["H"].sum()
        sum_D = counts["D"].sum()

        counts["H"] /= sum_H
        counts["D"] /= sum_D

        counts["Error_H"] = np.sqrt(
            (counts["Error_H"] / sum_H) ** 2
            + (counts["H"] * np.sqrt(sum_H) / sum_H**2) ** 2
        )
        counts["Error_D"] = np.sqrt(
            (counts["Error_D"] / sum_D) ** 2
            + (counts["D"] * np.sqrt(sum_D) / sum_D**2) ** 2
        )

        plt.figure(figsize=self.figsize, dpi=self.dpi)

        plt.bar(
            index - bar_width / 2,
            counts["H"].values,
            bar_width,
            yerr=counts["Error_H"].values,
            label="H",
            color="dimgray",
        )
        plt.bar(
            index + bar_width / 2,
            counts["D"].values,
            bar_width,
            yerr=counts["Error_D"].values,
            label="D",
            color="darkgray",
        )

        plt.xlabel("Dissociation position")
        plt.ylabel("Dissociation events")
        plt.xticks(index, counts.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.legend()
        plt.show()

        plt.clf()

        """ Plot Occurence divided by time spent on position """
        # result = pd.DataFrame()
        # for mol in ["H", "D"]:
        #     result[f"Error_{mol}"] = counts[f"Error_{mol}"]
        #     for pos in counts.index:
        #         result.loc[pos, mol] = (
        #             counts.loc[pos, mol] / self.position_times(mol).mean().loc[str(pos)]
        #         )
        #     sum = result[mol].sum()
        #     result[mol] /= sum
        #     # print(counts[f"Error_{mol}"])
        #     counts[f"Error_{mol}"] = np.sqrt(
        #         (counts[f"Error_{mol}"] / sum) ** 2
        #         + (result[mol] * np.sqrt(sum) / sum**2) ** 2
        #     )
            # print(counts[f"Error_{mol}"])

        # plt.figure(figsize=self.figsize, dpi=self.dpi)
        # plt.bar(
        #     index - bar_width / 2,
        #     result["H"].values,
        #     bar_width,
        #     yerr=result["Error_H"].values,
        #     label="H",
        #     color="dimgray",
        # )
        # plt.bar(
        #     index + bar_width / 2,
        #     result["D"].values,
        #     bar_width,
        #     yerr=result["Error_D"].values,
        #     label="D",
        #     color="darkgray",
        # )
        # plt.xlabel("Dissociation position")
        # plt.ylabel("Occurrence divided by time spent (s$^{-1}$)")
        # plt.xticks(index, result.index)
        # ax = plt.gca()
        # ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        # plt.legend()
        # plt.show()

        # plt.clf()

        """ Plot Occurence divided by reaction rate """
        # occurence_rate = pd.DataFrame()

        # for mol in ["H", "D"]:
        #     for pos in counts.index:
        #         occurence_rate.loc[pos, mol] = (
        #             counts.loc[pos, mol] / self.ratios.loc[mol, pos]
        #         )
        #     occurence_rate[mol] /= occurence_rate[mol].sum()

        # plt.figure(figsize=self.figsize, dpi=self.dpi)
        # plt.bar(
        #     index - bar_width / 2,
        #     occurence_rate["H"].values,
        #     bar_width,
        #     label="H",
        #     color="dimgray",
        # )
        # plt.bar(
        #     index + bar_width / 2,
        #     occurence_rate["D"].values,
        #     bar_width,
        #     label="D",
        #     color="darkgray",
        # )
        # plt.xlabel("Dissociation position")
        # plt.ylabel("Occurrence divided by reaction rate")
        # plt.xticks(index, occurence_rate.index)
        # ax = plt.gca()
        # ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        # plt.legend()
        # plt.show()

        # plt.clf()

        """ Plot Occurence divided by reaction rate and time spent on position """
        # Divide by time spent
        # print(counts)
        for mol in ["H", "D"]:
            max_position_time = self.position_times(mol).sum().max()
            for pos in counts.index:
                value = self.position_times(mol).sum().loc[str(pos)]
                # counts.loc[pos, mol] /= (value / max_position_time)
                value_error = np.sqrt(value)
                max_position_time_error = np.sqrt(max_position_time)
                error = (value / max_position_time) * np.sqrt((value_error / value) ** 2 + (max_position_time_error / max_position_time) ** 2)
                counts.loc[pos, mol] /= (value / max_position_time)
                # counts.loc[pos, f"Error_{mol}"] = counts.loc[pos, f"Error_{mol}"] * error / (value / max_position_time)
                # counts.loc[pos, f"Error_{mol}"] = counts.loc[pos, mol] * np.sqrt((error/))
                counts.loc[pos, mol] /= self.ratios.loc[mol, pos]
            counts[mol] /= counts[mol].sum()
        # print(counts)


        # def baseline(x):
        #     return x - 1/len(occurence_rate.index)

        # result = result.applymap(baseline)

        # print(result)

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.axhline(y=1/len(counts["D"].values), color='black', linestyle='--', linewidth=1)
        plt.bar(
            index - bar_width / 2,
            counts["H"].values,
            bar_width,
            yerr=counts["Error_H"].values,
            label="H",
            color="dimgray",
        )
        plt.bar(
            index + bar_width / 2,
            counts["D"].values,
            bar_width,
            yerr=counts["Error_D"].values,
            label="D",
            color="darkgray",
        )
        plt.xlabel("Dissociation position")
        plt.ylabel("Dissociation events divided by rate and time")
        plt.xticks(index, counts.index)
        ax = plt.gca()
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.legend()
        plt.show()

        plt.clf()

    def plot_dissociation_details(self, ax_1=('linear', (0,0.5)), ax_2=('linear', (0,1.1))):
        """Compute data"""
        data = self.data.copy()
        data = data[data["Diss atom"].isin(["H", "D"])]

        counts = (
            data.groupby(["Diss pos", "Diss atom"])
            .size()
            .reset_index(name="Occurences")
        )
        counts["Occurences"] = counts["Occurences"].astype(float)
        counts = counts.pivot(
            index="Diss pos", columns="Diss atom", values="Occurences"
        ).fillna(0)

        for column in counts.columns:
            for key, value in self.symmetries.items():
                if isinstance(value, list):
                    counts.loc[key, column] += counts.loc[value, column].sum()
                    counts.loc[key, column] /= 1 + len(value)
                else:
                    counts.loc[key, column] += counts.loc[value, column]
                    counts.loc[key, column] /= 2

        symmetry_values = []
        for value in self.symmetries.values():
            if isinstance(value, list):
                symmetry_values.extend(value)
            else:
                symmetry_values.append(value)
        counts.drop(symmetry_values, inplace=True)

        index = np.arange(len(counts.index))
        bar_width = 0.4

        """ Plot histogram of Occurences / Dissociation position"""  # TODO: Add error bars
        counts = counts.reindex(self.symmetries.keys())
        counts["Error_H"] = np.sqrt(counts["H"])
        counts["Error_D"] = np.sqrt(counts["D"])

        sum_H = counts["H"].sum()
        sum_D = counts["D"].sum()

        counts["H"] /= sum_H
        counts["D"] /= sum_D

        counts["Error_H"] = np.sqrt(
            (counts["Error_H"] / sum_H) ** 2
            + (counts["H"] * np.sqrt(sum_H) / sum_H**2) ** 2
        )
        counts["Error_D"] = np.sqrt(
            (counts["Error_D"] / sum_D) ** 2
            + (counts["D"] * np.sqrt(sum_D) / sum_D**2) ** 2
        )
        
        ratios = {
            "H": self.ratios.loc["H"],# / self.ratios.loc["H"].sum(),
            "D": self.ratios.loc["D"],# / self.ratios.loc["D"].sum(),
        }

        position_times = {}
        maximum = max(max(self.position_times("H", True).sum()), max(self.position_times("D", True).sum()))
        # print(maximum)
        for mol in ["H", "D"]:
            position_times[mol] = []
            mol_data = self.position_times(mol).sum()
            for pos in self.symmetries.keys():
                position_times[mol].append(
                    mol_data.loc[str(pos)] / maximum # / maximum
                )

        position_occurences = {}
        maximum = max(max(self.position_times("H", False).mean()), max(self.position_times("D", False).mean()))
        for mol in ["H", "D"]:
            position_occurences[mol] = []
            mol_data = self.position_times(mol, False).mean()
            for pos in self.symmetries.keys():
                position_occurences[mol].append(
                    mol_data.loc[str(pos)] / maximum
                )


        fig, ax = plt.subplots(
            figsize=(self.figsize[0], self.figsize[1]),
            dpi=self.dpi,
            layout="constrained",
        )
        # ax = fig.add_subplot(111)

        ax.bar(
            index - bar_width / 2,
            counts["H"].values,
            bar_width,
            yerr=counts["Error_H"].values,
            label="Hydrogen dissocication",
            color="dimgray",
        )
        ax.bar(
            index + bar_width / 2,
            counts["D"].values,
            bar_width,
            yerr=counts["Error_D"].values,
            label="Deuterium dissociation",
            color="darkgray",
        )

        ax2 = ax.twinx()
        ax2.errorbar(
            index - bar_width / 2,
            ratios["H"],
            label="Relative RRKM dissociation rate",
            fmt="v",
            color="#690000",
        )
        ax2.errorbar(
            index + bar_width / 2,
            ratios["D"],
            # label="Deuterium relative RRKM rates",
            fmt="v",
            color="#a90000",
        )

        ax2.errorbar(
            index - bar_width / 2,
            position_times["H"],
            label="Relative time spent at position",
            fmt="s",
            color="#006900",
        )
        ax2.errorbar(
            index + bar_width / 2,
            position_times["D"],
            # label="Deuterium time spent at positions",
            fmt="s",
            color="#00a900",
        )

        # ax2.errorbar(
        #     index - bar_width / 2,
        #     position_occurences["H"],
        #     label="Hops away from position",
        #     fmt="o",
        #     color="#000069",
        # )
        # ax2.errorbar(
        #     index + bar_width / 2,
        #     position_occurences["D"],
        #     # label="Deuterium time spent at positions",
        #     fmt="o",
        #     color="#0000a9",
        # )

        ax.set_xlabel("Dissociation position")
        ax.set_ylabel("Dissociation events")
        ax2.set_ylabel("Relative time spent / Relative RRKM rate")
        plt.xticks(index, counts.index)
        fig.legend(loc="outside lower center", ncol=2)
        ax.set_yscale(ax_1[0])
        ax2.set_yscale(ax_2[0])
        ax.set_ylim(ax_1[1])
        ax2.set_ylim(ax_2[1])
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        # ax2.yaxis.set_major_formatter(PercentFormatter(xmax=1))
        plt.show()

        plt.clf()

    def data_tables(self):
        print(f"Data analysis for {self.name}")
        data = self.dissociation_counts
        data["Dissociated"] = (data["H"] + data["D"]) / (
            data["H"] + data["D"] + data["None"]
        )
        data["%H"] = data["H"] / (data["H"] + data["D"])
        data["%D"] = data["D"] / (data["H"] + data["D"])

        print(
            f"Percentage dissociated = {data["Dissociated"].mean():.1%} ± {data["Dissociated"].std():.1%}"
        )
        combined_std = (data["%H"].std() ** 2 + data["%D"].std() ** 2) ** 0.5
        print(
            f"H/D loss ratio = {data["%H"].mean():.1%}/{data["%D"].mean():.1%} ± {combined_std:.1%}"
        )
        ratio = "5%/95%" if self.perdeuterated else "95%/5%"
        print(f"H/D loss ratio (assuming full scrambling and 50%/50% rates) = {ratio}")
        if self.perdeuterated:
            H_alpihatic_chance = 0.1
            RRKM_ratio = self.ratios.loc["H/D ratio"].mean()
            H_diss_chance = RRKM_ratio / (RRKM_ratio + 1)
            H_loss_chance = H_alpihatic_chance * H_diss_chance
            print(
                f"H/D loss ratio (assuming full scrambling and RRKM ratios) = {H_loss_chance:.2%}/{1-H_loss_chance:.2%}\n"
            )
        else:
            D_alpihatic_chance = 0.1
            RRKM_ratio = self.ratios.loc["H/D ratio"].mean()
            D_diss_chance = 1 - (RRKM_ratio / (RRKM_ratio + 1))
            D_loss_chance = D_alpihatic_chance * D_diss_chance
            print(
                f"H/D loss ratio (assuming full scrambling and RRKM ratios) = {1-D_loss_chance:.2%}/{D_loss_chance:.2%}\n"
            )

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
            print(f"{(data["HH time"] * 9).mean()/data["Diss time"].mean()}")
            print(
                f"HH mean relative time: {data["%HH time"].mean():.1%} ± {data["%HH time"].std():.1%} / Median: {data["%HH time"].median():.1%}"
            )
        else:
            print(
                f"DD mean relative time: {data["%DD time"].mean():.1%} ± {data["%DD time"].std():.1%} / Median: {data["%DD time"].median():.1%}"
            )
        print(f"{(data["HD time"]).mean()/data["Diss time"].mean()}")
        print(
            f"HD mean relative time: {data["%HD time"].mean():.1%} ± {data["%HD time"].std():.1%} / Median: {data["%HD time"].median():.1%}"
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

    def histogram(self):

        alpha_positions = self.alpha.dissociation_raw.sum()
        beta_positions = self.beta.dissociation_raw.sum()

        symmetries = {}
        for key, value in self.alpha.symmetries.items():
            if isinstance(value, list):
                for v in value:
                    symmetries[v] = key
            else:
                symmetries[value] = key

        df = pd.DataFrame(list(alpha_positions.items()), columns=["Position", "Value"])
        df["Symmetrical Position"] = (
            df["Position"].map(symmetries).fillna(df["Position"]).astype(int)
        )
        alpha = df.groupby("Symmetrical Position")["Value"].agg(["mean", "std"])

        print(df)

        df = pd.DataFrame(list(beta_positions.items()), columns=["Position", "Value"])
        df["Symmetrical Position"] = (
            df["Position"].map(symmetries).fillna(df["Position"]).astype(int)
        )
        beta = df.groupby("Symmetrical Position")["Value"].agg(["mean", "std"])

        bar_width = 0.4  # Width of the bars
        index = np.arange(len(alpha["mean"].index))  # The label locations

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.bar(
            index - bar_width / 2,
            alpha["mean"],
            bar_width,
            yerr=alpha["std"],
            label=self.alpha.name,
            capsize=5,
            color="dimgray",
        )
        plt.bar(
            index + bar_width / 2,
            beta["mean"],
            bar_width,
            yerr=beta["std"],
            label=self.beta.name,
            capsize=5,
            color="darkgray",
        )

        plt.xlabel("Dissociation position")
        plt.ylabel("Occurrences")
        plt.xticks(index, alpha["mean"].index)
        plt.legend()
        plt.show()

    def histogram_time(self):
        alpha_data = self.alpha.data.dropna(subset=["Diss pos"])["Diss time"]
        beta_data = self.beta.data.dropna(subset=["Diss pos"])["Diss time"]

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.hist(alpha_data, bins=128, color="dimgray", alpha=0.5, label=self.alpha.name)
        plt.hist(beta_data, bins=128, color="darkgray", alpha=0.5, label=self.beta.name)
        plt.xlabel("Dissociation time [s]")
        plt.ylabel("Frequency")
        plt.legend()
        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        plt.show()

    def histogram_hops(self):
        alpha_data = self.alpha.data.dropna(subset=["Diss pos"])["# hops"]
        beta_data = self.beta.data.dropna(subset=["Diss pos"])["# hops"]

        plt.figure(figsize=self.figsize, dpi=self.dpi)
        plt.hist(alpha_data, bins=128, color="dimgray", alpha=0.5, label=self.alpha.name)
        plt.hist(beta_data, bins=128, color="darkgray", alpha=0.5, label=self.beta.name)
        plt.xlabel("Total number of scrambling hops")
        plt.ylabel("Occurrences")
        plt.legend()
        plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
        plt.show()