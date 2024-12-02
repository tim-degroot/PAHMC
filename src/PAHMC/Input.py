import os, re, sys
import numpy as np
import logging
import yaml


logger = logging.getLogger(__name__)


class Input_reader:
    """Takes the input file and initializes a parameter class for the simulation to be run"""

    def __init__(self, filename):
        try:
            file = open(filename, "r")
        except FileNotFoundError:
            logger.error("The provided input file '" + filename + "'  was not found.")
            sys.exit(2)

        with open(filename, "r") as file:
            data = yaml.load(file, Loader=yaml.CLoader)

        self.mol_name = data.get("Name")

        self.mol_edge = data.get("Initial edge").split()
        for i, edge in enumerate(self.mol_edge):
            edge = edge[1:-1].split(",")
            if len(edge) not in (1, 2, 3, 4):
                logger.error("Individual edges can only have up to 4 components. ")
                sys.exit(2)

            for j, atoms in enumerate(edge):

                if atoms not in ("H", "D", "HH", "HD", "DD"):
                    logger.error(
                        "Incorrect atoms defined, only no atom (0), hydrogen (H), deuterium (D), or aliphatic groups (HH, HD, DD) are currently supported. "
                    )
                    sys.exit(2)

                edge[j] = atoms
            self.mol_edge[i] = edge

            self.init_edge = [0] * len(self.mol_edge)

            for n in range(0, len(self.init_edge)):
                self.init_edge[n] = tuple(self.mol_edge[n])

            self.init_edge = tuple(self.init_edge)

        self.mol_edge_numbers = data.get("Initial edge numbers", None).split()
        for i, edge in enumerate(self.mol_edge_numbers):
            edge = edge[1:-1].split(",")
            if len(edge) is not len(self.mol_edge[i]):
                logger.error("Edges and edge numbering need to have the same shape. ")
                sys.exit(2)

            for j, number in enumerate(edge):
                edge[j] = number
            self.mol_edge_numbers[i] = edge

        if len(self.mol_edge_numbers) is not len(self.mol_edge):
            logger.error("Edges and edge numbering need to have the same size. ")
            sys.exit(2)

        self.init_edge_numbers = [0] * len(self.mol_edge_numbers)

        for n in range(0, len(self.init_edge_numbers)):
            self.init_edge_numbers[n] = tuple(self.mol_edge_numbers[n])

        self.init_edge_numbers = tuple(self.init_edge_numbers)

        self.mol_links = data.get("Initial links").split()
        for i, link in enumerate(self.mol_links):
            link = link[1:-1].split(",")

            if len(link) > 2:
                logger.error(
                    "Individual links can only have up to two components currently. "
                )
                sys.exit(2)
            for j, link_ in enumerate(link):
                if link_ != "0":
                    logger.error(
                        "Incorrect link specified, only empty tertiary carbons (0) are accepted. "
                    )
                    sys.exit(2)
                link[j] = link_
            self.mol_links[i] = link
        if len(self.mol_links) != len(self.mol_edge):
            logger.error("Edges and links need to have the same size. ")
            sys.exit(2)

        self.init_links = tuple(self.mol_links)

        self.mol_links_numbers = data.get("Initial link numbers").split()
        for i, link in enumerate(self.mol_links_numbers):
            link = link[1:-1].split(",")

            if len(link) is not len(self.mol_links[i]):
                logger.error("Links and link numbering need to have the same shape. ")
                sys.exit(2)

            for j, link_number in enumerate(link):
                link[j] = link_number
            self.mol_links_numbers[i] = link
        if len(self.mol_links_numbers) is not len(self.mol_links):
            logger.error("Links and link numbering need to have the same size. ")
            sys.exit(2)

        self.init_links_numbers = tuple(self.mol_links_numbers)

        self.cross_links = {}
        if data.get("Initial cross-links") is not None:
            cross_links = data.get("Initial cross-links").split()
            cross_links_tuples = [
                tuple(map(str, item.strip("()").split(","))) for item in cross_links
            ]

            for a, b in cross_links_tuples:
                self.cross_links[a] = b
                self.cross_links[b] = a

        try:
            self.energy = float(data.get("Energy"))
        except ValueError:
            logger.error(
                "Cannot understand the simulation energy given, it should be a single energy. "
            )

        try:
            self.iterations = int(data.get("Iterations"))
        except ValueError:
            logger.error("Invalid number of iterations given. ")
            sys.exit(2)

        try:
            self.t_max = float(data.get("Maximum time"))
        except ValueError:
            logger.error("Invalid maximum time given.")
            sys.exit(2)

        self.handling = data.get("Error handling").lstrip()
        if self.handling not in ("w", "q"):
            logger.error(
                "Invalid error handling mode given. Supported are: display warnings (w), display warnings and abort (q)."
            )
            sys.exit(2)

        self.rate_dir = data.get("Rates path")

        self.reactionrates = {}
        self.dE = {}
        self.reactionkeys = {}

        if not os.path.isdir(self.rate_dir):
            logger.error(
                f"The provided rates directory '{self.rate_dir}' was not found."
            )
            sys.exit(2)

        rates = data.get("Rates")

        for ratelist, ratefile in rates.items():
            filepath = os.path.join(self.rate_dir, ratefile)

            if not os.path.isfile(filepath):
                logger.error(f"The provided rate file '{filepath}' was not found.")
                sys.exit(2)

            rates = np.loadtxt(filepath, unpack=True, skiprows=2)

            with open(filepath) as f:
                firstline = next(f)
                delta = float(
                    re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*\b", firstline)[-1]
                )  #  (this regex finds all numbers in the first line, then saves only the last one)

            for reac in ratelist.split(","):
                self.reactionkeys[reac] = ratelist.split(",")[0]
                self.reactionrates[reac] = rates
                self.dE[reac] = delta

        N_rates = len(self.reactionrates.keys())
        logger.info(f"{N_rates} rate files successfully read.")
