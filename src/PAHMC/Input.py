# Import relevant libraries
import os, re, sys
import numpy as np
import logging
import yaml


logger = logging.getLogger(__name__)

# Import other relevant files
# (none currently)


# Create a class object to contain the parameters of the simulation
class Input_reader:
    """Takes the input file and initializes a parameter class for the simulation to be run"""

    def __init__(self, filename):
        # Try to open file, raise an error if file was not found
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

                # need to check for correct atom weights in edges
                # print(atoms)
                if atoms not in ("H", "D", "HH", "HD", "DD"):
                    logger.error(
                        "Incorrect atoms defined, only no atom (0), hydrogen (H), deuterium (D), or aliphatic groups (HH, HD, DD) are currently supported. "
                    )
                    sys.exit(2)

                # Assemble the seperate parts back together into a filled edge
                edge[j] = atoms
            # Edges to use for the simulation
            self.mol_edge[i] = edge

            # keep a record of the initial edge structure just in case
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

                # Assemble the seperate parts back together into integer filled edge
                edge[j] = number
            # Edges to use for the simulation
            self.mol_edge_numbers[i] = edge

        if len(self.mol_edge_numbers) is not len(self.mol_edge):
            logger.error("Edges and edge numbering need to have the same size. ")
            sys.exit(2)

        # keep a record of the initial edge structure just in case
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

        # Read in the range of the simulation, and convert to float/integer
        if isinstance(data.get("Energy range"), float):
            self.energy_range = data.get("Energy range")
        else:
            self.energy_range = data.get("Energy range").split(",")
            if len(self.energy_range) == 1:
                try:
                    self.energy_range = [float(self.energy_range[0])]
                except ValueError:
                    logger.error("The single energy should be a floating point value. ")
                    sys.exit(2)
            elif len(self.energy_range) == 3:
                try:
                    self.energy_range[0] = float(self.energy_range[0])
                    self.energy_range[1] = float(self.energy_range[1])
                    self.energy_range[2] = int(self.energy_range[2])
                except ValueError:
                    logger.error("Cannot understand the simulation range given. ")
                    sys.exit(2)
            else:
                logger.error(
                    "Cannot understand the simulation range given, it should be either a single energy or an energy range (min, max, nsteps). "
                )

        # Read in the number of simulations for each energy
        try:
            self.iterations = int(data.get("Iterations per energy"))
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

        if not os.path.isdir(self.rate_dir):
            logger.error(
                "The provided rates directory '" + self.rate_dir + "'  was not found."
            )
            sys.exit(2)

        rates = data.get("Rates")

        for ratelist, ratefile in rates.items():
            filepath = os.path.join(self.rate_dir, ratefile)

            if not os.path.isfile(filepath):
                logger.error(
                    "The provided rate file '" + filepath + "'  was not found."
                )
                sys.exit(2)

            rates = np.loadtxt(filepath, unpack=True, skiprows=2)

            with open(filepath) as f:
                firstline = next(f)  # Read in first line of the file
                # Then select the value of Delta out of the rate file by regular expression
                # (this regex finds all numbers in the first line, then saves only the last one)
                delta = float(
                    re.findall(r"[-+]?\d+[\.]?\d*[eE]?[-+]?\d*\b", firstline)[-1]
                )

            for reac in ratelist.split(","):
                self.reactionrates[reac] = rates
                self.dE[reac] = delta

            N_rates = len(self.reactionrates.keys())

            logger.info(str(N_rates) + " rate files sucessfully read.")
