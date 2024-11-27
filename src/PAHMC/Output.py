import numpy as np
from matplotlib import pyplot as plt
import logging
import os

logger = logging.getLogger(__name__)


def Data_Output(
    outfilename,
    E,
    dissociation_atom,
    dissociation_position,
    dissociation_time,
    N_scramble_hops,
    N_D_hops,
    HH_time,
    HD_time,
    DD_time,
    mc,
    cross_hops,
):

    file_name = outfilename.split(".")[0]
    file_path = f"{file_name}_data.csv"
    file_exists = os.path.isfile(file_path)

    with open(file_path, "a") as data_file:
        headers = [
            "MC#",
            "Diss atom",
            "Diss pos",
            "Diss time",
            "# hops",
            "# D hops",
            "# cross hops",
            "HH time",
            "HD time",
            "DD time",
        ]

        if not file_exists:
            data_file.write(",".join(headers) + "\n")

        values = [
            f"MC{mc}",
            dissociation_atom,
            dissociation_position,
            round(dissociation_time, 10),
            N_scramble_hops,
            N_D_hops,
            cross_hops,
            round(HH_time, 10),
            round(HD_time, 10),
            round(DD_time, 10),
        ]

        data_file.write(",".join([str(value) for value in values]) + "\n")


def Structure_Output(outfilename, E, iteration, molecule):
    file_name = outfilename.split(".")[0]
    file_path = f"{file_name}_iteration_{iteration}_mol_structures.log"

    with open(file_path, "a") as struct_file:
        struct_file.write(f"{str(molecule.edges)}\n{str(molecule.links)}\n\n")


def End_Structures_Output(outfilename, E, edge, mc):
    file_name = outfilename.split(".")[0]
    file_path = f"{file_name}_end_structures.out"

    with open(file_path, "a") as endstruct_file:
        endstruct_file.write(f"MC{mc}\t{edge}\n")


def STD_Output(
    outfilename,
    dissociation_atoms,
    dissociation_positions,
    dissociation_times,
    N_scramble_hops,
    N_D_hops,
):
    """Creates an output file with a summary of the results"""

    # TODO: Currently hardcoded, needs to be incorporated better in future
    time_bin_size = 1e-06
    hops_bin_size = 50000
    max_time = 0.02

    logger.info("Writing output file...")
    outfile = open(outfilename, "w")

    with open(outfilename, "w") as outfile:
        energies = list(dissociation_times.keys())

        for E in energies:
            outfile.write(f"Dissociations for Energy {E}:\n")
            outfile.write(
                f"H: {dissociation_atoms[E]['H']}, D: {dissociation_atoms[E]['D']}, None: {dissociation_atoms[E]['None']}\n"
            )

            outfile.write(f"Dissociation positions for Energy {E}:\n")
            for key in list(dissociation_positions[E].keys()):
                outfile.write(f"{key}: {dissociation_positions[E][key]}\n")

    file_name = outfilename.split(".")[0]

    time_bins = np.arange(0, max_time, time_bin_size)
    hops_bins = np.arange(0, np.amax(N_scramble_hops[E]), hops_bin_size)

    logger.info("Plotting results...")
    plt.hist(dissociation_times[E], bins=time_bins)
    plt.xlabel("Dissociation Time [s]")
    plt.ylabel("Frequency")
    plt.title("Dissociation times for energy " + str(E) + " eV")
    plt.savefig(file_name + "_diss_time.png")
    plt.clf()

    plt.hist(N_scramble_hops[E], bins=hops_bins)
    plt.xlabel("Number of hops")
    plt.ylabel("Frequency")
    plt.title("Number of hops for energy " + str(E) + " eV")
    plt.savefig(file_name + "_hops.png")


def Hops_Output(outfilename, mc, key_hops, E):
    file_name = outfilename.split(".")[0]
    file_path = f"{file_name}_key_hops.csv"
    file_exists = os.path.isfile(file_path)

    with open(file_path, "a") as hops_file:
        sorted_keys = sorted(key_hops.keys())
        headers = ["MC#"] + sorted_keys

        if not file_exists:
            hops_file.write(",".join(headers) + "\n")

        sorted_keys = sorted(key_hops.keys())
        values = [f"MC{str(mc)}"] + [str(key_hops[key]) for key in sorted_keys]
        hops_file.write(",".join(values) + "\n")
