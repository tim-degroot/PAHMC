import numpy as np
from matplotlib import pyplot as plt
import logging

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

    fname = outfilename.split(".")[0]

    data_file = open(f"{fname}_{E}_data.log", "a")

    data_file.write(
        "MC{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            mc,
            dissociation_atom,
            dissociation_position,
            dissociation_time,
            N_scramble_hops,
            N_D_hops,
            cross_hops,
            HH_time,
            HD_time,
            DD_time,
        )
    )

    data_file.close()


def Structure_Output(outfilename, E, iteration, molecule):
    fname = outfilename.split(".")[0]

    struct_file = open(f"{fname}_{E}_iteration_{iteration}_mol_structures.log", "a")

    struct_file.write(str(molecule.edges))
    struct_file.write("\n")
    struct_file.write(str(molecule.links))
    struct_file.write("\n")
    struct_file.write("\n")

    struct_file.close()


def End_Structures_Output(outfilename, E, edge, mc):
    fname = outfilename.split(".")[0]

    endstruct_file = open(f"{fname}_{E}_end_structures.out", "a")

    endstruct_file.write(f"MC{mc}\t{edge}")
    endstruct_file.write("\n")

    endstruct_file.close()


def STD_Output(
    outfilename,
    dissociation_atoms,
    dissociation_positions,
    dissociation_times,
    N_scramble_hops,
    N_D_hops,
):
    """Creates an output file with a summary of the results"""

    # Currently hardcoded, needs to be incorporated better in future
    time_bin_size = 1e-06
    hops_bin_size = 50000
    max_time = 0.001

    logger.info("Writing output file...")
    outfile = open(outfilename, "w")

    fname = outfilename.split(".")[0]

    energies = list(dissociation_times.keys())

    for E in energies:

        outfile.write(f"Dissociations for Energy {E}:\n")
        outfile.write(
            f"H: {dissociation_atoms[E]['H']}, D: {dissociation_atoms[E]['D']}, None: {dissociation_atoms[E]['None']}\n"
        )

        outfile.write(f"Dissociation positions for Energy {E}:\n")
        for key in list(dissociation_positions[E].keys()):
            outfile.write(f"{key}: {dissociation_positions[E][key]}\n")

        # make bins for histogram
        time_bins = np.arange(0, max_time, time_bin_size)
        hops_bins = np.arange(0, np.amax(N_scramble_hops[E]), hops_bin_size)

        # Save to file

        # time_file = open(fname+'_'+str(E)+'_diss_time.out', 'w')
        # hops_file = open(fname+'_'+str(E)+'_hops.out', 'w')
        # D_hops_file = open(fname+'_'+str(E)+'_D_hops.out', 'w')

        # np.savetxt(time_file, dissociation_times[E])
        # np.savetxt(hops_file, N_scramble_hops[E])
        # np.savetxt(D_hops_file, N_D_hops[E])

        # time_file.close()
        # hops_file.close()

        logger.info("Plotting results...")
        plt.hist(dissociation_times[E], bins=time_bins)
        plt.xlabel("Dissociation Time [s]")
        plt.ylabel("Frequency")
        plt.title("Dissociation times for energy " + str(E) + " eV")
        plt.savefig(fname + "_diss_time.png")
        plt.clf()

        plt.hist(N_scramble_hops[E], bins=hops_bins)
        plt.xlabel("Number of hops")
        plt.ylabel("Frequency")
        plt.title("Number of hops for energy " + str(E) + " eV")
        plt.savefig(fname + "_hops.png")

    outfile.close()
