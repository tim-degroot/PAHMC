#! /usr/bin/env python3

import sys, getopt
import logging
import argparse

import numpy as np
import multiprocessing as mp
import copy as cp
import random

from PAHMC.Functions import Possible_reactions
import PAHMC.Reactions as React
from PAHMC.Functions import Choose_reaction
from PAHMC.Input import Input_reader
from PAHMC.Molecule import Molecule
from PAHMC.Ratecheck import Check_available_rates
from PAHMC.Output import STD_Output
from PAHMC.Output import Data_Output
from PAHMC.Output import Structure_Output
from PAHMC.Output import End_Structures_Output


def run_iterations(start, end, input, value, molecule, queue, outputfile):
    logger.info(f"Running iterations {start + 1}-{end} out of {input.iterations}...")

    processes = [
        mp.Process(
            target=Parallel_Single_MC,
            args=(
                value,
                input.t_max,
                cp.deepcopy(molecule) if end == input.iterations else molecule,
                input,
                queue,
                j_iter,
                outputfile,
            ),
        )
        for j_iter in range(start, end)
    ]

    for process in processes:
        process.start()

    for process in processes:
        process.join()


def Parallel_Single_MC(E, max_time, molecule, rates, queue, j_iter, outfilename):
    """Run a single MC"""

    print(f"MC {j_iter} starting", flush=True)

    # Make a list of all reaction keys that have rates specified
    specified_rates = list(rates.reactionrates.keys())

    energy = cp.copy(E)
    time = 0
    total_hops = 0
    D_hops = 0
    diss_atom = None
    diss_position = None

    while time < max_time:

        if j_iter == 1:
            Structure_Output(outfilename, E, j_iter, molecule)
        # Determine the possible reactions (and aliphatic site surroundings)
        molecule.possible_reactions = Possible_reactions(molecule, specified_rates)

        # Choose reaction from the possibilities
        reactionkey, dt = Choose_reaction(energy, molecule.possible_reactions, rates)

        # Error if there's no dissociation before the rates go to 0
        if reactionkey == None:
            logger.error("Ran out of reaction rates before dissociation.")
            break

        # Update time
        time += dt

        # Update energy
        energy -= rates.dE[reactionkey]

        if molecule.al_place == "e":
            if len(molecule.index[0]) == 1:
                m = molecule.index[0][0]
                n = 0
            elif len(molecule.index[0]) == 2:
                m, n = molecule.index[0]

            if molecule.edges[m][n] == "HH":
                molecule.HH_time += dt
            elif molecule.edges[m][n] in ("HD", "DH"):
                molecule.HD_time += dt
            elif molecule.edges[m][n] == "DD":
                molecule.DD_time += dt

        # Carry out the reaction
        if "to" in reactionkey:
            # First some bookkeeping (saving the position of the aliphatic site and such)
            # Copy key and remove the atom to move (first character of string)
            key = reactionkey.replace(reactionkey[0], "")
            # Split the remaining key into the current site number, and the next site number
            current = key.split("to")[0]
            molecule.positions.append(current)

            if reactionkey[0] == "D":
                D_hops += 1

            # Do the scrambling
            React.Do_scramble(reactionkey, molecule)
        elif "diss" in reactionkey:
            # Save the atom that dissociates
            diss_atom = reactionkey[0]
            # Save the site the dissociation happened
            diss_position = reactionkey.replace(reactionkey[0], "")
            diss_position = diss_position.replace("diss", "")

            # Do the dissociation
            React.Do_dissociation(diss_atom, molecule)
            break

        # Keep track of how many 'hops' are done by the aliphatic site
        total_hops += 1

        if total_hops % 500000 == 0:
            print(f"MC {str(j_iter)} hops: {str(total_hops)}", flush=True)

    print(f"MC {j_iter} ending", flush=True)
    queue.put(
        [
            diss_atom,
            diss_position,
            time,
            total_hops,
            D_hops,
            molecule.edges,
            molecule.HH_time,
            molecule.HD_time,
            molecule.DD_time,
            j_iter,
        ]
    )


def Do_MC(inputfile, outputfile, cores):
    """Function to perform multiple Monte Carlo simulations"""

    queue = mp.Queue()

    input = Input_reader(inputfile)
    warn_setting = input.handling

    molecule = Molecule(input)

    Check_available_rates(input, molecule, warn_setting)

    if isinstance(input.energy_range, float):
        Energy = [input.energy_range]
    else:
        Energy = np.linspace(
            input.energy_range[0], input.energy_range[1], num=input.energy_range[2] + 1
        )

    dissociation_atoms = {}
    dissociation_times = {}
    dissociation_positions = {}
    N_scramble_hops = {}
    N_D_hops = {}

    # Run multiple MC events
    for value in Energy:
        dissociation_atoms[value] = {"H": 0, "D": 0, "None": 0}
        N_scramble_hops[value] = []
        N_D_hops[value] = []
        dissociation_times[value] = []
        dissociation_positions[value] = {}

        for i in range(len(molecule.edge_numbers)):
            dissociation_positions[value][i] = 0

        if __name__ == "__main__":
            for iter in range(0, input.iterations, cores):
                end = (
                    iter + cores
                    if iter + cores < input.iterations
                    else input.iterations
                )
                run_iterations(iter, end, input, value, molecule, queue, outputfile)

                while True:
                    (
                        diss_atom,
                        diss_position,
                        time,
                        hops,
                        D_hops,
                        end_struct,
                        HH_time,
                        HD_time,
                        DD_time,
                        mc,
                    ) = queue.get()

                    if diss_atom == None:
                        dissociation_atoms[value]["None"] += 1
                    else:
                        logger.info(f"diss_atom={diss_atom}, diss_position={diss_position}, value={value}, time={time}, hopes={hops}, D_hops={D_hops}")
                        dissociation_atoms[value][diss_atom] += 1
                        dissociation_positions[value][diss_position] = dissociation_positions[value].get(diss_position, 0) + 1
                        dissociation_times[value].append(time)
                    N_scramble_hops[value].append(hops)
                    N_D_hops[value].append(D_hops)

                    Data_Output(
                        outputfile,
                        value,
                        diss_atom,
                        diss_position,
                        time,
                        hops,
                        D_hops,
                        HH_time,
                        HD_time,
                        DD_time,
                        mc,
                    )

                    if diss_atom != None:
                        End_Structures_Output(outputfile, value, end_struct, mc)

                    if (
                        len(N_scramble_hops[value]) == iter + cores
                        or len(N_scramble_hops[value]) == input.iterations
                    ):
                        break

    if __name__ == "__main__":
        STD_Output(
            outputfile,
            dissociation_atoms,
            dissociation_positions,
            dissociation_times,
            N_scramble_hops,
            N_D_hops,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="PAHMC",
        description="Process input file and rate definition file."
    )
    parser.add_argument("inputfile", type=str, help="Input YAML file")
    parser.add_argument("cores", type=int, help="Number of parallel processes to run")
    parser.add_argument("-o", "--output", type=str, help="Output file", default=None)
    parser.add_argument("-l", "--log", type=str, help="Log file", default=None)

    args = parser.parse_args()

    sim_name = args.inputfile.split(".")[0]

    outputfile = sim_name + ".out" if args.output is None else args.output
    logfile = sim_name + ".log" if args.log is None else args.log

    logging.basicConfig(
        filename=logfile,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    Do_MC(args.inputfile, outputfile, args.cores)

    logger.info("Simulation done.")
