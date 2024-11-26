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
from PAHMC.Output import (
    STD_Output,
    Data_Output,
    Structure_Output,
    End_Structures_Output,
    Hops_Output,
)


def Parallel_Single_MC(E, max_time, molecule, rates, queue, j_iter, outfilename, debug):
    print(f"MC {j_iter} starting", flush=True)

    specified_rates = list(rates.reactionrates.keys())
    energy = cp.copy(E)
    time = 0
    total_hops = 0
    D_hops = 0
    cross_hops = 0
    diss_atom = None
    diss_position = None
    key_hops = {key: 0 for key in set(molecule.reactionkeys.values()) if "to" in key}

    while time < max_time:
        if debug and j_iter == 0:
            Structure_Output(outfilename, E, j_iter, molecule)

        molecule.possible_reactions = Possible_reactions(molecule, specified_rates)
        reactionkey, dt = Choose_reaction(energy, molecule.possible_reactions, rates)

        if reactionkey == None:
            logger.error("Ran out of reaction rates before dissociation.")
            break

        time += dt
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

        if "to" in reactionkey:
            key = reactionkey.replace(reactionkey[0], "")
            current, next = key.split("to")
            molecule.positions.append(current)

            if reactionkey[0] == "D":
                D_hops += 1

            if next == molecule.cross_links.get(current):
                cross_hops += 1

            key_hops[molecule.reactionkeys.get(reactionkey)] += 1

            React.Do_scramble(reactionkey, molecule)

        elif "diss" in reactionkey:
            diss_atom = reactionkey[0]
            diss_position = reactionkey.replace(reactionkey[0], "")
            diss_position = diss_position.replace("diss", "")
            React.Do_dissociation(diss_atom, molecule)
            break

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
            cross_hops,
            key_hops,
        ]
    )


def worker(iter, input, value, molecule, queue, outputfile, debug):
    Parallel_Single_MC(
        value,
        input.t_max,
        cp.deepcopy(molecule),
        input,
        queue,
        iter,
        outputfile,
        debug,
    )


def process_results(queue, input, value, iter):
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
            cross_hops,
            key_hops,
        ) = queue.get()

        if diss_atom is None:
            dissociation_atoms[value]["None"] += 1
        else:
            dissociation_atoms[value][diss_atom] += 1
            dissociation_positions[value][diss_position] = (
                dissociation_positions[value].get(diss_position, 0) + 1
            )
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
            cross_hops,
        )

        Hops_Output(outputfile, mc, key_hops, value)

        if diss_atom is not None:
            End_Structures_Output(outputfile, value, end_struct, mc)

        if (
            len(N_scramble_hops[value]) == iter + 1
            or len(N_scramble_hops[value]) == input.iterations
        ):
            break


def main(inputfile, outputfile, cores, debug):
    logger.info(f"Running on {cores} cores")
    manager = mp.Manager()
    queue = manager.Queue()
    logger.info(f"Reading data from: {inputfile}")
    input = Input_reader(inputfile)
    warn_setting = input.handling
    molecule = Molecule(input)

    Check_available_rates(input, molecule, warn_setting)

    Energy = input.energy

    global dissociation_atoms, dissociation_times, dissociation_positions, N_scramble_hops, N_D_hops
    dissociation_atoms = {}
    dissociation_times = {}
    dissociation_positions = {}
    N_scramble_hops = {}
    N_D_hops = {}

    dissociation_atoms[Energy] = {"H": 0, "D": 0, "None": 0}
    N_scramble_hops[Energy] = []
    N_D_hops[Energy] = []
    dissociation_times[Energy] = []
    dissociation_positions[Energy] = {}

    for i in molecule.edge_numbers:
        for j in i:
            dissociation_positions[Energy][j] = 0

    pool = mp.Pool(cores)

    logger.info("Starting simulations")
    tasks = [
        (iter, input, Energy, molecule, queue, outputfile, debug)
        for iter in range(input.iterations)
    ]
    pool.starmap(worker, tasks)
    pool.close()
    pool.join()
    logger.info("Finished simulations")

    process_results(queue, input, Energy, input.iterations)

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
        description="Perform Monte Carlo simulation of scrambling and photodissociation reactions on PAHs.",
    )
    parser.add_argument("inputfile", type=str, help="Input YAML file")
    parser.add_argument(
        "-c",
        "--cores",
        type=int,
        help="Number of parallel processes to run",
        default=mp.cpu_count(),
    )
    parser.add_argument("-o", "--output", type=str, help="Output file", default=None)
    parser.add_argument("-l", "--log", type=str, help="Log file", default=None)
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debugging")

    args = parser.parse_args()

    sim_name = args.inputfile.split(".")[0]

    outputfile = sim_name + ".out" if args.output is None else args.output
    logfile = sim_name + ".log" if args.log is None else args.log

    debug = args.debug

    logging.basicConfig(
        filename=logfile,
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    main(args.inputfile, outputfile, args.cores, args.debug)

    logger.info("Simulation done.")
