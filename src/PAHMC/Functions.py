# Import important libraries
import numpy as np
import logging
import sys

logger = logging.getLogger(__name__)


def Remove_missing_rate(reactions, specified_rates):
    return [key for key in reactions if key in specified_rates]


def Find_aliphatic_site(edges):
    """Finds the index of the aliphatic site(s) from the initial edge structure"""

    index = []
    al_deut = False  # is deuterium on the aliphatic site

    for i in range(len(edges)):
        if len(edges[i]) == 1:
            if edges[i][0] in (
                "HH",
                "HD",
                "DD",
                "DH",
            ):  # second brackets: [0] needed due to lists in list object
                index.append([i])
                if "D" in edges[i][0]:
                    al_deut = True
        elif len(edges[i]) in (2, 3, 4):
            for j in range(len(edges[i])):
                if edges[i][j] in ("HH", "HD", "DD", "DH"):
                    index.append([i, j])
                    if "D" in edges[i][j]:
                        al_deut = True
        else:
            logger.error("Invalid edge found.")
            sys.exit(2)

    if len(index) == 0:
        logger.error("No aliphatic site detected on edge of PAH, check input.")
        sys.exit(2)

    elif len(index) > 1:
        logger.error("Currently only one aliphatic site is implemented.")
        sys.exit(2)

    return index, al_deut


def Possible_reactions(molecule, specified_rates):
    """Determine the possible reactions that can happen from the site of the aliphatic group."""

    reactions = []
    deut = 0

    if len(molecule.index[0]) == 1:
        i = molecule.index[0][0]
        j = 0

    elif len(molecule.index[0]) == 2:
        i, j = molecule.index[0]

    if molecule.al_place == "e":
        # count how many deuterium in the aliphatic site
        deut = molecule.edges[i][j].count("D")

        current = molecule.edge_numbers[i][j]

        if j == 0:
            prev = molecule.link_numbers[i - 1][-1]
            if len(molecule.edge_numbers[i]) > 1:
                next = molecule.edge_numbers[i][j + 1]
            else:
                next = molecule.link_numbers[i][0]
        elif j != 0 and j < len(molecule.edge_numbers[i]) - 1:
            prev = molecule.edge_numbers[i][j - 1]
            next = molecule.edge_numbers[i][j + 1]
        elif j == len(molecule.edge_numbers[i]) - 1:
            prev = molecule.edge_numbers[i][j - 1]
            next = molecule.link_numbers[i][0]

    if molecule.al_place == "l":
        deut = molecule.links[i][j].count("D")
        current = molecule.link_numbers[i][j]

        if j == 0:
            prev = molecule.edge_numbers[i][-1]
            if len(molecule.link_numbers[i]) > 1:
                next = molecule.link_numbers[i][j + 1]
            else:
                if i < len(molecule.edge_numbers) - 1:
                    next = molecule.edge_numbers[i + 1][0]
                else:
                    next = molecule.edge_numbers[0][0]
        elif j != 0 and j < len(molecule.link_numbers[i]) - 1:
            prev = molecule.link_numbers[i][j - 1]
            next = molecule.link_numbers[i][j + 1]
        elif j == len(molecule.edge_numbers[i]) - 1:
            prev = molecule.link_numbers[i][j - 1]
            if i < len(molecule.edge_numbers) - 1:
                next = molecule.edge_numbers[i + 1][0]
            else:
                next = molecule.edge_numbers[0][0]

    if deut == 0:
        # Scrambling
        reactions.append(f"H{current}to{prev}")
        reactions.append(f"H{current}to{next}")
        # Dissociation
        reactions.append(f"H{current}diss")
        # if on edge two hydrogens that can move, so append again
        if molecule.al_place == "e":
            # Scrambling
            reactions.append(f"H{current}to{prev}")
            reactions.append(f"H{current}to{next}")
            # Dissociation
            reactions.append(f"H{current}diss")

    elif deut == 1:
        # Scrambling
        reactions.append(f"D{current}to{prev}")
        reactions.append(f"D{current}to{next}")
        # Dissociation
        reactions.append(f"D{current}diss")
        if molecule.al_place == "e":
            # Scrambling
            reactions.append(f"H{current}to{prev}")
            reactions.append(f"H{current}to{next}")
            # Dissociation
            reactions.append(f"H{current}diss")

    elif deut == 2:
        # Scrambling
        reactions.append(f"D{current}to{prev}")
        reactions.append(f"D{current}to{next}")
        # Dissociation
        reactions.append(f"D{current}diss")
        # if on edge two hydrogens that can move, so append again
        if molecule.al_place == "e":
            # Scrambling
            reactions.append(f"D{current}to{prev}")
            reactions.append(f"D{current}to{next}")
            # Dissociation
            reactions.append(f"D{current}diss")

    reactions = Remove_missing_rate(reactions, specified_rates)

    return reactions


def Probability(dt, weight, rate):
    """Determine the reaction probabilities from reaction rates through P = 1-exp(-w*k*dt)"""
    P = 1 - np.exp(-dt * weight * rate)

    return P


def Weights_degeneracy(possible_reac):
    """Determines the weights according to the degeneracies of each reaction type."""
    degeneracy = {}
    nreac = 0

    for item in possible_reac:
        if item not in degeneracy:
            x = possible_reac.count(item)
            degeneracy[item] = x
            nreac += 1

    return degeneracy, nreac


def Choose_reaction(E, possible_reactions, rates):
    """Chooses a random reaction from the possible reactions according to the probabilities"""

    np.random.seed()

    degeneracy, nreactions = Weights_degeneracy(possible_reactions)

    # split reaction keys and their weights (degeneracies)
    reactionkeys = list(degeneracy.keys())
    weights = list(degeneracy.values())

    # Prepare empty arrays for reaction rates and reaction probabilities
    r_probabilities = np.zeros(nreactions)
    r_rate = np.zeros(nreactions)

    for i, r_key in enumerate(reactionkeys):

        # Determine the closes match for the rate for a certain energy:
        rate_idx = np.argmin((np.abs(rates.reactionrates[r_key][0, :] - E)))

        # Save the rate
        r_rate[i] = rates.reactionrates[r_key][1, rate_idx]

    # Make sure that if all rates are 0 the function quits (as there is no reaction able to be chosen)
    if np.sum(r_rate) == 0:
        reaction = None
        dt = None
        return reaction, dt

    # Determine size of timestep
    dt = 1 / np.sum(r_rate)

    # Determine probabilities of the reactions
    for j in range(nreactions):
        r_probabilities[j] = Probability(dt, weights[j], r_rate[j])

    # Normalize the probabilities (otherwise numpy.random.choice gives an error)
    norm_r_prob = r_probabilities / np.sum(r_probabilities)

    # Choose a reaction according to the probabilities
    reaction = np.random.choice(reactionkeys, p=norm_r_prob)

    return reaction, dt
