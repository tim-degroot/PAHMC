import logging
import sys


logger = logging.getLogger(__name__)


def Determine_needed_rates(deuterium, edge_numbering, link_numbering, cross_links):
    # deuterium variable is a boolean for if deuterium is present or not
    all_reactions = []
    length = len(
        edge_numbering
    )  # both edges and links are same length (checked in input function)

    for i in range(length):
        for j in range(len(edge_numbering[i]) - 1):
            # Append reaction keys for the scrambling on the individual edges
            all_reactions.append(
                "H" + edge_numbering[i][j] + "to" + edge_numbering[i][j + 1]
            )
            all_reactions.append(
                "H" + edge_numbering[i][j + 1] + "to" + edge_numbering[i][j]
            )
            if deuterium:
                all_reactions.append(
                    "D" + edge_numbering[i][j] + "to" + edge_numbering[i][j + 1]
                )
                all_reactions.append(
                    "D" + edge_numbering[i][j + 1] + "to" + edge_numbering[i][j]
                )

        # Append the reaction keys for individual edge to links and back
        all_reactions.append("H" + edge_numbering[i][-1] + "to" + link_numbering[i][0])
        all_reactions.append("H" + link_numbering[i][0] + "to" + edge_numbering[i][-1])
        if deuterium:
            all_reactions.append(
                "D" + edge_numbering[i][-1] + "to" + link_numbering[i][0]
            )
            all_reactions.append(
                "D" + link_numbering[i][0] + "to" + edge_numbering[i][-1]
            )
        # take last element of link for moving anticlockwise (implemented in case of bay regions)
        k = len(link_numbering[i - 1])
        all_reactions.append(
            "H" + edge_numbering[i][0] + "to" + link_numbering[i - 1][k - 1]
        )
        all_reactions.append(
            "H" + link_numbering[i - 1][k - 1] + "to" + edge_numbering[i][0]
        )
        if deuterium:
            all_reactions.append(
                "D" + edge_numbering[i][0] + "to" + link_numbering[i - 1][k - 1]
            )
            all_reactions.append(
                "D" + link_numbering[i - 1][k - 1] + "to" + edge_numbering[i][0]
            )

        for l in range(len(edge_numbering[i])):
            all_reactions.append("H" + edge_numbering[i][l] + "diss")
            if deuterium:
                all_reactions.append("D" + edge_numbering[i][l] + "diss")

        for m in range(len(link_numbering[i])):
            all_reactions.append("H" + link_numbering[i][m] + "diss")
            if deuterium:
                all_reactions.append("D" + link_numbering[i][m] + "diss")

    for a, b in cross_links.items():
        all_reactions.append(f"D{a}to{b}")
        all_reactions.append(f"H{a}to{b}")

    return all_reactions


def Check_available_rates(rates, molecule, warn_setting):
    """Checks if all rates are specified for every available edge, and acts according to setting"""

    needed_rates = Determine_needed_rates(
        molecule.Deuterium,
        molecule.edge_numbers,
        molecule.link_numbers,
        molecule.cross_links,
    )

    for key in needed_rates:
        if key not in rates.reactionrates.keys():
            if warn_setting in ("w", "W"):
                logger.warning(
                    f"No rates found for {key}, continuing without rates for {key}."
                )
            elif warn_setting in ("q", "Q"):
                logger.error(f"No rates found for {key}, stopping run now.")
                sys.exit(2)
