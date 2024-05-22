"""
Description:    Visual summary of the performance results. 
Dependencies:   matplotlib==3.8.2
Usage:          python3 visual_summary_performancy.py -i path/to/performance/results.tsv -o path/to/out/dir
"""

import argparse
import math
from collections import defaultdict

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def cli() -> argparse.Namespace:
    """
    Command line interface for this script.

    :return: Namespace of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to input TSV file.")
    parser.add_argument("-o", type=str, required=True, help="Path to output directory.")
    return parser.parse_args()


def main() -> None:
    """
    Driver function.
    """
    args = cli()

    fontsize = 18

    # Parse performance results.
    data = defaultdict(list)
    with open(args.i, "r") as file_open:
        file_open.readline()  # Skip header.
        for line in file_open:
            # ind, num_heavy_atoms, num_bonds, style, look, runtime, file_size = line.strip().split("\t")
            ind, *rest = line.strip().split("\t")
            data[ind].append(rest)

    # Histogram of number of heavy atoms.
    num_heavy_atoms = [int(data[ind][0][0]) for ind in data]
    num_bonds = [int(data[ind][0][1]) for ind in data]
    bins = max(max(max(num_heavy_atoms), max(num_bonds)), 60)
    plt.hist(
        num_heavy_atoms,
        bins=bins,
        range=(0, bins),
        edgecolor="red",
        alpha=0.7,
        label="Number of heavy atoms",
        histtype="step",
    )
    plt.hist(
        num_bonds,
        bins=bins,
        range=(0, bins),
        edgecolor="blue",
        alpha=0.7,
        label="Number of bonds",
        histtype="step",
    )
    plt.ylim(0, 300)
    plt.grid(axis="y", linestyle="--", color="black", alpha=0.3)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.xlabel("Count", fontsize=fontsize)
    plt.ylabel("Frequency", fontsize=fontsize)
    legend = plt.legend(loc="upper left", fontsize=14)
    legend.get_frame().set_alpha(0.25)
    plt.savefig(
        f"{args.o}/atom_and_bound_counts.png", dpi=300, bbox_inches="tight", transparent=True
    )
    plt.clf()

    # Plot average runtime per number of heavy atoms for different styles for cartoon look.
    def plot_runtime(style: str, look: str, label: str, color: str, linestyle: str):
        spacefilling = defaultdict(list)
        for ind in data:
            for item in data[ind]:
                if item[2] == style and item[3] == look:
                    spacefilling[int(item[0])].append(float(item[4]))

        num_bins = 50
        bins = [i for i in range(1, num_bins + 1, 1)]
        bin_heights = [0 for _ in range(num_bins)]
        bin_yerrs = [0 for _ in range(num_bins)]
        for i in range(0, num_bins):
            bin = i + 1
            values = spacefilling[bin]
            if len(values) == 0:
                continue
            avg = sum(values) / len(values)
            std = math.sqrt(sum([(value - avg) ** 2 for value in values]) / len(values))
            bin_heights[i] = avg
            bin_yerrs[i] = std

        # Plot avg runtime vs. number of heavy atoms for every style for Cartoon.
        non_zero = [i for i in range(len(bin_heights)) if bin_heights[i] != 0]
        bins = [bins[i] for i in non_zero]
        bin_heights = [bin_heights[i] for i in non_zero]
        alpha = 0.7 if linestyle == "-" else 1.0
        plt.plot(bins, bin_heights, color=color, linestyle=linestyle, alpha=alpha)
        plt.fill_between(
            bins,
            [bin_heights[i] - bin_yerrs[i] for i in range(len(bin_heights))],
            [bin_heights[i] + bin_yerrs[i] for i in range(len(bin_heights))],
            color=color,
            alpha=0.1,
        )

    plot_runtime("SPACEFILLING", "CARTOON", "Space-filling", "red", "-")
    plot_runtime("BALL_AND_STICK", "CARTOON", "Ball-and-stick", "blue", "-")
    plot_runtime("TUBE", "CARTOON", "Tube", "green", "-")
    plot_runtime("WIREFRAME", "CARTOON", "Wireframe", "orange", "-")

    # Plot average runtime per number of heavy atoms for different styles for glossy look.
    plot_runtime("SPACEFILLING", "GLOSSY", "Space-filling", "red", ":")
    plot_runtime("BALL_AND_STICK", "GLOSSY", "Ball-and-stick", "blue", ":")
    plot_runtime("TUBE", "GLOSSY", "Tube", "green", ":")
    plot_runtime("WIREFRAME", "GLOSSY", "Wireframe", "orange", ":")

    plt.grid(axis="y", linestyle="--", color="black", alpha=0.3)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylim(0, 1200)
    plt.xlim(0, 50)
    plt.xlabel("Number of heavy atoms", fontsize=fontsize)
    plt.ylabel("Runtime (milliseconds)", fontsize=fontsize)

    blue_patch = mpatches.Patch(color="red", label="Space-filling")
    green_patch = mpatches.Patch(color="blue", label="Ball-and-stick")
    red_patch = mpatches.Patch(color="green", label="Tube")
    yellow_patch = mpatches.Patch(color="orange", label="Wireframe")
    dotted_line = mlines.Line2D([], [], linestyle="-", color="black", label="Cartoon")
    solid_line = mlines.Line2D([], [], linestyle=":", color="black", label="Glossy")
    legend_handles = [blue_patch, green_patch, red_patch, yellow_patch, dotted_line, solid_line]
    legend = plt.legend(loc="upper left", fontsize=14, handles=legend_handles)
    legend.get_frame().set_alpha(0.25)

    plt.savefig(f"{args.o}/speed_per_molecule.png", dpi=300, bbox_inches="tight", transparent=True)
    plt.clf()

    # Plot average file size per number of heavy atoms for different styles for cartoon look.
    def plot_runtime(style: str, look: str, label: str, color: str, linestyle: str):
        spacefilling = defaultdict(list)
        for ind in data:
            for item in data[ind]:
                if item[2] == style and item[3] == look:
                    spacefilling[int(item[0])].append(float(item[5]))

        num_bins = 50
        bins = [i for i in range(1, num_bins + 1, 1)]
        bin_heights = [0 for _ in range(num_bins)]
        bin_yerrs = [0 for _ in range(num_bins)]
        for i in range(0, num_bins):
            bin = i + 1
            values = spacefilling[bin]
            if len(values) == 0:
                continue
            avg = sum(values) / len(values)
            std = math.sqrt(sum([(value - avg) ** 2 for value in values]) / len(values))
            bin_heights[i] = avg
            bin_yerrs[i] = std

        # Plot avg runtime vs. number of heavy atoms for every style for Cartoon.
        non_zero = [i for i in range(len(bin_heights)) if bin_heights[i] != 0]
        bins = [bins[i] for i in non_zero]
        bin_heights = [bin_heights[i] for i in non_zero]
        alpha = 0.7 if linestyle == "-" else 1.0
        plt.plot(bins, bin_heights, color=color, linestyle=linestyle, alpha=alpha)
        plt.fill_between(
            bins,
            [bin_heights[i] - bin_yerrs[i] for i in range(len(bin_heights))],
            [bin_heights[i] + bin_yerrs[i] for i in range(len(bin_heights))],
            color=color,
            alpha=0.1,
        )

    plot_runtime("SPACEFILLING", "CARTOON", "Space-filling", "red", "-")
    plot_runtime("BALL_AND_STICK", "CARTOON", "Ball-and-stick", "blue", "-")
    plot_runtime("TUBE", "CARTOON", "Tube", "green", "-")
    plot_runtime("WIREFRAME", "CARTOON", "Wireframe", "orange", "-")

    # Plot average file size per number of heavy atoms for different styles for glossy look.
    plot_runtime("SPACEFILLING", "GLOSSY", "Space-filling", "red", ":")
    plot_runtime("BALL_AND_STICK", "GLOSSY", "Ball-and-stick", "blue", ":")
    plot_runtime("TUBE", "GLOSSY", "Tube", "green", ":")
    plot_runtime("WIREFRAME", "GLOSSY", "Wireframe", "orange", ":")

    plt.grid(axis="y", linestyle="--", color="black", alpha=0.3)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylim(0, 150)
    plt.xlim(0, 50)
    plt.xlabel("Number of heavy atoms", fontsize=fontsize)
    plt.ylabel("File size (kilobytes)", fontsize=fontsize)

    blue_patch = mpatches.Patch(color="red", label="Space-filling")
    green_patch = mpatches.Patch(color="blue", label="Ball-and-stick")
    red_patch = mpatches.Patch(color="green", label="Tube")
    yellow_patch = mpatches.Patch(color="orange", label="Wireframe")
    dotted_line = mlines.Line2D([], [], linestyle="-", color="black", label="Cartoon")
    solid_line = mlines.Line2D([], [], linestyle=":", color="black", label="Glossy")
    legend_handles = [blue_patch, green_patch, red_patch, yellow_patch, dotted_line, solid_line]
    legend = plt.legend(loc="upper left", fontsize=14, handles=legend_handles)
    legend.get_frame().set_alpha(0.25)

    plt.savefig(
        f"{args.o}/file_size_per_molecule.png", dpi=300, bbox_inches="tight", transparent=True
    )
    plt.clf()

    exit(0)


if __name__ == "__main__":
    main()
