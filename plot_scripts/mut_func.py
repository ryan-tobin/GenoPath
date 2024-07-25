import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import os
import re
import pandas as pd
import sys
from graphviz import Digraph
from matplotlib.colors import CSS4_COLORS
from PIL import Image, ImageDraw, ImageFont, ImageColor
from math import sqrt
import time
import glob as glob


def parse_csv_files(directory):
    mutation_data = {}
    for file_path in glob.glob(os.path.join(directory, "*.csv")):
        try:
            df = pd.read_csv(
                file_path, header=None, sep=",", na_values="NA", keep_default_na=False
            )
            if df.shape[1] < 3:
                print(
                    f"File {os.path.basename(file_path)} has less than 3 columns. Skipping."
                )
                continue
            for index, row in df.dropna().iterrows():
                positions = []
                raw_positions = str(row[2]).split(";")
                for raw_position in raw_positions:
                    if "-" in raw_position:
                        base, range_part = raw_position.split(":")
                        start_pos, end_pos = range_part.split("-")
                        positions.append(f"{base}:{start_pos}")
                        positions.append(f"{base}:{end_pos}")
                    else:
                        positions.append(raw_position)

                for position in positions:
                    if position:
                        file_identifier = os.path.basename(file_path).split("_")[0]
                        mutation_data.setdefault(position, []).append(
                            (file_identifier, row[0], row[1])
                        )
        except Exception as e:
            print(f"Error processing file {os.path.basename(file_path)}: {e}")
    return mutation_data


def process_driver_file(file_path):
    driver_mutations = pd.read_csv(file_path, sep="\t")
    position_to_gene = {}
    for index, row in driver_mutations.iterrows():
        raw_positions = row["Mutation Summary"].split(" ")[0]
        gene_name = row["Driver Gene"]
        if "-" in raw_positions:
            start, end = raw_positions.split("-")
            base, start_pos = start.rsplit(":", 1)
            end_pos = end
            position_to_gene[f"{base}:{start_pos}"] = gene_name
            position_to_gene[f"{base}:{end_pos}"] = gene_name
        else:
            position_to_gene[raw_positions] = gene_name

    return driver_mutations, position_to_gene


def match_positions(mutation_data, driver_mutations, position_to_gene):
    matched_positions = {}
    for position, gene_name in position_to_gene.items():
        if position in mutation_data:
            matched_positions[position] = mutation_data[position]
            for file_identifier, mutation_type, count in mutation_data[position]:
                print(
                    f"Gene {gene_name} Position {position} has mutation count {count} in file {file_identifier} with mutation {mutation_type} and is a driver mutation."
                )
        else:
            print(
                f"No match found for position: {position} belonging to gene {gene_name}"
            )
    return matched_positions


def mutation_plots(
    target_dir,
    csv_files_dir,
    matched_positions,
    position_to_gene,
    position_to_gene_name,
):
    print("Creating mutation plots...")
    time.sleep(1)
    csv_files = [f for f in os.listdir(csv_files_dir) if f.endswith(".csv")]
    alphabet = "CDEFGHIJKLMNOPQRSTUVWXYZ"
    plots = []

    color_map = {
        "C>A": "blue",
        "C>G": "black",
        "C>T": "red",
        "T>A": "silver",
        "T>C": "lightgreen",
        "T>G": "pink",
    }

    driver_mutation_color = "gold"
    color_map_with_driver = {**color_map, "Driver": driver_mutation_color}

    legend_patches = [
        Patch(color=color, label=label)
        for label, color in color_map_with_driver.items()
    ]

    for index, file in enumerate(csv_files):
        file_path = os.path.join(csv_files_dir, file)
        title_prefix = alphabet[index % len(alphabet)]
        title_branch = file.split("_")[0]

        df = pd.read_csv(
            file_path,
            header=None,
            usecols=[0, 1],
            names=["Substitution", "Variant Count"],
        )
        df["Substitution Type"] = df["Substitution"].str.extract(r"\[(.*?)\]")
        df["Variant Count"] = (
            pd.to_numeric(df["Variant Count"], errors="coerce").fillna(0).astype(int)
        )

        fig, ax = plt.subplots(figsize=(10, 5))
        unique_types = df["Substitution Type"].unique()
        bar_positions = []
        bar_width = 0.4

        current_position = 0
        for i, mutation_type in enumerate(unique_types):
            mutation_df = df[df["Substitution Type"] == mutation_type]
            midpoint = current_position + (len(mutation_df) - 1) / 2
            bar_positions.append(midpoint)
            current_position += len(mutation_df)
            for _, row in mutation_df.iterrows():
                color = color_map[mutation_type]  
                gene_name = None  
                for position, matches in matched_positions.items():
                    for match in matches:
                        filename = match[0].split("_")[0]
                        substitution = match[1]
                        updated_match = f"{filename}:{substitution}"
                        if (
                            filename == title_branch
                            and substitution == row["Substitution"]
                        ):
                            base_position = (
                                f"{position.split(':')[0]}:{position.split(':')[1]}"
                            )
                            gene_name = position_to_gene_name.get(base_position)
                            break
                    if gene_name:
                        break

                if gene_name:
                    ax.bar(
                        row["Substitution"],
                        row["Variant Count"],
                        color="gold",
                        width=bar_width,
                    )
                    ax.text(
                        row["Substitution"],
                        row["Variant Count"] + np.max(df["Variant Count"]) * 0.05,
                        gene_name,
                        ha="center",
                        va="bottom",
                        rotation=90,
                        color="black",
                    )
                else:
                    ax.bar(
                        row["Substitution"],
                        row["Variant Count"],
                        color=color,
                        width=bar_width,
                    )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1)
        ax.spines["bottom"].set_linewidth(1)
        max_count = df["Variant Count"].max()
        ax.yaxis.set_major_locator(plt.MultipleLocator(10))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(5))

        ticks = np.arange(0, max_count + 1, 5)
        ax.set_yticks(ticks)
        ax.set_yticklabels([str(int(x)) for x in ticks], fontweight="bold")
        ax.set_xticks(bar_positions)
        ax.set_xticklabels(df["Substitution Type"].unique(), rotation=90, fontsize=10)

        ax.set_xlabel("Substitution Type", fontsize=12, fontweight="bold")
        ax.set_ylabel("Variant Count", fontsize=12, fontweight="bold")
        ax.set_title(
            f"({title_prefix}) Branch {title_branch}",
            loc="left",
            fontsize=14,
            fontweight="bold",
        )
        ax.legend(handles=legend_patches)

        plt.tight_layout()
        plot_file_name = f"{title_branch}_{title_prefix}_bar_plot.png"
        plt.savefig(os.path.join(target_dir, plot_file_name))
        plt.close()
        plots.append(os.path.join(target_dir, plot_file_name))

        images = [Image.open(plot) for plot in plots]
        num_columns = 5
        num_rows = int(np.ceil(len(images) / num_columns))

        default_image = Image.new("RGB", (0, 0))

        max_widths_per_column = [
            max(
                (images[i::num_columns] if images[i::num_columns] else [default_image]),
                key=lambda img: img.width,
            ).width
            for i in range(num_columns)
        ]
        max_heights_per_row = [
            max(
                [img.size[1] for img in images[i * num_columns : (i + 1) * num_columns]]
                or [default_image.size[1]]
            )
            for i in range(num_rows)
        ]

        total_width = sum(max_widths_per_column)
        total_height = sum(max_heights_per_row)

        new_im = Image.new("RGB", (total_width, total_height), "white")

        y_offset = 0
        for row in range(num_rows):
            x_offset = 0
            for col in range(num_columns):
                index = row * num_columns + col
                if index < len(images):
                    im = images[index]
                    new_im.paste(im, (x_offset, y_offset))
                    x_offset += im.size[0]
                else:
                    x_offset += max_widths_per_column[col]
            y_offset += max_heights_per_row[row]

        combined_plot_file_name = "combined_bar_plots.png"
        new_im.save(os.path.join(target_dir, combined_plot_file_name))
        print("Mutation plots generated.")
        print(
            f"Combined bar chart image saved to {os.path.join(target_dir, combined_plot_file_name)}"
        )
