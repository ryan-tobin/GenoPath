import matplotlib.pyplot as plt
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


def is_color_dark(rgb):
    r, g, b = int(rgb[1:3], 16), int(rgb[3:5], 16), int(rgb[5:7], 16)
    brightness = 0.299 * r + 0.587 * g + 0.114 * b
    return brightness < 128


def generate_color_map(signatures):
    dark_colors = [
        color
        for color, hex in CSS4_COLORS.items()
        if is_color_dark(hex)
        and "light" not in color
        and "white" not in color
        and "gray" not in color
        and color != "black"
    ]

    color_map = {
        signature: dark_colors[i % len(dark_colors)]
        for i, signature in enumerate(signatures)
    }
    return color_map


def read_signatures_from_control_file(control_file_path):
    with open(control_file_path, "r") as file:
        for line in file:
            parts = line.strip().split("\t") 
            if parts[0] == "Signature List":
                signatures = parts[1].split(",")
                return signatures  

    return []


def signature_pie_plots(target_dir, signature_file_dir, control_file_path):
    print("Creating pie plots for signature distribution...")
    time.sleep(1)
    signatures = read_signatures_from_control_file(control_file_path)
    print(f"Signature names: {signatures}")
    color_map = generate_color_map(signatures)
    print(f"Color map: {color_map}")

    df = pd.read_csv(signature_file_dir, sep="\t")

    relevant_columns = [col for col in df.columns if "Relative activity" in col]
    if relevant_columns:
        activity_col = relevant_columns[0]
        print(f"The relevant column is named: {activity_col}")
    else:
        print("Relative activity column not found.")
        return

    branches = df["Branch"].unique()
    print(f"Unique Branches Are: {branches}")
    for branch in branches:
        branch_data = df[df["Branch"] == branch][["Signature", activity_col]]
        print(f"Relevant branch data: {branch_data}")
        sizes = branch_data[activity_col].values
        print(f"Sizes: {sizes}")
        labels = branch_data["Signature"].values
        print(f"Labels: {labels}")
        colors = [color_map[label] for label in labels if label in color_map]
        print(f"Colors: {colors}")

        fig, ax = plt.subplots()
        wedges, texts = ax.pie(sizes, colors=colors, startangle=90)

        for i, p in enumerate(wedges):
            ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
            x = np.cos(np.deg2rad(ang))
            y = np.sin(np.deg2rad(ang))
            label_x = x * 0.72 * p.r
            label_y = y * 0.8 * p.r
            ax.text(
                label_x,
                label_y,
                labels[i],
                ha="center",
                va="center",
                fontsize=20,
                fontname="Arial",
            )

        ax.axis("equal")

        output_filename = os.path.join(target_dir, f"{branch}_signature_pie.png")
        plt.savefig(output_filename, transparent=True)
        plt.close()
        print("Pie plots created.")


def euclidean_distance(point1, point2):
    return sqrt((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)
