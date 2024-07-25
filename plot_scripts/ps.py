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

from plot_scripts.mut_func import mutation_plots
from plot_scripts.pie_func import (
    signature_pie_plots,
    read_signatures_from_control_file,
    generate_color_map,
)
from plot_scripts.tree_func import plot_tree, plot_tree_with_driver


def combine(
    target_dir,
    tree_output,
    tree_output_driver,
    bar_plots_output_dir,
    input_file_base_name,
    control_file,
):
    tree_image = Image.open(tree_output)
    tree_image_driver = Image.open(tree_output_driver)
    signatures = read_signatures_from_control_file(control_file)
    color_map = generate_color_map(signatures)
    print(f"Color map: {color_map}")
    scale_factor = 3
    scale_factor_font = 2
    scale_factor_tree = 1.25
    scale_factor_driver = 1.15
    new_width = int(tree_image.width * scale_factor_tree)
    new_height = int(tree_image.height * scale_factor_tree)
    new_width_driver = int(tree_image_driver.width * scale_factor_driver)
    new_height_driver = int(tree_image_driver.height * scale_factor_driver)
    resized_tree_image = tree_image.resize(
        (new_width, new_height), Image.Resampling.LANCZOS
    )
    resized_tree_image_driver = tree_image_driver.resize(
        (new_width_driver, new_height_driver), Image.Resampling.LANCZOS
    )
    base_width = 250
    base_height_per_item = 20
    font_size = 12
    title_height = 30 * scale_factor_font
    legend_width = base_width * scale_factor
    legend_height = (
        title_height + (base_height_per_item *
                        len(color_map) + 10) * scale_factor
    )
    new_font_size = font_size * scale_factor_font
    legend_image = Image.new(
        "RGBA", (legend_width, legend_height), (255, 255, 255, 0))
    draw = ImageDraw.Draw(legend_image)
    font = ImageFont.load_default(new_font_size)
    title_text = "Signature Color Map"
    title_x = 0
    title_y = 0
    draw.text((title_x, title_y), title_text, fill="black", font=font)
    for i, (label, color) in enumerate(color_map.items()):
        rect_start_x = 10 * scale_factor
        rect_start_y = (
            title_height + (i * base_height_per_item + 5) * scale_factor - 10
        )
        rect_end_x = 30 * scale_factor
        rect_end_y = title_height + \
            (i * base_height_per_item + 15) * scale_factor - 10
        text_start_x = 35 * scale_factor
        text_start_y = rect_start_y

        rgba_color = ImageColor.getcolor(color, "RGBA")
        draw.rectangle(
            [rect_start_x, rect_start_y, rect_end_x, rect_end_y], fill=rgba_color
        )
        draw.text((text_start_x, text_start_y), label, fill="black", font=font)

    bar_plots_image = Image.open(bar_plots_output_dir)

    total_width = max(tree_image.width + 500,
                      bar_plots_image.width, legend_image.width)
    total_height = (
        tree_image.height + 50 + legend_height + bar_plots_image.height
    )

    combined_image = Image.new("RGB", (total_width, total_height), "white")

    combined_image.paste(resized_tree_image, (10, 20))
    combined_image.paste(resized_tree_image_driver,
                         (10 + resized_tree_image.width, 20))
    margin = 20
    x_coordinate = resized_tree_image.width + margin
    y_coordinate = 10

    combined_image.paste(
        legend_image, (x_coordinate, y_coordinate), legend_image)

    combined_image.paste(
        bar_plots_image, (0, tree_image.height + legend_image.height + 150))

    output_file_name = f"Phylo_bar_final.png"
    combined_image.save(os.path.join(target_dir, output_file_name))

    print("Combined image with tree, legend, and bar plots saved.")


def run(
    target_dir,
    csv_files_dir,
    signature_files_dir,
    control_file,
    summary_file_path,
    input_file_base_name,
    tree_output,
    tree_output_driver,
    bar_plots_output_dir,
    matched_positions,
    driver_data,
    position_to_gene,
    position_to_gene_name,
):
    print("Plot creation starting...")
    time.sleep(2)
    mutation_plots(
        target_dir,
        csv_files_dir,
        matched_positions,
        position_to_gene,
        position_to_gene_name,
    )
    signature_pie_plots(target_dir, signature_files_dir, control_file)
    plot_tree(target_dir, summary_file_path, input_file_base_name)
    plot_tree_with_driver(
        target_dir,
        summary_file_path,
        input_file_base_name,
        matched_positions,
        driver_data,
    )
    combine(
        target_dir,
        tree_output,
        tree_output_driver,
        bar_plots_output_dir,
        input_file_base_name,
        control_file,
    )
    print("Plot creation complete.")
    time.sleep(2)
