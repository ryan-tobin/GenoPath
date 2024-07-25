from PIL import Image, ImageDraw, ImageFont
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
import glob
import math


def debug_print_tree_lines(summary_file_path):
    tree_structure_started = False
    print("Contents under #Tree:")
    with open(summary_file_path, "r") as file:
        for line in file:
            if line.startswith("#Tree"):
                tree_structure_started = True
                continue
            if tree_structure_started:
                if line.strip():
                    print(line.strip())
                else:
                    break


def plot_tree(
    target_dir, summary_file_path, input_file_base_name, ranksep=2, nodesep=2
):
    print("Creating Phylogenetic Tree...")
    time.sleep(1)
    find_and_append_root_node(summary_file_path)

    df = pd.read_csv(summary_file_path, sep="\t", comment="#")
    dot = Digraph(comment="Phylogenetic Tree")
    nodes = {}

    if "Normal" not in df["New branch ID"].values:
        dot.node("Normal", "Normal")

    for index, row in df.iterrows():
        node_name = row["New branch ID"]
        mutation_count = row["Number of mutations"]
        if pd.notnull(mutation_count):
            node_label = f"{node_name}\n({int(mutation_count)})"
            dot.node(node_name, node_label)
            nodes[node_name] = int(mutation_count)
            print(f"Node added: {node_name} with mutation count {mutation_count}")
        else:
            print(f"Node {node_name} skipped due to NaN mutation count")

    with open(summary_file_path, "r") as file:
        tree_structure_started = False
        for line in file:
            if line.startswith("#Tree"):
                tree_structure_started = True
                continue
            if tree_structure_started and "->" in line:
                parent, daughter = line.strip().split("->")
                if parent.strip() in nodes and daughter.strip() in nodes:
                    dot.edge(parent.strip(), daughter.strip())
                    print(
                        f"Established relationship: {parent.strip()} -> {daughter.strip()}"
                    )
                elif parent.strip() == "Normal":
                    dot.edge(parent.strip(), daughter.strip())
                    print(
                        f"Established relationship: {parent.strip()} -> {daughter.strip()}"
                    )

    dot.attr(ranksep=str(ranksep), nodesep=str(nodesep))
    graph_title = "(A) Phylogenetic Tree"
    dot.attr("graph", label=graph_title, labelloc="t", fontsize="20", fontcolor="black")
    dot_file_path = os.path.join(target_dir, "Phylogenetic_Tree")
    print(f"Rendering DOT file: {dot_file_path}")
    dot.render(dot_file_path, view=False, format="dot")

    output_file_name = "Phylogenetic_Tree"
    output_file_path = os.path.join(target_dir, output_file_name)
    print(f"Rendering PNG file: {output_file_path}")
    dot.render(output_file_path, view=False, format="png")

    dot_file_path = os.path.join(target_dir, "Phylogenetic_Tree.dot")
    txt_file_path = os.path.join(target_dir, "Phylogenetic_Tree.txt")
    edges_txt_path = os.path.join(target_dir, "Phylogenetic_Tree_Edges.txt")

    save_edge_info_to_txt(dot_file_path, edges_txt_path)

    edge_info = parse_edge_info(edges_txt_path)
    print(f"Edge information: {edge_info}")
    output_file_path = output_file_path
    mark_midpoints_on_png(output_file_path + ".png", edge_info, target_dir)
    print("Midpoints marking completed.")
    print("Phylogenetic Tree created.")


def plot_tree_with_driver(
    target_dir,
    summary_file_path,
    input_file_base_name,
    matched_positions,
    driver_data,
    ranksep=2,
    nodesep=2,
):
    print("Creating Phylogenetic Tree...")
    time.sleep(1)
    find_and_append_root_node(summary_file_path)

    df = pd.read_csv(summary_file_path, sep="\t", comment="#")
    dot = Digraph(comment="Phylogenetic Tree")
    nodes = {}

    if "Normal" not in df["New branch ID"].values:
        dot.node("Normal", "Normal")

    for index, row in df.iterrows():
        node_name = row["New branch ID"]
        mutation_count = row["Number of mutations"]
        if pd.notnull(mutation_count):
            node_label = f"{node_name}\n({int(mutation_count)})"
            dot.node(node_name, node_label)
            nodes[node_name] = int(mutation_count)
            print(f"Node added: {node_name} with mutation count {mutation_count}")
        else:
            print(f"Node {node_name} skipped due to NaN mutation count")

    with open(summary_file_path, "r") as file:
        tree_structure_started = False
        for line in file:
            if line.startswith("#Tree"):
                tree_structure_started = True
                continue
            if tree_structure_started and "->" in line:
                parent, daughter = line.strip().split("->")
                if parent.strip() in nodes and daughter.strip() in nodes:
                    dot.edge(parent.strip(), daughter.strip())
                    print(
                        f"Established relationship: {parent.strip()} -> {daughter.strip()}"
                    )
                elif parent.strip() == "Normal":
                    dot.edge(parent.strip(), daughter.strip())
                    print(
                        f"Established relationship: {parent.strip()} -> {daughter.strip()}"
                    )

    dot.attr(ranksep=str(ranksep), nodesep=str(nodesep))
    graph_title = "(B) Phylogenetic Tree with Driver Mutations"
    dot.attr("graph", label=graph_title, labelloc="t", fontsize="20", fontcolor="black")
    dot_file_path = os.path.join(target_dir, "Phylogenetic_Tree")
    print(f"Rendering DOT file: {dot_file_path}")
    dot.render(dot_file_path, view=False, format="dot")

    output_file_name = "Phylogenetic_Tree"
    output_file_path = os.path.join(target_dir, output_file_name)
    print(f"Rendering PNG file: {output_file_path}")
    dot.render(output_file_path, view=False, format="png")

    dot_file_path = os.path.join(target_dir, "Phylogenetic_Tree.dot")
    txt_file_path = os.path.join(target_dir, "Phylogenetic_Tree.txt")
    edges_txt_path = os.path.join(target_dir, "Phylogenetic_Tree_Edges.txt")

    save_edge_info_to_txt(dot_file_path, edges_txt_path)

    edge_info = parse_edge_info(edges_txt_path)
    print(f"Edge information: {edge_info}")
    output_file_path = output_file_path
    mark_midpoints_on_png_driver(
        output_file_path + ".png", edge_info, target_dir, matched_positions, driver_data
    )
    print("Midpoints marking completed.")
    print("Phylogenetic Tree created.")


def find_and_append_root_node(summary_file_path):
    parent_nodes = set()
    child_nodes = set()
    normal_relationship_exists = False

    with open(summary_file_path, "r") as file:
        for line in file:
            if "->" in line:
                parent, child = line.strip().split("->")
                parent = parent.strip()
                child = child.strip()
                if parent == "Normal":
                    normal_relationship_exists = True
                parent_nodes.add(parent)
                child_nodes.add(child)

    if not normal_relationship_exists:
        root_nodes = parent_nodes - child_nodes

        root_nodes.discard("Normal")

        for root_node in root_nodes:
            append_normal_relationship(summary_file_path, root_node)
            print(
                f"[INFO] Added 'Normal->{root_node}' to ensure 'Normal' is the ultimate root."
            )
    else:
        print(
            "[INFO] A 'Normal' relationship already exists in the file. No new relationship added."
        )


def append_normal_relationship(summary_file_path, root_node):
    with open(summary_file_path, "a") as file:
        file.write(f"Normal->{root_node}\n")
    print(f"[INFO] Added new root relationship to summary file: 'Normal->{root_node}'")


def dot_to_text(dot_path, txt_path):
    try:
        with open(dot_path, "r") as dot_file, open(txt_path, "w") as txt_file:
            for line in dot_file:
                txt_file.write(line)
        print(f"Conversion successful: {dot_path} to {txt_path}")
    except Exception as e:
        print(f"[ERROR] Conversion failed: {e}")


def save_edge_info_to_txt(dot_file_path, edges_txt_path):
    with open(dot_file_path, "r") as dot_file, open(edges_txt_path, "w") as edges_file:
        for line in dot_file:
            if "->" in line:
                edges_file.write(line)
    print(f"[INFO] Edge information saved to: {edges_txt_path}")


def parse_edge_info(edges_txt_path):
    edge_info = {}
    edge_pattern = re.compile(r'(\S+)\s*->\s*(\S+)\s*\[pos="([^"]+)"\];')

    with open(edges_txt_path, "r") as edges_file:
        for line in edges_file:
            match = edge_pattern.search(line)
            if match:
                parent, child, pos_data = match.groups()
                control_points = [
                    (float(x), float(y))
                    for x, y in (
                        pt.split(",") for pt in pos_data.replace("e,", "").split()
                    )
                ]
                edge_info[(parent, child)] = control_points
                print(f"Control points for edge {parent} -> {child}: {control_points}")
    return edge_info


def mark_midpoints_on_png(png_path, edge_info, target_dir, graphviz_dpi=96.0):
    try:
        img = Image.open(png_path)
    except FileNotFoundError:
        print(f"[ERROR] PNG file not found: {png_path}")
        return None

    for (parent, child), control_points in edge_info.items():
        child_prefix = child.split("_")[0]
        pie_chart_path = os.path.join(target_dir, f"{child_prefix}_signature_pie.png")

        pie_chart_image = load_and_resize_image(pie_chart_path, (210, 180))
        if pie_chart_image:
            if len(control_points) >= 4:
                midpoint_graphviz = (
                    (control_points[1][0] + control_points[4][0]) / 2,
                    (control_points[1][1] + control_points[4][1]) / 2,
                )
                midpoint_png = graphviz_to_png_coordinates(
                    midpoint_graphviz, img.height, dpi=graphviz_dpi
                )

                top_left = (
                    int(midpoint_png[0] - pie_chart_image.width / 2),
                    int(midpoint_png[1] - pie_chart_image.height / 2),
                )
                img.paste(pie_chart_image, top_left, pie_chart_image)

    modified_png_path = png_path.replace(".png", "_complete.png")
    img.save(modified_png_path)
    print(f"Specific midpoints with pie charts drawn and saved to {modified_png_path}")


def mark_midpoints_on_png_driver(
    png_path, edge_info, target_dir, matched_positions, driver_data, graphviz_dpi=96.0
):
    try:
        img = Image.open(png_path)
    except FileNotFoundError:
        print(f"[ERROR] PNG file not found: {png_path}")
        return None

    draw = ImageDraw.Draw(img)
    try:
        font_size = 18
        font = ImageFont.load_default(font_size)
    except IOError:
        print("[WARNING] Custom font not found. Using default font.")
        font = None
    average_char_width = (
        6 
    )
    spacing = 5 

    driver_data["Preprocessed Mutation Summary"] = (
        driver_data["Mutation Summary"].str.split().str[0]
    )

    driver_data[["chr", "position"]] = driver_data[
        "Preprocessed Mutation Summary"
    ].str.split(":", expand=True)

    for (parent, child), control_points in edge_info.items():
        gene_names_set = set()
        for mutation_summary, mutations in matched_positions.items():
            chromosome, position = mutation_summary.split(":")

            matching_driver_data = driver_data[
                (driver_data["chr"] == chromosome)
                & (driver_data["position"] == position)
            ]

            for mutation in mutations:
                if len(mutation) < 2:
                    continue

                file_name = mutation[0]
                if file_name == child:
                    gene_names = matching_driver_data["Driver Gene"].values
                    gene_names_set.update(gene_names)

        if gene_names_set and len(control_points) >= 5:
            midpoint_graphviz = (
                (control_points[1][0] + control_points[4][0]) / 2,
                (control_points[1][1] + control_points[4][1]) / 2,
            )
            midpoint_png = graphviz_to_png_coordinates(
                midpoint_graphviz, img.height, dpi=graphviz_dpi
            )

            shift_right = 5

            dot_radius = 6
            dot_center = (midpoint_png[0] + shift_right, midpoint_png[1])
            dot_bbox = [
                dot_center[0] - dot_radius,
                dot_center[1] - dot_radius,
                dot_center[0] + dot_radius,
                dot_center[1] + dot_radius,
            ]
            draw.ellipse(dot_bbox, fill="blue")

            gene_names_str = ", ".join(sorted(gene_names_set))
            text_position = (dot_center[0] + dot_radius + 5, dot_center[1] - dot_radius)
            draw.text(text_position, gene_names_str, fill="black", font=font)

    modified_png_path = os.path.join(
        target_dir, os.path.basename(png_path).replace(".png", "_driver.png")
    )
    img.save(modified_png_path)
    print(f"Driver mutations annotated and saved to {modified_png_path}")


def graphviz_to_png_coordinates(point, img_height, dpi=96):
    scaling_factor = dpi / 72.0
    x_pixels = point[0] * scaling_factor
    y_pixels = img_height - point[1] * scaling_factor
    return x_pixels, y_pixels


def load_and_resize_image(image_path, new_size):
    try:
        image = Image.open(image_path)
        resized_image = image.resize(new_size, Image.Resampling.LANCZOS)
        return resized_image
    except FileNotFoundError:
        print(f"[ERROR] Image file not found: {image_path}")
        return None


def calculate_edge_midpoint(control_points, i1, i2):
    if (
        not control_points
        or i1 < 0
        or i2 >= len(control_points)
        or i1 >= len(control_points)
        or i2 < 0
    ):
        return None

    i1, i2 = sorted([i1, i2])

    x1, y1 = control_points[i1]
    x2, y2 = control_points[i2]

    midpoint_x = (x1 + x2) / 2
    midpoint_y = (y1 + y2) / 2

    return (midpoint_x, midpoint_y)


def find_best_midpoint(control_points, actual_midpoint):
    min_distance = float("inf")
    best_midpoint = None
    best_pair = None
    for i in range(len(control_points)):
        for j in range(i + 1, len(control_points)):
            midpoint = (
                (control_points[i][0] + control_points[j][0]) / 2,
                (control_points[i][1] + control_points[j][1]) / 2,
            )

            distance = euclidean_distance(midpoint, actual_midpoint)

            if distance < min_distance:
                min_distance = distance
                best_midpoint = midpoint
                best_pair = (control_points[i], control_points[j])

    return best_pair, best_midpoint, min_distance


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

