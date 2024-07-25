import argparse
import os
import pandas as pd
from graphviz import Digraph
from PIL import Image, ImageDraw, ImageFont
import glob
import re


def parse_meltos_output_for_clones(meltos_file_path):
    node_clone_mapping = {}
    with open(meltos_file_path, "r") as file:
        for line in file:
            if line.startswith("Valid"):
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                node_id = parts[2]
                clone_ids = []
                array_end_found = False
                for part in parts[4:]:
                    if array_end_found:
                        clone_ids.extend(
                            [int(x) for x in part.split(",") if x.isdigit()]
                        )
                    elif part.startswith("[") and part.endswith("]"):
                        array_end_found = True
                    elif part.endswith("]"):
                        array_end_found = True
                node_clone_mapping[node_id] = clone_ids
    return node_clone_mapping


def save_node_clone_mapping_to_tsv(node_clone_mapping, output_tsv_path):
    print(f"Saving node-clone mapping to {output_tsv_path}...")
    with open(output_tsv_path, "w") as file:
        file.write("Node_ID\tClone_IDs\n")
        for node_id, clone_ids in node_clone_mapping.items():
            file.write(f"{node_id}\t{','.join(map(str, clone_ids))}\n")
    print("Saving completed.")


def find_most_recent_meltos_file(target_dir):
    search_pattern = os.path.join(target_dir, "*meltos*.txt")
    meltos_files = glob.glob(search_pattern)

    if not meltos_files:
        print("No .meltos files found in the specified directory")
        return None

    most_recent_file = max(meltos_files, key=os.path.getctime)
    return most_recent_file


def read_driver_mutation_file(driver_file_path):
    driver_mutations = []
    with open(driver_file_path, "r") as file:
        next(file)
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                print(f"Skipping line due to unexpected format: {line.strip()}")
                continue
            mutation_summary_part = parts[4].split(" ")[0]
            driver_gene = parts[0]

            if "-" in mutation_summary_part:
                base, positions = mutation_summary_part.split(":")
                start_pos, end_pos = positions.split("-")
                mutation_summary = (
                    f"{base}:{start_pos}"
                    if start_pos == end_pos
                    else [f"{base}:{start_pos}", f"{base}:{end_pos}"]
                )
                if isinstance(mutation_summary, list):
                    for summary in mutation_summary:
                        driver_mutations.append((summary, driver_gene))
                    continue
            else:
                mutation_summary = mutation_summary_part

            driver_mutations.append((mutation_summary, driver_gene))
    return driver_mutations


def match_driver_mutations_to_nodes(
    driver_mutations, matches_file_path, output_file_path
):
    with open(matches_file_path, "r") as matches_file, open(
        output_file_path, "w"
    ) as output_file:
        next(matches_file)
        output_file.write("Node_ID\tClone_ID\tPosition\tDriver_Gene\n")

        for line in matches_file:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            node_id, clone_id, position_1, position_2 = parts

            for position in [position_1, position_2]:
                chrom, pos = re.split("[:-]", position)[:2]

                for driver_summary, driver_gene in driver_mutations:
                    driver_chrom, driver_pos = re.split("[:-]", driver_summary)[:2]

                    if chrom == driver_chrom and pos == driver_pos:
                        output_file.write(
                            f"{node_id}\t{clone_id}\t{position}\t{driver_gene}\n"
                        )
                        break


def match_node_id_with_sv_and_save(tsv_file_path, sv_file_path, output_matches_path):
    node_to_clones = {}
    with open(tsv_file_path, "r") as file:
        next(file)
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            node_id = parts[0]
            clone_ids = parts[1].split(",") if parts[1] else []
            node_to_clones[node_id] = clone_ids

    sv_details = {}
    with open(sv_file_path, "r") as file:
        next(file)
        for line in file:
            parts = line.strip().split("\t")
            sv_id = parts[0]
            sv_positions = parts[1].split("-")
            if len(sv_positions) == 2:
                position_1, position_2 = sv_positions
            else:
                position_1, position_2 = parts[1], ""
            sv_details[sv_id] = (position_1, position_2)

    with open(output_matches_path, "w") as output_file:
        output_file.write("Node_ID\tClone_ID\tPosition_1\tPosition_2\n")
        for node_id, clone_ids in node_to_clones.items():
            for clone_id in clone_ids:
                if clone_id in sv_details:
                    position_1, position_2 = sv_details[clone_id]
                    output_file.write(
                        f"{node_id}\t{clone_id}\t{position_1}\t{position_2}\n"
                    )

    print(f"Matches successfully written to {output_matches_path}")


def read_node_clone_mapping(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep="\t")

    node_clone_mapping = pd.Series(
        mapping_df.Clade.values, index=mapping_df.NodeID
    ).to_dict()

    return node_clone_mapping


def plot_meltos_tree(
    meltos_output_path,
    node_clone_mapping_file,
    target_dir,
    driver_file_path,
    output_filename="meltos_output_tree",
    ranksep=1,
    nodesep=1,
):
    dot = Digraph(comment="Meltos Tree")
    graph_title = "Meltos Tree Visualization with Clone Annotation"
    dot.attr(ranksep=str(ranksep), nodesep=str(nodesep))
    dot.attr("graph", label=graph_title, labelloc="t", fontsize="20", fontcolor="black")
    relationships = []
    with open(meltos_output_path, "r") as file:
        for line in file:
            if "->" in line:
                parts = line.strip().split("->")
                parent = (
                    parts[0].strip().split()[0]
                )  
                child = (
                    parts[1].strip().split()[0]
                )  
                relationships.append((parent, child))

    node_clone_mapping = read_node_clone_mapping(node_clone_mapping_file)

    for parent, child in relationships:
        if parent == "0":
            parent_label = "Normal"
        else:
            parent_label = f"{parent} ({node_clone_mapping.get(int(parent), 'Unknown').replace(';', ',')})"

        child_label = f"{child} ({node_clone_mapping.get(int(child), 'Unknown').replace(';', ',')})"

        if parent_label == "Normal":
            dot.node(
                parent, parent_label, shape="ellipse"
            ) 
        else:
            dot.node(parent, parent_label)

        dot.node(child, child_label)
        dot.edge(parent, child)

    dot_file_path = os.path.join(target_dir, "SV_tree")
    print(f"Rendering DOT file: {dot_file_path}")
    dot.render(dot_file_path, view=False, format="dot")

    # Generate the Graph Output
    output_file_path = os.path.join(target_dir, output_filename)
    dot.render(output_file_path, format="png", view=False)
    print(f"Tree diagram saved to: {output_file_path}.png")

    dot_file_path = os.path.join(target_dir, "SV_tree.dot")
    txt_file_path = os.path.join(target_dir, "SV_tree.txt")
    edges_txt_path = os.path.join(target_dir, "SV_tree_Edges.txt")

    save_edge_info_to_txt(dot_file_path, edges_txt_path)

    edge_info = parse_edge_info(edges_txt_path)
    print(f"Edge information: {edge_info}")
    output_file_path = output_file_path
    driver_gene_by_node = read_driver_matches(driver_file_path)
    mark_midpoints_on_png(
        output_file_path + ".png", edge_info, driver_gene_by_node, target_dir
    )
    print("Midpoints marking completed.")
    print("Phylogenetic Tree created.")


def mark_midpoints_on_png(
    png_path, edge_info, driver_gene_by_node, target_dir, graphviz_dpi=96.0
):
    try:
        img = Image.open(png_path)
        draw = ImageDraw.Draw(img)
        font = ImageFont.truetype("arial.ttf", 15) 
    except FileNotFoundError:
        print(f"[ERROR] PNG file not found: {png_path}")
        return

    for (parent, child), control_points in edge_info.items():
        if (
            child in driver_gene_by_node
        ): 
            if len(control_points) >= 4:
                midpoint_graphviz = (
                    (control_points[1][0] + control_points[4][0]) / 2,
                    (control_points[1][1] + control_points[4][1]) / 2,
                )
                midpoint_png = graphviz_to_png_coordinates(
                    midpoint_graphviz, img.height, dpi=graphviz_dpi
                )

                radius = 5  
                top_left = (midpoint_png[0] - radius, midpoint_png[1] - radius)
                bottom_right = (midpoint_png[0] + radius, midpoint_png[1] + radius)
                draw.ellipse([top_left, bottom_right], fill="blue")

                gene_name = driver_gene_by_node[child]
                draw.text(
                    (midpoint_png[0] + radius + 5, midpoint_png[1] - radius),
                    gene_name,
                    fill="black",
                    font=font,
                )

    modified_png_path = png_path.replace(".png", "_driver_annotated.png")
    img.save(modified_png_path)
    print(f"Driver mutations annotated and saved to {modified_png_path}")


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


def graphviz_to_png_coordinates(point, img_height, dpi=96):
    scaling_factor = dpi / 72.0
    x_pixels = point[0] * scaling_factor
    y_pixels = img_height - point[1] * scaling_factor
    return x_pixels, y_pixels


def read_driver_matches(driver_matches_path):
    driver_gene_by_node = {}
    df = pd.read_csv(driver_matches_path, sep="\t")
    print(df)
    for _, row in df.iterrows():
        node_id = str(row["Node_ID"])
        driver_gene = row["Driver_Gene"]
        driver_gene_by_node[node_id] = driver_gene
    return driver_gene_by_node


def plot_meltos(target_dir, sv_file_path, driver_file_path,input_file_base_name):
    meltos_file_path = find_most_recent_meltos_file(target_dir)
    
    if not meltos_file_path:
        print("Failed to find a Meltos output file. Please check the directory.")
        return
    
    print(f"Processing Meltos output file: {meltos_file_path}")
    
    node_clone_mapping = parse_meltos_output_for_clones(meltos_file_path)
    output_tsv_path = os.path.join(target_dir, "node_clone_mapping.tsv")
    save_node_clone_mapping_to_tsv(node_clone_mapping, output_tsv_path)
    
    output_matches_path = os.path.join(target_dir, "matches.tsv")
    match_node_id_with_sv_and_save(output_tsv_path, sv_file_path, output_matches_path)
    
    driver_mutations = read_driver_mutation_file(driver_file_path)
    driver_matches_path = os.path.join(target_dir, "driver_matches.tsv")
    match_driver_mutations_to_nodes(driver_mutations, output_matches_path, driver_matches_path)
    
    print("Integration and matching process complete.")
    node_clone_mapping_file = os.path.join(target_dir, f"{input_file_base_name}snv_CloneFinderCloneID.txt")

    
    plot_meltos_tree(meltos_file_path, node_clone_mapping_file, target_dir, driver_matches_path)
    
    print("Phylogenetic Tree Creation Complete.")
