import pandas as pd
from ete3 import Tree, TreeStyle, TextFace
import sys

def clone_tree(presence, tree_file,clone_path):
    cpdf = pd.read_csv(presence, sep='\t', index_col='Tumor')
    tree = Tree(tree_file)

    tcn = cpdf.columns[0]
    fsn = None
    for leaf in tree.iter_leaves():
        if leaf.name in cpdf.columns:
            prd = cpdf[leaf.name].dropna()

            fv = prd.iloc[0]
            sn = prd.name

            print("First value:", fv)
            print("Series name:", sn)

            if fsn is None:
                fsn = sn

            print("First series name encountered:", fsn)

            tcn = fsn
            i = 0
            for t, p in prd.items():
                color = "blue" if p > 0 else "black"

                l = f"{p:.3f}"

                if leaf.name == tcn:
                    l = f"{t}\n{l}"

                pf = TextFace(l, fgcolor=color)
                pf.margin_bottom = 5
                pf.margin_top = 5
                pf.margin_right = 10

                leaf.add_face(pf, column=i,
                              position="aligned")
                i += 1

    ts = TreeStyle()
    ts.show_scale = False
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    
    tree.render(clone_path, tree_style=ts)


def tumor_tree(presence, tree_file,tree_path):
    cpdf = pd.read_csv(presence, sep='\t', index_col='Tumor')
    cpdf = cpdf.T
    tree = Tree(tree_file)
    fsn = None
    for leaf in tree.iter_leaves():
        if leaf.name in cpdf.columns:
            prd = cpdf[leaf.name].dropna()

            fv = prd.iloc[0]
            sn = prd.name

            print("First value:", fv)
            print("Series name:", sn)

            if fsn is None:
                fsn = sn

            print("First series name encountered:", fsn)

            tcn = fsn
            i = 0
            for t, p in prd.items():
                color = "blue" if p > 0 else "black"

                l = f"{p:.3f}"

                if leaf.name == tcn:
                    l = f" {t}\n{l}"

                pf = TextFace(l, fgcolor=color)
                pf.margin_bottom = 5
                pf.margin_top = 5
                pf.margin_right = 10

                leaf.add_face(pf, column=i,
                              position="aligned")
                i += 1

    ts = TreeStyle()
    ts.show_scale = False
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    tree.render(tree_path, tree_style=ts)

