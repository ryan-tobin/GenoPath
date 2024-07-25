import subprocess
import os


def run_make_phylosignare_input(aln_phylosig, target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    original_dir = os.getcwd()
    phylosignare_dir = os.path.join(target_dir, "CloneFinder")
    os.chdir(phylosignare_dir)

    try:
        subprocess.run(
            ["python", "make_PhyloSigFinder_Input.py", aln_phylosig], check=True
        )
    except subprocess.CalledProcessError as e:
        print(f"Error creating PhyloSignare Input file: {e}")
    finally:
        os.chdir(original_dir)


def makePhyloInput(aln_phylosig, target_dir):
    run_make_phylosignare_input(aln_phylosig, target_dir)
