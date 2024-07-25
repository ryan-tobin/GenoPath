import subprocess
import os


def clone_git_repo(repo_url, target_dir):
    path_dir = os.path.join(target_dir, "PathFinder")
    if not os.path.isdir(path_dir):
        subprocess.run(["git", "clone", repo_url, path_dir], check=True)


def create_output_directory(target_dir):
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


def run_pathfinder(
    aln_phylosig,
    processed_clone_presence_output,
    primary,
    max_graphs_per_tree,
    target_dir,
):
    create_output_directory(target_dir)

    pathfinder_repo_url = "https://github.com/SayakaMiura/PathFinder"
    clone_git_repo(pathfinder_repo_url, target_dir)

    original_dir = os.getcwd()
    pathfinder_dir = os.path.join(target_dir, "PathFinder")
    os.chdir(pathfinder_dir)

    try:
        command = [
            "python",
            "pathfinder.py",
            aln_phylosig,
            processed_clone_presence_output,
            "-o",
            target_dir,
        ]

        if primary is not None:
            command.extend(["--primary", primary])

        if max_graphs_per_tree is not None:
            command.extend(["--max_graphs_per_tree", str(max_graphs_per_tree)])

        subprocess.run(command, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error in running PathFinder: {e}")
    finally:
        os.chdir(original_dir)
