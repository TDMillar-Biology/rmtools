
from .rm_track import load_data
from pathlib import Path
import matplotlib.pyplot as plt

def run_from_cli(args):
    df = load_data(Path(args.rm))


    sizes = df["end"] - df["start"]

    plt.hist(sizes, bins=500)
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("Distribution of Sizes")
    plt.show()