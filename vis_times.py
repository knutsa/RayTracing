import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


import numpy as np

def parse_data(file_path):
    lines = Path(file_path).read_text().split("\n")

    num_runs = int(lines[0].split(":")[1].strip())
    ranges_times = []
    fixed_avg_time = float(lines[-1].split(":")[1].replace("s", "").strip())

    for line in lines[2:-1]:
        parts = line.split(',')
        num_ranges = int(parts[0].split('=')[1])
        avg_time = float(parts[1].split('=')[1].replace('s', '').strip())
        ranges_times.append((num_ranges, avg_time))

    return np.array(ranges_times)


plot_folder = Path(__file__).parent / "benchmark_plots"

Ts = 535.3 #Serial computation time
Ps = [8, 32, 128, 512]
file_paths = [f"data{p}" for p in Ps]

ranges_timess = [parse_data(file_path + ".txt") for file_path in file_paths]


for fn, ranges_times in zip(file_paths, ranges_timess):
    plt.clf()
    plt.plot(ranges_times[:, 0], ranges_times[:, 1])
    plt.scatter(ranges_times[:, 0], ranges_times[:, 1])
    plt.gca().set_xscale("log")
    plt.xlabel("number of ranges")
    plt.ylabel("Execution time (s)")
    plt.savefig(plot_folder / f"{fn}_time_plot.png", bbox_inches="tight")
    plt.clf()
    plt.plot(ranges_times[:, 0], Ts / ranges_times[:, 1])
    plt.scatter(ranges_times[:, 0], Ts / ranges_times[:, 1])
    plt.gca().set_xscale("log")
    plt.xlabel("number of ranges")
    plt.ylabel(r"$S$")
    plt.savefig(plot_folder / f"{fn}_speedup_plot.png")

plt.clf()
for p, ranges_times in zip(Ps, ranges_timess):
    plt.plot(ranges_times[:, 0], Ts / ranges_times[:, 1] / p, label=f"P={p}")
    plt.scatter(ranges_times[:, 0], Ts / ranges_times[:, 1]/ p)
plt.gca().set_xscale("log")
plt.legend()
plt.xlabel("number of ranges")
plt.ylabel(r"$\eta$")
plt.savefig(plot_folder / "efficiencies.png", bbox_inches="tight")