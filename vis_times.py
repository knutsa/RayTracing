import sys
import re
import numpy as np
import matplotlib.pyplot as plt

ranges_times = np.array([
    (8, 158.574),
    (80, 99.0333),
    (800, 89.9568),
    (8000, 89.542),
    (80000, 91.7084),
    (800000, 125.055),
])

plt.plot(ranges_times[:, 0], ranges_times[:, 1])
plt.scatter(ranges_times[:, 0], ranges_times[:, 1])
plt.gca().set_xscale("log")
plt.xlabel("number of ranges")
plt.ylabel("Execution time (s)")
plt.savefig("benchmarks.png", bbox_inches="tight")

plt.figure()
plt.plot(ranges_times[:, 0], 535.3 / ranges_times[:, 1])
plt.scatter(ranges_times[:, 0], 535.3 / ranges_times[:, 1])
plt.gca().set_xscale("log")
plt.xlabel("number of ranges")
plt.ylabel("Sppedup factor")
plt.savefig("speedup.png", bbox_inches="tight")