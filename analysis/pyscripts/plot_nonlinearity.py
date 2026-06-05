#!/usr/bin/env python3
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot non_linear histogram for W (PbWO4) modules")
parser.add_argument("json", nargs="?", default="new_calibration_NonLinearity.json",
                    help="Path to calibration JSON file")
parser.add_argument("-o", "--output", default=None, help="Save plot to file (e.g. nl.png)")
args = parser.parse_args()

with open(args.json) as f:
    data = json.load(f)

nl_values = [entry["non_linear"] for entry in data
             if entry.get("name", "").startswith("W") and entry["non_linear"] != 0.0]

print(f"W modules with non-zero nl: {len(nl_values)}")
print(f"  mean  = {np.mean(nl_values):.4f}")
print(f"  std   = {np.std(nl_values):.4f}")
print(f"  min   = {np.min(nl_values):.4f}")
print(f"  max   = {np.max(nl_values):.4f}")

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(nl_values, bins=50, color="steelblue", edgecolor="black", linewidth=0.5)
ax.set_xlabel("nl  (non-linearity parameter)", fontsize=13)
ax.set_ylabel("Counts", fontsize=13)
ax.set_title("PbWO$_4$ Module Non-linearity Distribution", fontsize=14)
ax.axvline(np.mean(nl_values), color="red", linestyle="--", linewidth=1.5,
           label=f"mean = {np.mean(nl_values):.3f}")
ax.legend(fontsize=11)
plt.tight_layout()

if args.output:
    plt.savefig(args.output, dpi=150)
    print(f"Saved to {args.output}")
else:
    plt.show()
