import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import csv

from pathlib import Path

DATA_PATH = Path( "../data")

data = pd.read_csv(DATA_PATH / "sample-2-ATR.csv")

plt.figure(figsize=(10, 6))

plt.plot(
    data["Wavenumber_cm-1"],
    data["Absorbance_au"],
    color="black",
    linewidth=1
)

# FTIR convention: high wavenumber on the left
plt.gca().invert_xaxis()

plt.xlabel("Wavenumber (cm$^{-1}$)")
plt.ylabel("Absorbance (a.u.)")

plt.title("sample-2 ATR 1 cm-1")

plt.tight_layout()

plt.savefig(
    "sample-2-ATR.png",
    dpi=300,
    bbox_inches="tight"
)
