import matplotlib.pyplot as plt
import pandas as pd
import csv

from pathlib import Path

DATA_PATH = Path( "../data")

# load data
data = pd.read_csv(DATA_PATH / "beetroot-ATR-1.csv")

# plot
plt.plot(data["wavenumber cm-1"], data["Absorbance (a.u.)"])
