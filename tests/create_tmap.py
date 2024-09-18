import sys
import pickle
import numpy as np
import scipy.stats as ss

from timeit import default_timer as timer

import tmap as tm
import pandas as pd
from faerun import Faerun

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Coniguration for the tmap layout
CFG_TMAP = tm.LayoutConfiguration()
CFG_TMAP.k = 50
CFG_TMAP.kc = 50
CFG_TMAP.sl_scaling_min = 1.0
CFG_TMAP.sl_scaling_max = 1.0
CFG_TMAP.sl_repeats = 1
CFG_TMAP.sl_extra_scaling_steps = 2
CFG_TMAP.placer = tm.Placer.Barycenter
CFG_TMAP.merger = tm.Merger.LocalBiconnected
CFG_TMAP.merger_factor = 2.0
CFG_TMAP.merger_adjustment = 0
CFG_TMAP.fme_iterations = 1000
CFG_TMAP.sl_scaling_type = tm.ScalingType.RelativeToDesiredLength
CFG_TMAP.node_size = 1 / 55
CFG_TMAP.mmm_repeats = 1

DATA = []
JOURNAL = []
LABEL = []
ISSNS = []
DIMS = 2048
ENC = tm.Minhash(DIMS)
LF = tm.LSHForest(DIMS, 128, store=True)

# with open("fingerprints.csv", "r") as f:
#     for line in f:
#         DATA.append(list(map(int, line.split(","))))

DF = pd.read_csv("data.csv", header=None)
for i, row in DF.iterrows():
    LABEL.append(
        row[4].replace("'", "")
        + "__"
        + row[4].replace("'", "")
        + "__"
        + f'<a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/{row[0]}">Read</a>'
    )
    JOURNAL.append(row[5])
    ISSNS.append(str(row[6]).replace("-", "").split(","))

# JOURNAL_MAP = dict([(y, x + 1) for x, y in enumerate(sorted(set(JOURNAL)))])
# JOURNAL = [JOURNAL_MAP[x] for x in JOURNAL]

SJR_MAP = pickle.load(open("issn_sjr_map.pickle", "rb"))

SJR = []
for i, issns in enumerate(ISSNS):
    found = False
    for issn in issns:
        if issn in SJR_MAP:
            SJR.append(float(str(SJR_MAP[issn][0]).replace(",", ".")))
            found = True
            break
    if not found:
        SJR.append(0.0)

JOURNAL_VALUES = []
SIZES = []
for i, value in enumerate(JOURNAL):
    if value == "Nature":
        JOURNAL_VALUES.append(1)
        SIZES.append(1)
    elif value == "Cell":
        JOURNAL_VALUES.append(2)
        SIZES.append(1)
    elif "Angewandte Chemie" in value:
        JOURNAL_VALUES.append(3)
        SIZES.append(1)
    elif value == "Science (New York, N.Y.)":
        JOURNAL_VALUES.append(4)
        SIZES.append(1)
    else:
        JOURNAL_VALUES.append(0)
        SIZES.append(0)

pickle.dump(
    (JOURNAL, LABEL, SJR, ISSNS, JOURNAL_VALUES, SIZES),
    open("properties.pickle", "wb+"),
    protocol=pickle.HIGHEST_PROTOCOL,
)


# JOURNAL, LABEL, SJR, ISSNS, JOURNAL_VALUES, SIZES = pickle.load(
#     open("properties.pickle", "rb")
# )

TAB10 = plt.get_cmap("Set1").colors
COLORS = [TAB10[8], TAB10[0], TAB10[1], TAB10[2], TAB10[3]]
MY_CM = LinearSegmentedColormap.from_list("my_map", COLORS, N=len(COLORS))

SJR = ss.rankdata(np.array(SJR) / max(SJR)) / len(SJR)


def main():
    print("Loaded data ...")
    start = timer()
    # fps = []

    # for row in DATA:
    #     fps.append(tm.VectorUint(list(row)))

    # print("Converted fingerprints ...")

    # LF.batch_add(ENC.batch_from_sparse_binary_array(fps))

    print("Added fingerprints ...")

    LF.index()

    LF.restore("pmc.dat")

    x_tmap, y_tmap, s, t, _ = tm.layout_from_lsh_forest(LF, CFG_TMAP)

    LF.clear()

    print("Created map ...")

    legend_labels = [
        (1, "Nature"),
        (2, "Cell"),
        (3, "Angewandte Chemie"),
        (4, "Science"),
        (0, "Other"),
    ]

    faerun = Faerun(view="front", coords=False, legend_title="")
    faerun.add_scatter(
        "PMC",
        {"x": x_tmap, "y": y_tmap, "c": JOURNAL_VALUES, "labels": LABEL, "s": SIZES},
        colormap=MY_CM,
        # {"x": x_tmap, "y": y_tmap, "c": SJR, "labels": LABEL},
        # colormap="rainbow",
        point_scale=2.5,
        max_point_size=10,
        shader="smoothCircle",
        has_legend=True,
        categorical=True,
        legend_labels=legend_labels,
        legend_title="Journal",
        interactive=True,
    )
    faerun.add_tree(
        "PMC_tree", {"from": s, "to": t}, point_helper="PMC", color="#222222"
    )
    faerun.plot("pmc")

    print(timer() - start)


if __name__ == "__main__":
    main()
