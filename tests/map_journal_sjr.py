import pickle
import pandas as pd


def main():
    """ The main function """
    df = pd.read_csv("scimagojr_2018.csv", sep=";")
    data = {}

    for i, row in df.iterrows():
        ids = row[4].replace(" ", "").split(",")
        sjr = row[5]
        categories = [
            s.strip()
            for s in (
                row[18]
                .replace(" (Q1)", "")
                .replace(" (Q2)", "")
                .replace(" (Q3)", "")
                .replace(" (Q4)", "")
                .split(";")
            )
        ]

        for i in ids:
            data[i] = [sjr, categories]

    pickle.dump(
        data, open("issn_sjr_map.pickle", "wb+"),
        protocol=pickle.HIGHEST_PROTOCOL
    )


if __name__ == "__main__":
    main()
