from collections import Counter
import pandas as pd


def main():
    """ The main function """

    # Extract the journal names

    # df = pd.read_csv("data.csv")
    # with open("data_journals.csv", "w+") as f:
    #     for i, row in df.iterrows():
    #         f.write(row[5] + "\n")

    # Count Journals

    journals = []
    with open("data_journals.csv", "r") as f:
        for line in f:
            line = line.strip()
            journals.append(line)

    print(Counter(journals).most_common(50))


if __name__ == "__main__":
    main()
