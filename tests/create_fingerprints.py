import pickle
import pandas as pd

word_dict = dict(pickle.load(open("word_count.pickle", "rb")))
print(len(word_dict))

all_words = {}
for i, (key, value) in enumerate(word_dict.items()):
    all_words[key] = i

data = []
labels = []


df = pd.read_csv("data.csv", header=None)

with open("fingerprints.csv", "w+") as f_out:
    for i, row in df.iterrows():
        if i % 100 == 0:
            print(i)

        fp = []
        words = set(str(row[7]).split(" "))

        for word in words:
            fp.append(all_words[word])

        # for i, word in enumerate(all_words):
        #     if word in words:
        #         fp.append(i)
        #     if len(fp) == len(words):
        #         break

        f_out.write(",".join(map(str, sorted(fp))) + "\n")


# with open('books_labels.csv', 'w+') as f_out:
#     for row in labels:
#         f_out.write(','.join(row) + '\n')
