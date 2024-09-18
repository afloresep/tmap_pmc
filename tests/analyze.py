import pickle

with open("word_count.pickle", "rb") as f:
    word_counts = pickle.load(f)
    print(word_counts.most_common(1000))
