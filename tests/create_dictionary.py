import pickle
from collections import Counter

import pandas as pd

total = Counter({})

df = pd.read_csv('data.csv', header=None)


with open('word_count_part.pickle', 'wb+') as f_out_part:
    for i, row in df.iterrows():
        if i % 1000 == 0:
            print(i, len(total))
            pickle.dump(total, f_out_part, protocol=pickle.HIGHEST_PROTOCOL)

        word_count = Counter(str(row[7]).split(' '))
        total += word_count

with open('word_count.pickle', 'wb+') as f:
    pickle.dump(total, f, protocol=pickle.HIGHEST_PROTOCOL)
