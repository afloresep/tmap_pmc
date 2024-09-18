# fingerprint_generation.py

import pandas as pd
import nltk
from nltk.tokenize import word_tokenize
from rdkit.Chem import rdMHFPFingerprint
from tqdm import tqdm
import pickle
import os
import numpy as np

nltk.download('punkt', quiet=True)

def load_data(file_path):
    return pd.read_csv(file_path)

def preprocess_text(df):
    relevant_columns = ['target_name', 'target_protein_class', 'Target_Taxonomy', 'Target_organism', 'Target_type']
    return df[relevant_columns].apply(lambda x: ' '.join(x.astype(str)), axis=1)

def tokenize_text(text):
    return word_tokenize(text.lower())

def generate_fingerprint(tokenized_text, n_permutations=2048):
    encoder = rdMHFPFingerprint.MHFPEncoder(n_permutations)
    rdkit_fp = encoder.FromStringArray(tokenized_text)
    return np.array(rdkit_fp, dtype=np.uint32)  # Convert to numpy array

def generate_and_save_fingerprints(file_path, output_path):
    
    tqdm.pandas()

    print("Loading data...")
    df = load_data(file_path)
    df = df.drop_duplicates(subset='Target_ID', keep='first')
    
    print("Preprocessing and tokenizing text...")
    preprocessed_text = preprocess_text(df)
    tokenized_texts = preprocessed_text.progress_apply(tokenize_text)
    
    print("Generating fingerprints...")
    fingerprints = tokenized_texts.progress_apply(generate_fingerprint)
    
    numpy_fingerprints = []

    # Transform from RdKit fingerprint object to np.array
    for i in fingerprints:
        numpy_fingerprints.append(np.array(i))

    # # Transform list to array to remove duplicates
    # numpy_fingerprints = np.array(numpy_fingerprints)

    # # Remove duplicated targets 
    # numpy_fingerprints = np.unique(numpy_fingerprints, axis=0)

    # Transform back to list (it makes my life easier in later steps)
    # numpy_fingerprints = numpy_fingerprints.tolist()

    print(f"Saving {len(numpy_fingerprints)} fingerprints to {output_path}")
    with open(output_path, 'wb') as f:
        pickle.dump(numpy_fingerprints, f)
    
    return numpy_fingerprints

if __name__ == "__main__":
    file_path = r"C:\Users\aflor\OneDrive\Work\tmap_pmc\data\dataset.csv"
    output_path = r"C:\Users\aflor\OneDrive\Work\tmap_pmc\data\fingerprints.pkl"
    
    fingerprints = generate_and_save_fingerprints(file_path, output_path)
    print(f"Generated and saved {len(fingerprints)} fingerprints.")