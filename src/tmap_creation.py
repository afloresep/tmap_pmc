# tmap_creation.py

import os
import pickle
import tmap as tm
import numpy as np
from tqdm import tqdm
from fingerprint_generation import generate_and_save_fingerprints
from toolbox import list_to_vectorUint, safe_create_categories, map_protein_class, map_target_organism, map_target_taxonomy
import pandas as pd
import logging
import warnings
from tqdm import tqdm
import tmap as tm
import pandas as pd
import numpy as np
from timeit import default_timer as timer
from faerun import Faerun
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from mapchiral.mapchiral import encode as mapc_enc
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import pickle
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import html

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def load_or_generate_fingerprints(fingerprint_path, csv_path):
    if os.path.exists(fingerprint_path):
        print(f"Loading existing fingerprints from {fingerprint_path}")
        with open(fingerprint_path, 'rb') as f:
            return pickle.load(f)
    else:
        print(f"No existing fingerprints found. Generating new fingerprints...")
        return generate_and_save_fingerprints(csv_path, fingerprint_path)


def plot_faerun(x, y, s, t, df):

    """
    Plot the data using Faerun.
    
    Args:
    x (list): X coordinates.
    y (list): Y coordinates.
    s (list): Source nodes for tree plot.
    t (list): Target nodes for tree plot.
    df (pandas.DataFrame): DataFrame with target data.
    """
    f = Faerun(view="front", coords=False, clear_color="#FFFFFF")

    # Create categories
    protein_class_labels, protein_class_data = safe_create_categories(df['mapped_protein_class'])
    taxonomy_labels, taxonomy_data = safe_create_categories(df['map_target_taxonomy'])
    organism_labels, organism_data = safe_create_categories(df['map_target_organism'])

    labels = []

    for i, row in df.iterrows():

        target_name = str(row["target_name"]).strip()
        target_name = html.escape(target_name)
        if not target_name:
            target_name = "N/A"

        target_type_name = str(row['Target_type']).strip()
        target_type_name = html.escape(target_type_name)
        if not target_type_name:
            target_type_name = "N/A"

        target_organism_name = str(row["Target_organism"]).strip()  # Convert to string and remove leading/trailing whitespace
        target_organism_name = html.escape(target_organism_name)  # Escape special characters
        if not target_organism_name:
            target_organism_name = "N/A"  # Provide a default value if empty

        target_taxonomy_name = str(row["Target_Taxonomy"]).strip()  # Convert to string and remove leading/trailing whitespace
        target_taxonomy_name = html.escape(target_taxonomy_name)  # Escape special characters
        if not target_taxonomy_name:
            target_taxonomy_name = "N/A"  # Provide a default value if empty

        mapped_protein_class = str(row["target_protein_class"]).strip()  # Convert to string and remove leading/trailing whitespace
        mapped_protein_class = html.escape(mapped_protein_class)  # Escape special characters
        if not mapped_protein_class:
            mapped_protein_class = "N/A"  # Provide a default value if empty
        
        labels.append(
            row['canonical_smiles']
            + '__'
            + f'<a target="_blank" href="https://www.ebi.ac.uk/chembl/target_report_card/{row["Target_ID"]}">{target_name}</a><br>'
            + '__'
            + f'<small style="font-size:14px;color:grey;">Target Type: {target_type_name}</small>'
            + '__'
            + f'<small style="font-size:14px;color:grey;">Target organism: {target_organism_name}</small>'
            + '__'
            + f'<small style="font-size:14px;color:grey;">Target taxonomy: {target_taxonomy_name}</small>'
            + '__'
            + f'<small style="font-size:14px;color:grey;">Target protein class: {mapped_protein_class}</small>'
        )

    # Add scatter plot
    f.add_scatter(
        "mapc_targets",
        {
            "x": x,
            "y": y,
            "c": [protein_class_data, taxonomy_data, organism_data],
            "labels": labels,
        },
        shader="smoothCircle",
        point_scale=4.0,
        max_point_size=20,
        interactive=True,
        legend_labels=[protein_class_labels, taxonomy_labels, organism_labels],
        categorical=[True, True, True],
        colormap=['tab10', 'tab10', 'tab10'],
        series_title=['Protein Class', 'Target Taxonomy', 'Target Organism'],
        has_legend=True,
    )

    # Add tree
    f.add_tree("mapc_targets_tree", {"from": s, "to": t}, point_helper="mapc_targets", color="#222222")
    
    # Plot
    f.plot('mapc_targets', template='smiles')

def main():
    csv_path = r"C:\Users\biolab\Desktop\Alex\Alex's\OneDrive\Work\tmap_pmc\data\dataset.csv"
    fingerprint_path = r"C:\Users\biolab\Desktop\Alex\Alex's\OneDrive\Work\tmap_pmc\data\fingerprint.pkl"

    df = pd.read_csv(csv_path)

    fingerprints = pd.Series(load_or_generate_fingerprints(fingerprint_path, csv_path))
    fingerprints = fingerprints.apply(list_to_vectorUint)
    
    
    # Apply the mapping function to the column to reduce number of unique values and make it color codeable 
    df['mapped_protein_class'] = df['target_protein_class'].apply(map_protein_class)
    df['map_target_taxonomy'] = df['Target_Taxonomy'].apply(map_target_taxonomy)
    df['map_target_organism'] = df['Target_organism'].apply(map_target_organism)

    logger.info('Indexing...')
    lf = tm.LSHForest(512, 128, store=True)
    lf.batch_add(fingerprints)
    lf.index()

    logger.info('Setting layout')
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 40
    cfg.mmm_repeats = 2
    cfg.sl_extra_scaling_steps = 10
    cfg.k = 20
    cfg.sl_scaling_type = tm.RelativeToAvgLength
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
    
    logger.info('Plotting...')

    plot_faerun(x, y, s, t, df)
    
if __name__ == "__main__":
    start_time = timer()
    main()
    end_time = timer()
    logger.info('TMAP successfully generated.')
    logger.info(f"Total execution time: {(end_time - start_time)/60:.2f} minutes")