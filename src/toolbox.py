import os
import pickle
import tmap as tm
import numpy as np
from tqdm import tqdm
from fingerprint_generation import generate_and_save_fingerprints
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


def list_to_vectorUint(lst):
    """
    Convert list or numpy array to tm.VectorUint type.
    
    Args:
    lst (list or np.array): The list or array to convert.
    
    Returns:
    tm.VectorUint: Converted VectorUint object.
    """
    return tm.VectorUint(lst)


def safe_create_categories(series):
    """
    Create categories from a pandas Series, handling NaN values.
    
    Args:
    series (pandas.Series): The input data series.
    
    Returns:
    tuple: (labels, data) for Faerun plotting.
    """
    return Faerun.create_categories(series.fillna('Unknown').astype(str))

def map_protein_class(value):
    """Map protein class to a simplified category."""
    if pd.isna(value):
        return np.nan
    
    value = value.lower().strip()
    
    if 'enzyme' in value:
            return 'Enzyme'
    elif 'membrane receptor' in value: 
        return 'Membrane receptor'
    elif ' ion channel' in value:
        return 'Ion Channel'
    elif 'transcription factor' in value: 
        return 'Transcription factor'
    elif 'epigenetic regulator' in value:
        return 'Epigenetic regulator'
    elif 'cytosolic protein' in value:
        return 'cytosolic protein'
    else:
        return 'Other'

def map_target_taxonomy(value):
    if 'Eukaryotes' in value:
        if 'Oxidoreductase' in value:
            return 'Oxidoreductase'
        elif 'Transferase' in value:
            return 'Transferase' 
        elif 'Hydrolase' in value:
            return 'Hydrolase'
        else:
            return 'Eukaryotes'
    elif 'Bacteria' in value: 
        return 'Bacteria'
    elif 'Fungi' in value:
        return 'Fungi'
    elif 'Viruses' in value: 
        return 'Viruses'
    elif 'unclassified' in value:
        return 'Unclassified'
    else:
        return 'Other'

def map_target_organism(value):
    """Map target organism to a simplified category."""
    if 'sapiens' in value:
        return 'Homo sapiens'
    elif 'virus' in value:
        return 'Virus'
    elif any(organism in value for organism in ['rattus', 'Musculus']):
        return 'Rat'
    elif 'taurus' in value:
        return 'Bovid'
    elif any(organism in value for organism in ['scrofa', 'Macaca', 'porcellus', 'oryctolagus', 'canis', 'Cricetulus']):
        return 'Other mammals'
    elif any(bacteria in value for bacteria in ['Mycobacterium', 'Escherichia', 'Salmonella', 'Staphylococcus', 'Pseudomonas', 'Bacillus', 'Acinetobacter']):
        return 'Bacteria'
    elif any(parasite in value for parasite in ['Plasmodium', 'Trypanosoma', 'Schistosoma', 'Leishmania']):
        return 'Parasites'
    else:
        return 'Others'
