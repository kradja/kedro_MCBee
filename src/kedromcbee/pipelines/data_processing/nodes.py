import os
from typing import Dict, Tuple, List
from kedro.extras.datasets.text import TextDataSet

def walkthrough_prokka_bins(prokka_dir: str, gff_dir: str, fasta_dir: str) -> [TextDataSet, TextDataSet]:
    gff_files = []
    fasta_files = []
    for root, dirs, files in os.walk(prokka_dir): # Each bin in each experimental type
        for name in files:
            if name.endswith('gff'): # Main file containing annotations
                gff_files.append(os.path.join(root,name))
                fasta_files.append(os.path.join(root,name[:name.rfind('.')] +'.faa')) 
    gffdataset = TextDataSet(filepath = gff_dir)
    fastadataset = TextDataSet(filepath = fasta_dir)

    gffdataset.save(','.join(gff_files))
    fastadataset.save(','.join(fasta_files))
    return gffdataset, fastadataset
