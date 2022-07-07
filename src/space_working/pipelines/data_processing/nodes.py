import pandas as pd
import subprocess
import glob
import pdb
import os
from kedro.extras.datasets.text import TextDataSet
from kedro.extras.datasets.biosequence import BioSequenceDataSet
import re

def _is_true(x: pd.Series) -> pd.Series:
    return x == "t"


def _parse_percentage(x: pd.Series) -> pd.Series:
    x = x.str.replace("%", "")
    x = x.astype(float) / 100
    return x


def _parse_money(x: pd.Series) -> pd.Series:
    x = x.str.replace("$", "").str.replace(",", "")
    x = x.astype(float)
    return x


def walkthrough_prokka_bins(prokka_dir: str, gff_dir: str, fasta_dir: str) -> [TextDataSet, TextDataSet]:
    gff_files = []
    fasta_files = []
    for root, dirs, files in os.walk(prokka_dir): # Each bin in each experimental type
        for name in files:
            if name.endswith('gff'): # Main file containing annotations
                gff_files.append(os.path.join(root,name))
                fasta_files.append(os.path.join(root,name[:name.rfind('.')] +'.faa'))
    return ','.join(gff_files), ','.join(fasta_files)

def parse_gff(gff_files: str) -> pd.DataFrame: #, prokka_bins_file: str) -> pd.DataFrame:
    gff_cont = []
    unique_prokka_id = {}
    for gff in gff_files.split(','):
        with open(gff) as f:
            for line in f:
                if 'ID=' in line:
                    unique_id = re.findall(r'ID=(.+?);',line)[0]
                    if unique_id not in unique_prokka_id:
                        exp_type = gff[gff.rfind('/')+1:gff.rfind('.')]
                        unique_prokka_id[unique_id[:unique_id.find('_')]] = exp_type
                    # I wanted to initially check if these sequences were used in the fasta file
                    # if unique_id in prokka_sequences or unique_id in merged:
                    # Protein annotations
                    if 'UniProtKB:' in line:
                        prod_annot = re.findall(r'UniProtKB:(.+?)[;|,]',line)
                    elif  'HAMAP:' in line:
                        prod_annot = re.findall(r'HAMAP:(.+?)[;|,]',line)
                    else:
                        prod_annot = ''
                        #infer = re.findall(r'inference=(.+?);',line)
                        #if infer:
                        #    infer_all.add(infer[0]) # Check if I missed any annotations
                    if 'gene=' in line:
                        gene_id = re.findall(r'gene=(.+?);',line)[0]
                    else: gene_id = ''

                    scaf = line[:line.find('\t')]
                    product = re.findall(r'product=(.+?)\n',line)
                    gff_cont.append([scaf,unique_id,'/'.join(prod_annot),'/'.join(product),gene_id])
                    prod_annot = ''
                    #else:
                    #    continue
                elif line.startswith('>'):
                    break
                else:
                    continue
    prokka_gff = pd.DataFrame(gff_cont,columns=['scaffold','gid','annot','ainfo','gene'])
    #_write_dataframe(prokka_gff,prokka_gff_file)
    prokka_bins = pd.DataFrame(unique_prokka_id.items(), columns = ['unique_id','bin'])
    #_write_dataframe(prokka_bins,prokka_bins_file)
    return prokka_gff, prokka_bins

def merge_fasta_files(fasta_files: TextDataSet, merged_fasta: str):
    """Merging all the fasta files from every prokka bin into one large file
    """
    if not glob.glob(merged_fasta):
        cmd = f"cat {fasta_files.replace(',',' ')} > {merged_fasta}"
        _ = subprocess.run(cmd,shell=True)
        print(cmd)

def unique_seq_merged_ids(merged_fasta: BioSequenceDataSet) -> [BioSequenceDataSet, pd.DataFrame]:
    """Unique sequences from prokka bins with merged ids
    Saving information for each ID that is merged with another and the resulting merge
    """
    unique_merged = {}
    merged_ids = set()
    for record in merged_fasta:
        sequence = str(record.seq)
        if sequence in unique_merged:
            if unique_merged[sequence].id in merged_ids:
                merged_ids.remove(unique_merged[sequence].id)
            unique_merged[sequence].id = unique_merged[sequence].id + '/' + record.id
            merged_ids.add(unique_merged[sequence].id)
        else:
            unique_merged[sequence] = record

    dict_merged_ids = {}
    for merged_id in merged_ids:
        for original_id in merged_id.split('/'):
            dict_merged_ids[original_id] = merged_id
    df_merged_ids = pd.DataFrame(dict_merged_ids.items(),columns=['single','merged']).set_index('single')
    return list(unique_merged.values()),df_merged_ids

def hypo_prot_sequences(prokka_seq: BioSequenceDataSet, merged_ids: pd.DataFrame, prokka_gff: pd.DataFrame) -> BioSequenceDataSet:
    """Creates sequences of proteins annotated as hypothetical by prokka
    """
    
    hypo_prot = prokka_gff[prokka_gff.ainfo == 'hypothetical protein']
    hypo_prot_list = hypo_prot.gid.tolist() #this list won't work since I merged some'
    hypo_prot_seq= {}
    for x in hypo_prot_list:
        if x in prokka_sequences:
            hypo_prot_seq[x] = prokka_sequences[x]
        else:
            hypo_prot_seq[merged_ids[x].merged] = prokka_sequences[merged_ids[x].merged]
    return hypo_prot_seq

def run_cdhit(input_file,seq_id,output_file):
    seq_id = 0.90
    cdhit_name = f'cdhit_prokka_nometa_allbins_allexp_global_{seq_id}'
    cdhit_file = os.path.join(datasets,cdhit_name)
    if not os.path.isfile(cdhit_file):
        run_cdhit(hypo_prot_seq_file,seq_id,cdhit_file)
    cmd = f'../programs/cdhit/cd-hit -i {input_file} -c {seq_id} -aL {seq_id} -d 0 -o {output_file}'
    _ = subprocess.call(cmd, shell=True)

def preprocess_companies(companies: pd.DataFrame) -> pd.DataFrame:
    """Preprocesses the data for companies.

    Args:
        companies: Raw data.
    Returns:
        Preprocessed data, with `company_rating` converted to a float and
        `iata_approved` converted to boolean.
    """
    return companies.head(40)


def preprocess_shuttles(shuttles: pd.DataFrame) -> pd.DataFrame:
    """Preprocesses the data for shuttles.

    Args:
        shuttles: Raw data.
    Returns:
        Preprocessed data, with `price` converted to a float and `d_check_complete`,
        `moon_clearance_complete` converted to boolean.
    """
    shuttles["d_check_complete"] = _is_true(shuttles["d_check_complete"])
    shuttles["moon_clearance_complete"] = _is_true(shuttles["moon_clearance_complete"])
    shuttles["price"] = _parse_money(shuttles["price"])
    return shuttles.head(40)


def create_model_input_table(
    shuttles: pd.DataFrame, companies: pd.DataFrame, reviews: pd.DataFrame
) -> pd.DataFrame:
    """Combines all data to create a model input table.

    Args:
        shuttles: Preprocessed data for shuttles.
        companies: Preprocessed data for companies.
        reviews: Raw data for reviews.
    Returns:
        Model input table.

    """
    rated_shuttles = shuttles.merge(reviews, left_on="id", right_on="shuttle_id")
    model_input_table = rated_shuttles.merge(
        companies, left_on="company_id", right_on="id"
    )
    model_input_table = model_input_table.dropna()
    return model_input_table
