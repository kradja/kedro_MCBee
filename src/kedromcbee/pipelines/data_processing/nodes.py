import os
from Bio import SeqIO
import re
import glob
import subprocess
import pdb
from typing import Dict, Tuple, List
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.extras.datasets.text import TextDataSet
import pandas as pd
#from kedro.extras.datasets.pandas import CSVDataSet

def _write_dataframe(input_df: pd.DataFrame, file_location: str):
    input_df.to_csv(file_location,sep = ',')

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

def parse_gff(gff_files: TextDataSet, prokka_bins_file: str) -> pd.DataFrame: 
    gff_cont = []
    unique_prokka_id = {}
    for gff in gff_files.load().split(','):
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
    prokka_bins = pd.DataFrame(unique_prokka_id.items())
    _write_dataframe(prokka_bins,prokka_bins_file)
    return prokka_gff

























def merge_fasta_files(fasta_files: TextDataSet, merged_fasta: str) -> str:
    if not glob.glob(merged_fasta):
        tmp = fasta_files.load()
        cmd = f"cat {tmp.replace(',',' ')} > {merged_fasta}"
        _ = subprocess.run(cmd,shell=True)
        print(cmd)
    return cmd

def unique_seq_merged_ids(merged_fasta: BioSequenceDataSet) -> [BioSequenceDataSet, pd.DataFrame]:
    unique_merged = {}
    merged_ids = set()
    with open(merged_file) as f:
        for record in SeqIO.parse(f,'fasta'):
            sequence = str(record.seq)
            if sequence in unique_merged:
                if unique_merged[sequence].id in merged_ids:
                    merged_ids.remove(unique_merged[sequence].id)
                unique_merged[sequence].id = unique_merged[sequence].id + '/' + record.id
                merged_ids.add(unique_merged[sequence].id)
            else:
                unique_merged[sequence] = record
    #ftmp = os.path.join(inter_data,'prokka_sequences.fasta')
    #prokka_seq = BioSequenceDataSet(filepath = ftmp,
    #                                load_args={'format':'fasta'},
    #                                save_args={'format':'fasta'})
    #prokka_seq.save(list(unique_merged.values()))
    dict_merged_ids = {}
    for merged_id in merged_ids:
        for original_id in merged_id.split('/'):
            dict_merged_ids[original_id] = merged_id
    tmp = pd.DataFrame(dict_merged_ids.items())
    #_write_dataframe(tmp,'')
    return list(unique_merged.values()), tmp

def comparethings(umerged: BioSequenceDataSet, merged_ids: Dict, gff_cont: List[List[str]],inter_data: str) -> str:
    sequence_list = umerged.load()
    prokka_sequences = {}
    for sequence in sequence_list:
        prokka_sequences[sequence.id] = sequence
    gff_df = pd.DataFrame(gff_cont)
    gff_df.columns = ['scaffold','gid','annot','ainfo','gene']
    hypo_prot = gff_df[gff_df.ainfo == 'hypothetical protein']
    hypo_prot_list = hypo_prot.gid.tolist() #this list won't work since I merged some'
    # things in the hypo_prot_list got merged so I just need to get the latest merged from previous function
    hypo_prot_sequences = {}
    for x in hypo_prot_list:
        if x in prokka_sequences:
            hypo_prot_sequences[x] = prokka_sequences[x]
        else:
            hypo_prot_sequences[merged_ids[x]] = prokka_sequences[merged_ids[x]]
    file_name = os.path.join(inter_data,'prokka_hypo_prot_allbins_allexp.fa')

    #write_fasta(hypo_prot_sequences,temp_dir,file_name)
    SeqIO.write(prokka_sequences.values(),file_name,'fasta')
    return file_name
