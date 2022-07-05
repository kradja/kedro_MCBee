from Bio import SeqIO
import os
import itertools
import pdb
import pandas as pd
import glob
import re
import subprocess
from typing import Dict, Tuple, List
from kedro.extras.datasets.biosequence import BioSequenceDataSet

def walkthrough_prokka_bins(prokka_dir: str) -> [List[str], List[str]]:
    gff_files = []
    fasta_files = []
    for root, dirs, files in os.walk(prokka_dir): # Each bin in each experimental type
        for name in files:
            if name.endswith('gff'): # Main file containing annotations
                gff_files.append(os.path.join(root,name))
                fasta_files.append(os.path.join(root,name[:name.rfind('.')] +'.faa')) # faa file - bin
    return gff_files, fasta_files

def parse_gff(gff_files: List[str]) -> [List[List[str]], Dict]:
    gff_cont = []
    unique_prokka_id = {}
    for gff in gff_files:
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
    return gff_cont, unique_prokka_id

def _merge_fasta_files(fasta_files: List[str], inter_data: str) -> str:
    merged_file = os.path.join(inter_data,'merged_prokka_bins.fasta')
    if not glob.glob(merged_file):
        cmd = f"cat {' '.join(fasta_files)} > {merged_file}"
        _ = subprocess.run(cmd,shell=True)
        print(cmd)
    return merged_file

def parse_write_unique_fasta(fasta_files: List[str], inter_data: str) -> [BioSequenceDataSet,Dict]:
    merged_file = _merge_fasta_files(fasta_files, inter_data)
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
    ftmp = os.path.join(inter_data,'unique_merged.fasta')
    dataset = BioSequenceDataSet(filepath = ftmp,
                                 load_args={'format':'fasta'},
                                 save_args={'format':'fasta'})
    dataset.save(list(unique_merged.values()))
    dict_merged_ids = {}
    for merged_id in merged_ids:
        for original_id in merged_id.split('/'):
            dict_merged_ids[original_id] = merged_id
    return dataset, dict_merged_ids

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
