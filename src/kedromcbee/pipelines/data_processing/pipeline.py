from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline
from pathlib import Path
import pdb
#from kedro.config import ConfigLoader
import yaml
from kedro.io import DataCatalog
from kedro.framework.hooks.specs import DataCatalogSpecs
from .nodes import walkthrough_prokka_bins, parse_gff 

#, merge_fasta_files, unique_seq_merged_ids

def create_pipeline(**kwargs) -> Pipeline:
    #catalogfile = './conf/base/catalog.yml'
    #with open(catalogfile,'r') as f:
    #    tmp_catalog = yaml.safe_load(f)
    return pipeline([
                node(
                    func= walkthrough_prokka_bins,
                    inputs=["params:prokka_dir",'params:gff_dir','params:fasta_dir'],
                    outputs=['gff_files','fasta_files']
                    ),
                node(
                    func=parse_gff,
                    #inputs=['gff_files','params:prokka_gff_file','params:prokka_bins_file'],
                    inputs=['gff_files','params:prokka_bins_file'],
                    outputs='prokka_gff'
                    )
                #node(
                #    func=merge_fasta_files,
                #    inputs=['fasta_files','params:merged_fasta'],
                #    outputs='rand'
                #    ),
                #node(
                #    func=unique_seq_merged_ids,
                #    inputs='merged_fasta',
                #    outputs=['prokka_seq','merged_ids']
                #    )
                #node(
                #    func=comparethings,
                #    inputs=['umerged','merged_ids','gff_cont','inter_data'],
                #    outputs='hypo_prot')
            ])
