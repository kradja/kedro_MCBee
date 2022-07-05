from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import walkthrough_prokka_bins, parse_gff

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
                node(
                    func= walkthrough_prokka_bins,
                    inputs=["params:prokka_dir",'params:gff_dir','params:fasta_dir'],
                    outputs=['gff_files','fasta_files']
                    ),
                node(
                    func=parse_gff,
                    inputs=['gff_files'],
                    outputs=['gff_df','prokka_bins'])
                #node(
                #    func=parse_write_unique_fasta,
                #    inputs=['fasta_files','inter_data'],
                #    outputs=['umerged','merged_ids'])
                #node(
                #    func=comparethings,
                #    inputs=['umerged','merged_ids','gff_cont','inter_data'],
                #    outputs='hypo_prot')
            ])
