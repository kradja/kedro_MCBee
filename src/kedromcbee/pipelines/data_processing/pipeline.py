from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import walkthrough_prokka_bins

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
                node(
                    func= walkthrough_prokka_bins,
                    inputs=["params:prokka_dir",'params:gff_dir','params:fasta_dir'],
                    outputs=['gff_files','fasta_files'],
                    )
            ])
