from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import hypo_prot_sequences, prokka_bins_faa, prokka_bins_gff


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=prokka_bins_faa,
                inputs="partition_prokka_faa",
                outputs=["fasta_prokka", "merged_ids"],
                name="prokka_bins_faa_node",
            ),
            node(
                func=prokka_bins_gff,
                inputs="partition_prokka_gff",
                outputs=["gff_prokka", "prokka_bins"],
                name="prokka_bins_gff_node",
            ),
            node(
                func=hypo_prot_sequences,
                inputs=["fasta_prokka", "merged_ids", "gff_prokka"],
                outputs="hypo_prot",
                name="hypo_prot_seq",
            ),
        ],
        namespace="data_processing",
        inputs=["partition_prokka_faa", "partition_prokka_gff"],
        outputs=["hypo_prot", "prokka_bins"],
    )
