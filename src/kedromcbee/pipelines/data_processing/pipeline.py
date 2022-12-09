from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import (
    go_annotations,
    go_ontology,
    hypo_prot_sequences,
    prokka_bins_faa,
    prokka_bins_gff,
    prokka_edges,
)


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
                outputs=["gff_prokka", "annots_gff_prokka", "prokka_bins"],
                name="prokka_bins_gff_node",
            ),
            node(
                func=hypo_prot_sequences,
                inputs=["fasta_prokka", "merged_ids", "gff_prokka"],
                outputs="updated_gff_prokka",
                name="hypo_prot_seq",
            ),
            node(
                func=prokka_edges,
                inputs=["updated_gff_prokka", "prokka_bins"],
                outputs="prokka_edges",
                name="prokka_edges",
            ),
            node(
                func=go_annotations,
                inputs="updated_gff_prokka",
                outputs="go_gff_prokka",
                name="go_annots",
            ),
            node(
                func=go_ontology,
                inputs="go_gff_prokka",
                outputs="hierarchy_gff_prokka",
                name="go_onts",
            ),
        ],
        namespace="data_processing",
        inputs=["partition_prokka_faa", "partition_prokka_gff"],
        outputs=["hierarchy_gff_prokka", "prokka_edges"],
    )
