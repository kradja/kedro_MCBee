from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import (
    go_annotations,
    go_ontology,
    preprocess_prokka_annotations,
    preprocess_prokka_sequences,
    prokka_bins_gff,
    prokka_edges,
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=preprocess_prokka_sequences,
                inputs="partition_prokka_faa",
                outputs=["fasta_prokka", "merged_ids"],
                name="preprocess_prokka_seq",
            ),
            node(
                func=prokka_bins_gff,
                inputs="partition_prokka_gff",
                outputs=["gff_prokka", "annots_gff_prokka", "prokka_bins"],
                name="prokka_bins_gff_node",
            ),
            node(
                func=preprocess_prokka_annotations,
                inputs=["fasta_prokka", "merged_ids", "gff_prokka"],
                outputs="updated_gff_prokka",
                name="preprocess_prokka_annots",
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
