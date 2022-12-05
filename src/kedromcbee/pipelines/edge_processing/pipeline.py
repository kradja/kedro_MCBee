from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import clustered_hypo_prot_edges, filter_clusters, prokka_annotation_edges# , binned_edges


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=filter_clusters,
                inputs=["cdhit_cluster", "updated_gff_prokka"],
                outputs=["cluster_df","merged_gff_prokka"],
                name="filter_clusters_node",
            ),
            node(
                func=clustered_hypo_prot_edges,
                inputs=["cluster_df", "prokka_bins","merged_gff_prokka"],
                outputs="cdhit_edges",
                name="clustered_hypo",
            ),
            node(
                func=prokka_annotation_edges,
                inputs=["merged_gff_prokka", "prokka_bins"],
                outputs="prokka_edges",
                name="annotated_edges",
            ),
            #node(
            #    func=binned_edges,
            #    inputs=["merged_gff_prokka", "prokka_bins"],
            #    outputs="bin_edges",
            #    name="binned",
            #),
        ],
        namespace="edge_processing",
        inputs=["cdhit_cluster","updated_gff_prokka","prokka_bins"],
        outputs=["cdhit_edges","prokka_edges","merged_gff_prokka"],# "bin_edges",
    )
