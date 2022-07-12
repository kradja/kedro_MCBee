from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import filter_clusters, clustered_hypo_prot_edges


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=filter_clusters,
                inputs="cdhit_cluster",
                outputs='cluster_df',
                name="filter_clusters_node",
            ),
            node(
                func=clustered_hypo_prot_edges,
                inputs=['cluster_df','prokka_bins'],
                outputs='cdhit_edges',
                name="clustered_hypo"
            ),
        ],
        namespace="edge_processing",
        inputs="cdhit_cluster",
        outputs="cluster_df",
    )
