from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import build_multilayer_network, analyze_networks


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=build_multilayer_network,
                inputs=["cdhit_edges","prokka_edges","bin_edges","merged_gff_prokka"],
                outputs="bee_graph",#["bee_mlnetwork","bee_graph"],
                name="build_mlnetwork",
            ),
            node(
                func=analyze_networks,
                inputs=["bee_graph","merged_gff_prokka"],
                outputs="edges_info",
                name="analyze_networks",
            ),
        ],
        namespace="networks",
        inputs=["cdhit_edges", "prokka_edges","bin_edges","merged_gff_prokka"],
        outputs="bee_graph",
    )
