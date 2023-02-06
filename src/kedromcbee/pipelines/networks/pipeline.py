from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import analyze_networks, build_multilayer_network


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=build_multilayer_network,
                inputs=[
                    "prokka_edges",
                    "hierarchy_go",
                    "uni_go",
                    "go_uni",
                    "go_gff_prokka",
                ],
                outputs=["bee_graph","annot_graph"],  # ["bee_mlnetwork","bee_graph"],
                name="build_mlnetwork",
            ),
            node(
                func=analyze_networks,
                inputs=["bee_graph","annot_graph", "go_gff_prokka"],
                outputs="edges_info",
                name="analyze_networks",
            ),
        ],
        namespace="networks",
        inputs=["prokka_edges", "hierarchy_go", "uni_go", "go_uni", "go_gff_prokka"],
        outputs=["bee_graph","annot_graph"],
    )
