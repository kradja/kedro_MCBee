from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import build_multilayer_network


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=build_multilayer_network,
                inputs=["cdhit_edges","prokka_edges"],
                outputs="bee_mlnetwork",
                name="build_mlnetwork",
            ),
        ],
        namespace="networks",
        inputs=["cdhit_edges", "prokka_edges"],
        outputs="bee_mlnetwork",
    )
