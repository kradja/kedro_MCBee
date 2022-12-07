"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from kedromcbee.pipelines import data_processing as dp
from kedromcbee.pipelines import networks as nk


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    data_processing_pipeline = dp.create_pipeline()
    networks_processing_pipeline = nk.create_pipeline()

    return {
        "__default__": data_processing_pipeline + networks_processing_pipeline,
        "dp": data_processing_pipeline,
        "nk": networks_processing_pipeline,
    }
