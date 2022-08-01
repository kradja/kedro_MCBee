"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from kedromcbee.pipelines import data_processing as dp
from kedromcbee.pipelines import edge_processing as ep
from kedromcbee.pipelines import networks as nk
#from space_working.pipelines import data_science as ds


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    data_processing_pipeline = dp.create_pipeline()
    edge_processing_pipeline = ep.create_pipeline()
    networks_processing_pipeline = nk.create_pipeline()
    
    #data_science_pipeline = ds.create_pipeline()

    return {
        "__default__": data_processing_pipeline + edge_processing_pipeline + networks_processing_pipeline,
        "dp": data_processing_pipeline,
        "ep": edge_processing_pipeline,
        "nk": networks_processing_pipeline,
    }
