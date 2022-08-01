import re
from pathlib import PurePosixPath
from typing import Any, Dict, List

import fsspec
import pandas as pd
from itertools import chain
from kedro.io import AbstractDataSet
from kedro.io.core import get_filepath_str, get_protocol_and_path
from programs.py3plex.py3plex.core import multinet

class MLNetworkDataSet(AbstractDataSet):
    def __init__(self, filepath: str):
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> pd.DataFrame:
        load_path = get_filepath_str(self._filepath, self._protocol)
        input_edges = []
        with self._fs.open(load_path, mode="r") as f:
            for line in f:
                input_edges.append(line.strip().split(","))

        mlnetwork = multinet.multi_layer_network(directed=False)
        mlnetwork.add_edges(input_edges, input_type='list')
        return mlnetwork

    def _save(
        self, input_data: multinet.multi_layer_network, sep: str = ","
    ) -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        with self._fs.open(save_path, mode="w") as f:
            for line in input_data.get_edges(data=True):
                two_nodes = list(chain.from_iterable(line[:-1]))
                two_nodes.append(str(line[-1]["weight"]))
                f.write(sep.join(two_nodes) + "\n")

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath, protocol=self._protocol)
