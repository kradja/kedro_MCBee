import networkx as nx
import sys
import pdb
MULTINET_PATH = "/groups/mcbee/radja/kedro_MCBee/src/programs/py3plex"
if not MULTINET_PATH in sys.path:
    print("Adding path to sys")
    sys.path.append(MULTINET_PATH)
    print(sys.path)

import pandas as pd
#from programs.py3plex.core import multinet
from programs.py3plex.py3plex.core import multinet
from kedromcbee.extras.dataset.mlnetwork_dataset import MLNetworkDataSet

def build_multilayer_network(cdhit_edges: pd.DataFrame, prokka_edges: pd.DataFrame) -> MLNetworkDataSet:
    A = multinet.multi_layer_network(directed=False)
    A.add_edges(cdhit_edges.values.tolist(), input_type='list')
    A.add_edges(prokka_edges.values.tolist(),input_type='list')
    print(A.basic_stats())
    return A

