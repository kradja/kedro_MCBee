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
    #for subset_edge in cdhit_edges:
    #    n1 = subset_edge[1]
    #    n2 = subset_edge[3]
    # try to subset the edges and remove layer info
    A = multinet.multi_layer_network(directed=False)
    A.add_edges(cdhit_edges.values.tolist(), input_type='list')
    A.add_edges(prokka_edges.values.tolist(),input_type='list')
    print(A.basic_stats())
    return A

# ledge,unodes = A.get_unique_entity()
# core_net = A.core_network
# Amat = A.get_supra_adjacency_matrix()
# Anode_list = A.node_order_in_matrix
# network_analysis(A)
# test_data = pd.read_csv(os.path.join(datasets,'testing_3multilayer.csv'))
# test_datamap = test_data.values.tolist()
# 
# dum = multinet.multi_layer_network(directed=False)
# dum.add_edges(test_datamap, input_type='list')
# print(dum.basic_stats())
# 
# mat = dum.get_supra_adjacency_matrix()
# node_list = dum.node_order_in_matrix
# pdb.set_trace()
# xx = dum.numeric_core_network
# num_mat = mat.toarray()
