import os
import subprocess
import pdb
import pandas as pd
import numpy as np
import networkx as nx
from typing import Any, Dict
from collections import Counter, OrderedDict

from kedro.framework.hooks import hook_impl
from kedro.io import DataCatalog
from kedro.pipeline.node import Node


def say_hello(node: Node):
    """An extra behaviour for a node to say hello before running."""
    print(f"Hello from {node.name}")

def network_analysis(A):
    net_deg = dict(nx.degree(A.core_network))
    unique_layers = set({n[1] for n in A.core_network.nodes()})
    # for each layer give a centrality
    nodes = A.core_network.nodes
    edges = A.core_network.edges
    print(f'Number of nodes {len(nodes)} and the number of edges {len(edges)}')
    degrees = list(net_deg.values())
    res = Counter(degrees)
    print(res)
    dres = dict(res)
    oddres = OrderedDict(sorted(dres.items()))
    print(oddres)
    ones = [x for x in net_deg if net_deg[x] == 1]
    node_nbhr = A.core_network.edges._nodes_nbrs()
    return node_nbhr, ones, oddres, nodes, edges

class DataCatalogHooks:
    @hook_impl
    def after_node_run(
        self, node: Node, catalog: DataCatalog, outputs: Dict[str, Any]
    ) -> None:
        if node.short_name == "hypo_prot_seq":
            hypo_prot_file = catalog.datasets.hypo_prot._filepath
            seq_id = catalog.datasets.params__data_processing__sequence_id._data
            cdhit_file = catalog.datasets.params__data_processing__cdhit_file._data
            cdhit_file = cdhit_file + str(seq_id)
            if not os.path.isfile(cdhit_file):
                cmd = f"/groups/mcbee/radja/kedro_MCBee/src/programs/cdhit/cd-hit -i {hypo_prot_file} -c {seq_id} -aL {seq_id} -d 0 -o {cdhit_file}"
                _ = subprocess.call(cmd, shell=True)

class MLNetworkHooks:
    @hook_impl
    def after_node_run(
        self, node: Node, catalog: DataCatalog, outputs: Dict[str, Any]
    ) -> None:
        if node.short_name == "build_mlnetwork":
            input_net = outputs["bee_mlnetwork"]
            # What do I want to do? What do I want to check?
            mlnetwork = outputs["bee_mlnetwork"]
            xx = catalog.datasets.edge_processing__gff_prokka.load()
            node_nbhr, ones, oddres, nodes, edges = network_analysis(mlnetwork)
            print(xx[xx.gid == ones[3][0]])
            adj = mlnetwork.core_network.edges._adjdict
            node_list = []
            for x in nodes:
                node_list.append(x)
            df = pd.DataFrame(node_list)
            df.columns = ["id","layer"]
            res = df[df.duplicated("id")] # If a merged node was in different layers. It is a duplicate
            mlnetwork.basic_stats()
            print(f"The number of edges is {len(edges)}")
            print(f"The number of nodes is {len(nodes)}")
            print(f"The number of unique nodes is {len(set(df.id))}")
            print(f"The number of duplicate nodes across layers: {len(set(res.id))}")
            C1 = mlnetwork.subnetwork(['high'],subset_by="layers")
            C2 = mlnetwork.subnetwork(['low'],subset_by="layers")
            C3 = mlnetwork.subnetwork(['control'],subset_by="layers")
            c1edge = list(C1.get_edges())
            c2edge = list(C2.get_edges())
            c3edge = list(C3.get_edges())
            print(f"Number of edges in high {len(c1edge)}")
            print(f"Number of edges in low {len(c2edge)}")
            print(f"Number of edges in control {len(c3edge)}")
            print(f"The number of interlayer edges {len(c1edge) +len(c2edge) + len(c3edge) }")
            el = list(edges.data())
            edge_df = pd.DataFrame(el)
            edge_df.columns = ["tup1","tup2","info"]
            edge_df["node1"], edge_df["layer1"] = zip(*edge_df.tup1)
            edge_df["node2"], edge_df["layer2"] = zip(*edge_df.tup2)
            edgedf2 = edge_df.drop(["tup1","tup2","info"],axis=1)
            res2 = np.where(edgedf2.layer1 == edgedf2.layer2)[0]
            print(f"The number of intralayer edges {len(res2)}")
            #test = xx[xx.gid.str.contains("/")]['gid','scaf_level']
            xx = xx.set_index('gid')
            edgedf2['n1_length'] = edgedf2.node1.map(xx['length'])
            edgedf2['n1_scaf'] = edgedf2.node1.map(xx['scaf_level'])
            edgedf2['n2_length'] = edgedf2.node2.map(xx['length'])
            edgedf2['n2_scaf'] = edgedf2.node2.map(xx['scaf_level'])
            pdb.set_trace()
            scaf = []
            # I need to attach the high, low, and control information to the scaffold
            #  df[df.id == "ABJEJBNN_00204/NICNEPOM_00010/JKDGBJBC_00092"]
            # this node is in 3 different levels and has three scaffolds. I need to know which scaffold is which for the 
            # I only care about edges on the same layer
            for edge in edges:
                for node in edge: #Removing layer info
                    res = xx.loc[node[0]]
                    #res = xx[xx.gid == node[0]]
                    scaf.append(res.scaf_level)
                #if len(set(scaf)) == 1:
                #else:
                #scaf = []
            #node_list = []
            #say_hello(node)

class ProjectHooks:
    @hook_impl
    def before_node_run(self, node: Node):
        # adding extra behaviour to a single node
        if node.short_name == "prokka_bins_gff_node":
            say_hello(node)
