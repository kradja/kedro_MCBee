import os
import subprocess
import pdb
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
    nodes = A.core_network.nodes()
    edges = A.core_network.edges()
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
                cmd = f"/groups/mcbee/radja/kedro_MCBee/programs/cdhit/cd-hit -i {hypo_prot_file} -c {seq_id} -aL {seq_id} -d 0 -o {cdhit_file}"
                _ = subprocess.call(cmd, shell=True)

class MLNetworkHooks:
    @hook_impl
    def after_node_run(
        self, node: Node, catalog: DataCatalog, outputs: Dict[str, Any]
    ) -> None:
        if node.short_name == "build_mlnetwork":
            input_net = outputs["bee_mlnetwork"]
            # What do I want to do? What do I want to check?
            pdb.set_trace()
            mlnetwork = outputs["bee_mlnetwork"]
            print(mlnetwork.summary())
            node_nbhr, ones, oddres, nodes, edges = network_analysis(mlnetwork)
            adj = mlnetwork.core_network.edges._adjdict
            say_hello(node)

class ProjectHooks:
    @hook_impl
    def before_node_run(self, node: Node):
        # adding extra behaviour to a single node
        if node.short_name == "prokka_bins_gff_node":
            say_hello(node)
