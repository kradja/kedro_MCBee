import os
import pdb
import subprocess
from collections import Counter, OrderedDict
from typing import Any, Dict

import networkx as nx
import numpy as np
import pandas as pd
from kedro.framework.hooks import hook_impl
from kedro.io import DataCatalog
from kedro.pipeline.node import Node


def get_best_distribution(data):
    dist_names = [
        "norm",
        "exponweib",
        "weibull_max",
        "weibull_min",
        "pareto",
        "genextreme",
    ]
    dist_results = []
    params = {}
    for dist_name in dist_names:
        dist = getattr(st, dist_name)
        param = dist.fit(data)
        params[dist_name] = param
        # Applying the Kolmogorov-Smirnov test
        D, p = st.kstest(data, dist_name, args=param)
        print("p value for " + dist_name + " = " + str(p))
        dist_results.append((dist_name, p))

    # select the best fitted distribution
    best_dist, best_p = max(dist_results, key=lambda item: item[1])
    # store the name of the best fit and its p value

    print("Best fitting distribution: " + str(best_dist))
    print("Best p value: " + str(best_p))
    print("Parameters for the best fit: " + str(params[best_dist]))

    return best_dist, best_p, params[best_dist]


def say_hello(node: Node):
    """An extra behaviour for a node to say hello before running."""
    print(f"Hello from {node.name}")


def network_analysis(A):
    net_deg = dict(nx.degree(A.core_network))
    unique_layers = set({n[1] for n in A.core_network.nodes()})
    # for each layer give a centrality
    nodes = A.core_network.nodes
    edges = A.core_network.edges
    print(f"Number of nodes {len(nodes)} and the number of edges {len(edges)}")
    degrees = list(net_deg.values())
    res = Counter(degrees)
    print(res)
    dres = dict(res)
    oddres = OrderedDict(sorted(dres.items()))
    print(oddres)
    ones = [x for x in net_deg if net_deg[x] == 1]
    node_nbhr = A.core_network.edges._nodes_nbrs()
    return node_nbhr, ones, oddres, nodes, edges


class MLNetworkHooks:
    @hook_impl
    def after_node_run(
        self, node: Node, catalog: DataCatalog, outputs: Dict[str, Any]
    ) -> None:
        if node.short_name == "not_build_mlnetwork":
            input_net = outputs["bee_mlnetwork"]
            # What do I want to do? What do I want to check?
            mlnetwork = outputs["bee_mlnetwork"]
            gff_prokka = catalog.datasets.edge_processing__merged_gff_prokka.load()
            node_nbhr, ones, oddres, nodes, edges = network_analysis(mlnetwork)
            print(gff_prokka[gff_prokka.gid == ones[3][0]])
            adj = mlnetwork.core_network.edges._adjdict
            node_list = []
            for x in nodes:
                node_list.append(x)
            df = pd.DataFrame(node_list)
            df.columns = ["id", "layer"]
            res = df[
                df.duplicated("id")
            ]  # If a merged node was in different layers. It is a duplicate
            mlnetwork.basic_stats()
            print(f"The number of edges is {len(edges)}")
            print(f"The number of nodes is {len(nodes)}")
            print(f"The number of unique nodes is {len(set(df.id))}")
            print(f"The number of duplicate nodes across layers: {len(set(res.id))}")
            C1 = mlnetwork.subnetwork(["high"], subset_by="layers")
            C2 = mlnetwork.subnetwork(["low"], subset_by="layers")
            C3 = mlnetwork.subnetwork(["control"], subset_by="layers")
            c1edge = list(C1.get_edges())
            c2edge = list(C2.get_edges())
            c3edge = list(C3.get_edges())
            print(f"Number of edges in high {len(c1edge)}")
            print(f"Number of edges in low {len(c2edge)}")
            print(f"Number of edges in control {len(c3edge)}")
            print(
                f"The number of interlayer edges {len(c1edge) +len(c2edge) + len(c3edge) }"
            )
            el = list(edges.data())
            edge_df = pd.DataFrame(el)
            edge_df.columns = ["tup1", "tup2", "info"]
            edge_df["node1"], edge_df["layer1"] = zip(*edge_df.tup1)
            edge_df["node2"], edge_df["layer2"] = zip(*edge_df.tup2)
            edgedf2 = edge_df.drop(["tup1", "tup2", "info"], axis=1)
            res2 = np.where(edgedf2.layer1 == edgedf2.layer2)[0]
            print(f"The number of intralayer edges {len(res2)}")
            # test = xx[xx.gid.str.contains("/")]['gid','scaf_level']
            gff_prokka = gff_prokka.set_index("gid")
            node_df = pd.DataFrame(nodes.data(), columns=["infos", "empty"])
            node_df["gene"], node_df["level"] = zip(*node_df.infos)
            node_df = node_df.drop(["infos", "empty", "level"], axis=1)
            node_df["length"] = node_df.gene.map(gff_prokka["length"])
            test = (
                node_df[node_df.length.str.contains(r"\|")]
                .length.str.split(r"\|")
                .apply(lambda x: np.mean(list(map(int, x))))
            )
            node_df.loc[test.index, "length"] = test
            node_length_vec = node_df.length.to_numpy()
            nlv = node_length_vec.astype(float)
            print(np.var(nlv))
            print(st.variation(nlv))  # std/mean
            print(np.std(nlv))
            print(np.mean(nlv))
            print(np.min(nlv))
            print(np.max(nlv))
            get_best_distribution(nlv)
            # variance stabilization just sqrt of x
            # tnlv = np.sqrt(nlv)
            # print(np.var(tnlv))
            # print(st.variation(tnlv)) # std/mean
            # print(np.std(tnlv))
            # print(np.mean(tnlv))
            # print(np.min(tnlv))
            # print(np.max(tnlv))
            node_df["length"] = node_df.length.astype(float)
            node_df["stable_length"] = node_df.length ** (1 / 2)
            node_df["nb_trans_length"] = node_df.length.apply(
                lambda x: np.arcsinh(np.sqrt(x + (3 / 8)) / (54904 + (6 / 8)))
            )
            node_df = node_df.drop_duplicates()

            # new_ndf = node_df.groupby('gene').length.agg(set).reset_index()
            # new_ndf = pd.merge(new_ndf,pd.DataFrame(new_ndf.length.tolist()),left_index=True,right_index=True)
            # new_ndf.columns = ['gene',"length","length2"]
            # new_ndf['stable_length'] = new_ndf.length2.astype(float)**(1/2)
            # new_ndf['nb_trans_length'] = new_ndf.length2.astype(float).apply(lambda x: np.arcsinh(np.sqrt(x+(3/8))/(54904+(6/8))))
            # I want to make sure that my I know that my duplicate genes don't have different lengths

            # But they never should have different lengths

            new_ndf = node_df.set_index("gene")
            edgedf2["n1_length"] = edgedf2.node1.map(new_ndf["length"])
            edgedf2["n1_scaf"] = edgedf2.node1.map(gff_prokka["scaf_level"])
            edgedf2["n2_length"] = edgedf2.node2.map(new_ndf["length"])
            edgedf2["n2_scaf"] = edgedf2.node2.map(gff_prokka["scaf_level"])
            G = mlnetwork.core_network
            G_degs = list(G.degree)
            degree_all_nodes = np.array([x[1] for x in G_degs])
            print(f"The minimum degree of all the nodes is {np.min(degree_all_nodes)}")
            print(np.max(degree_all_nodes))
            print(np.mean(degree_all_nodes))
            all_cc = [
                c for c in sorted(nx.connected_components(G), key=len, reverse=True)
            ]
            seem = [
                Counter(pd.DataFrame(x, columns=["gid", "layer"]).layer) for x in all_cc
            ]
            df2 = pd.DataFrame(seem)
            df2 = df2.fillna(0).astype(np.int32)
            df2["length_cc"] = df2.sum(axis=1)
            cc_layers = (df2 == 0).sum(axis=1)
            df2["num_zeroes"] = cc_layers
            print(f"There are {len(all_cc)} connected components in the graph")
            print(Counter(cc_layers.to_numpy()))
            edge_diff = edgedf2.n1_length - edgedf2.n2_length
            edge_diff = np.abs(edge_diff.to_numpy())
            large_outliers = np.where(edge_diff > np.quantile(edge_diff, 0.75))[0]
            # pdb.set_trace()
            # scaf = []
            # I need to attach the high, low, and control information to the scaffold
            #  df[df.id == "ABJEJBNN_00204/NICNEPOM_00010/JKDGBJBC_00092"]
            # this node is in 3 different levels and has three scaffolds. I need to know which scaffold is which for the
            # I only care about edges on the same layer
            # for edge in edges:
            #    for node in edge: #Removing layer info
            #        res = G_degs.loc[node[0]]
            #        #res = xx[xx.gid == node[0]]
            #        scaf.append(res.scaf_level)
            # if len(set(scaf)) == 1:
            # else:
            # scaf = []
            # node_list = []
            # say_hello(node)


class ProjectHooks:
    @hook_impl
    def before_node_run(self, node: Node):
        # adding extra behaviour to a single node
        if node.short_name == "prokka_bins_gff_node":
            say_hello(node)
