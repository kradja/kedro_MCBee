import os
import pdb
import sys
import time
from collections import Counter, OrderedDict

import networkx as nx
import numpy as np
import pandas as pd
# import psutil
from kedro.extras.datasets.networkx import JSONDataSet

# from SPARQLWrapper import JSON, N3, SPARQLWrapper

# from programs.py3plex.py3plex.core import multinet
# from kedromcbee.extras.dataset.mlnetwork_dataset import MLNetworkDataSet
# MULTINET_PATH = "/groups/mcbee/radja/kedro_MCBee/src/programs/py3plex"
# if not MULTINET_PATH in sys.path:
#    print("Adding path to sys")
#    sys.path.append(MULTINET_PATH)
#    print(sys.path)


def _merge_node_layer(edges_df: pd.DataFrame, gff_prokka: pd.DataFrame):
    """There is no / in annotations in gff_prokka. No annotations were merged"""
    l1 = edges_df.node1.map(gff_prokka["length"])
    s1 = edges_df.node1.map(gff_prokka["scaf_level"])
    a1 = edges_df.node1.map(gff_prokka["annot"])

    l2 = edges_df.node2.map(gff_prokka["length"])
    s2 = edges_df.node2.map(gff_prokka["scaf_level"])
    a2 = edges_df.node2.map(gff_prokka["annot"])

    # start = time.time()
    # Takes a long long time
    # tmpl1 = edges_df.node1.apply(lambda x: gff_prokka.loc[x][['length','scaf_level','annot']])
    # tmpl2 = edges_df.node2.apply(lambda x: gff_prokka.loc[x][['length','scaf_level','annot']])
    # end = time.time()
    # print(end - start)

    # start = time.time()
    # one_merged = edges_df[['node1','layer1']].apply(tuple,axis=1)
    # two_merged = edges_df[['node2','layer2']].apply(tuple,axis=1)
    # end = time.time()
    # print(end - start)

    # start = time.time()
    # one_merged = pd.Series(zip(edges_df.node1,edges_df.layer1)).str.join("|")
    # two_merged = pd.Series(zip(edges_df.node2,edges_df.layer2)).str.join("|")
    # end = time.time()
    # print(end - start)

    start = time.time()
    tmp1 = zip(edges_df.node1, edges_df.layer1)
    tmp2 = zip(edges_df.node2, edges_df.layer2)
    one_merged = pd.Series(["|".join(x) for x in tmp1])
    two_merged = pd.Series(["|".join(x) for x in tmp2])
    end = time.time()
    print(end - start)

    # tmp = pd.concat([l1,s1,a1,l2,s2,a2],axis=1)
    # tmp.columns = ['length1','scaf1','annot1','length2','scaf2','annot2']
    # start = time.time()
    # edge_length = tmp[['length1','length2']].agg(np.mean,axis=1)
    # end = time.time()
    # print(f"edge length {end - start}")

    start = time.time()
    edge_length = pd.Series(np.mean(np.stack((l1, l2)), axis=0))
    end = time.time()
    print(f"faster numpy {end - start}")

    # start = time.time()
    # edge_scaf = []
    # for ind, scaf in enumerate(s1.values):
    #    edge_scaf.append(set(scaf.split('|')).intersection(set(s2.loc[2].split('|'))))

    # edge_scaf = pd.Series(edge_scaf)
    # s1 = s1.str.split('|').agg(set)
    # s2 = s2.str.split('|').agg(set)

    # end = time.time()
    # print(f"faster for loop {end - start}")

    start = time.time()
    edge_scaf = []
    edge_annot = []
    for ind, scaf in enumerate(s1.values):
        scaf2 = s2.loc[ind]
        if "/" in scaf or "/" in scaf2:
            inter = "|".join(set(scaf.split("|")).intersection(set(scaf2.split("|"))))
        elif scaf == scaf2:
            inter = scaf
        else:
            inter = ""
        edge_scaf.append(inter)
        if a1.loc[ind] == a2.loc[ind]:
            edge_annot = a1.loc[ind]
        else:
            # pdb.set_trace()
            edge_annot = ""

    edge_scaf = pd.Series(edge_scaf)
    edge_annot = pd.Series(edge_annot)
    end = time.time()
    print(f"more code faster for loop {end - start}")

    # start = time.time()
    # tmp['set_scaf1'] = tmp.scaf1.str.split('|').agg(set)
    # tmp['set_scaf2'] = tmp.scaf2.str.split('|').agg(set)
    # edge_scaf = tmp.apply(lambda x: x['set_scaf1'].intersection(x['set_scaf2']), axis=1)
    # end = time.time()
    # print(f"pandas {end - start}")

    # print(process.memory_info().rss)
    # print('after tmp contents')
    # process = psutil.Process(os.getpid())
    # print(process.memory_info().rss)

    # edge_annot = tmp.annot1 == tmp.annot2
    # edge_annot.loc[tmp[edge_annot].index] = tmp[edge_annot].annot1
    # edge_annot = edge_annot.replace(False,"")
    # edge_annot = tmp.apply(lambda x: x['annot1'].intersection(x['annot2']), axis=1)

    df = pd.concat(
        [one_merged, two_merged, edges_df.weight, edge_length, edge_scaf, edge_annot],
        axis=1,
    )
    df.columns = ["n1", "n2", "weight", "edge_length", "edge_scaf", "edge_annot"]
    # df.edge_scaf = df.edge_scaf.agg('|'.join)

    # end = time.time()
    # print(end - start)
    # print("old way")

    # start = time.time()
    # output_df = []
    # for row in edges_df.values:
    #    l1,s1,a1 = gff_prokka.loc[row[0]][['length','scaf_level','annot']].values
    #    l2,s2,a2 = gff_prokka.loc[row[2]][['length','scaf_level','annot']].values
    #    node1_merged = tuple(row[0:2])
    #    node2_merged = tuple(row[2:4])
    #    avg_length = np.mean([l1,l2])
    #    scaf_inter = set(s1.split('|')) & set(s2.split('|'))
    #    output_df.append([node1_merged,node2_merged,row[-1],avg_length,scaf_inter])
    # fdf = pd.DataFrame(output_df)
    # end = time.time()
    # print(end - start)
    return df


def build_multilayer_network(
    cdhit_edges: pd.DataFrame,
    prokka_edges: pd.DataFrame,
    bin_edges: pd.DataFrame,
    gff_prokka: pd.DataFrame,
) -> JSONDataSet:
    # gff_prokka = gff_prokka.set_index('gid')
    # xx = gff_prokka[gff_prokka.length.str.contains(r'\|')].length.str.split(r'\|').apply(lambda x: np.mean(list(map(int,x))))
    # gff_prokka.loc[xx.index,'length'] = xx
    # gff_prokka["length"] = gff_prokka.length.astype(float)
    # A = multinet.multi_layer_network(directed=False)
    # A.add_edges(cdhit_edges.values.tolist(), input_type='list')
    # A.add_edges(prokka_edges.values.tolist(),input_type='list')
    # print(A.basic_stats())
    G = nx.Graph()
    # Remove cdhit_edges2 only 3% of the total number of edges in the graph
    # We don't need cdhit since it looks like prokka does a good job
    # The hypothetical proteins that even have annotations, have some bad annotations
    # I don't care about the genes labelled hypothetical proteins. As isolated vertices they provide no information or usefulness

    prokka_edges2 = _merge_node_layer(prokka_edges, gff_prokka)
    print("Starting bin")
    # bin_edges = _merge_node_layer(bin_edges, gff_prokka)

    # graph_edges = pd.concat([cdhit_edges2, prokka_edges2])  # ,bin_edges])
    G = nx.from_pandas_edgelist(
        prokka_edges2, "n1", "n2", ["weight", "edge_length", "edge_scaf", "edge_annot"]
    )
    print(G.size())
    print(G.order())
    pdb.set_trace()
    sparql = SPARQLWrapper(
        "http://rdf.geneontology.org/blazegraph/sparql"
    )  # ../../../../data/01_raw/sparql')
    sparql.setReturnFormat(JSON)
    sparql.setQuery(
        """
    SELECT * WHERE {
        ?sub ?pred ?obj .
    }
    LIMIT 100
    OFFSET 0
    """
    )
    results = sparql.query().convert()
    res = [
        G.subgraph(x) for x in sorted(nx.connected_components(G), key=len, reverse=True)
    ]
    # print(G.get_edge_data(*list(G.edges(res[0][0]))[0]))
    # eee = list(G.edges(res[0][0]))
    # G.get_edge_data(eee[0])
    for cc in res:
        tmp = cc.nodes()
    node_cc = list(res[0].nodes())
    print(len(node_cc))
    print(len(res[0].edges()))
    print((240 * 241) / 2)
    # Every cc is strongly connected

    pdb.set_trace()
    # one = list(zip(cdhit_edges.node1,cdhit_edges.layer1))
    # two = list(zip(cdhit_edges.node2,cdhit_edges.layer2))
    # cweights = cdhit_edges.weight.tolist()
    # scaf_match = ''
    # for row_num in range(len(cdhit_edges)):
    #    slice_prokka = gff_prokka.loc[[one[row_num][0],two[row_num][0]]]
    #    length_arr = slice_prokka.length.to_numpy()
    #    edge_length = np.abs(np.diff(length_arr))[0]
    #    check = slice_prokka.scaffold.str.contains(r'\|').any()
    #    if check:
    #        scaf_match = list(set.intersection(*slice_prokka.scaf_level.str.split('|').apply(set).to_list()))
    #    else:
    #        set_scafs = set(slice_prokka.scaffold)
    #        if len(set_scafs) ==1:
    #            scaf_match = list(set_scafs)
    #        else:
    #            scaf_match = ''
    #    G.add_edge(
    #        one[row_num],two[row_num],
    #        weight = cweights[row_num],
    #        annot = 'cdhit',
    #        length = edge_length,
    #        scaffold_match = scaf_match)

    # one = list(zip(prokka_edges.node1,prokka_edges.layer1))
    # two = list(zip(prokka_edges.node2,prokka_edges.layer2))
    # pweights = prokka_edges.weight.tolist()
    # for row_num in range(len(prokka_edges)):
    #    slice_prokka = gff_prokka.loc[[one[row_num][0],two[row_num][0]]]
    #    length_arr = slice_prokka.length.to_numpy()
    #    edge_length = np.abs(np.diff(length_arr))[0]
    #    check = slice_prokka.scaffold.str.contains(r'\|').any()
    #    if check:
    #        scaf_match = list(set.intersection(*slice_prokka.scaf_level.str.split('|').apply(set).to_list()))
    #    else:
    #        set_scafs = set(slice_prokka.scaffold)
    #        if len(set_scafs) ==1:
    #            scaf_match = list(set_scafs)
    #        else:
    #            scaf_match = ''
    #    G.add_edge(
    #        one[row_num],two[row_num],
    #        weight = pweights[row_num],
    #        annot = gff_prokka.loc[one[row_num][0]].annot,
    #        length = edge_length,
    #        scaffold_match = scaf_match)
    return G


def analyze_networks(bee_graph: JSONDataSet, gff_prokka: pd.DataFrame) -> pd.DataFrame:
    nodes = bee_graph.nodes
    edges = bee_graph.edges
    adj = bee_graph.edges._adjdict
    pdb.set_trace()
    node_list = []
    for x in nodes:
        node_list.append(x)
    df = pd.DataFrame(node_list)
    df.columns = ["id", "layer"]
    res = df[
        df.duplicated("id")
    ]  # If a merged node was in different layers. It is a duplicate
    print(f"The number of edges is {len(edges)}")
    print(f"The number of nodes is {len(nodes)}")
    print(f"The number of unique nodes is {len(set(df.id))}")
    print(f"The number of duplicate nodes across layers: {len(set(res.id))}")
    # C1 = bee_mlnetwork.subnetwork(['high'],subset_by="layers")
    # C2 = bee_mlnetwork.subnetwork(['low'],subset_by="layers")
    # C3 = bee_mlnetwork.subnetwork(['control'],subset_by="layers")
    # c1edge = list(C1.get_edges())
    # c2edge = list(C2.get_edges())
    # c3edge = list(C3.get_edges())
    # print(f"Number of edges in high {len(c1edge)}")
    # print(f"Number of edges in low {len(c2edge)}")
    # print(f"Number of edges in control {len(c3edge)}")

    ## There was edge stuff here

    # gff_prokka = gff_prokka.set_index('gid')
    # node_df = pd.DataFrame(nodes.data(),columns = ['infos','empty'])
    # node_df['gene'], node_df['level'] = zip(*node_df.infos)
    # node_df = node_df.drop(["infos","empty","level"],axis=1)
    # node_df['length'] = node_df.gene.map(gff_prokka['length'])
    # test = node_df[node_df.length.str.contains(r'\|')].length.str.split(r'\|').apply(lambda x: np.mean(list(map(int,x))))
    # node_df.loc[test.index, 'length'] = test
    # node_length_vec = node_df.length.to_numpy()
    # nlv = node_length_vec.astype(float)
    # print(np.var(nlv))
    # print(st.variation(nlv)) # std/mean
    # print(np.std(nlv))
    # print(np.mean(nlv))
    # print(np.min(nlv))
    # print(np.max(nlv))
    # node_df["length"] = node_df.length.astype(float)
    # node_df['stable_length'] = node_df.length**(1/2)
    # node_df['nb_trans_length'] = node_df.length.apply(lambda x: np.arcsinh(np.sqrt(x+(3/8))/(54904+(6/8))))
    # node_df = node_df.drop_duplicates()
    # new_ndf = node_df.set_index('gene')
    # edgedf2['n1_length'] = edgedf2.node1.map(new_ndf['length'])
    # edgedf2['n1_scaf'] = edgedf2.node1.map(gff_prokka['scaf_level'])
    # edgedf2['n2_length'] = edgedf2.node2.map(new_ndf['length'])
    # edgedf2['n2_scaf'] = edgedf2.node2.map(gff_prokka['scaf_level'])
    G_degs = list(bee_graph.degree)
    degree_all_nodes = np.array([x[1] for x in G_degs])
    print(f"The minimum degree of all the nodes is {np.min(degree_all_nodes)}")
    print(np.max(degree_all_nodes))
    print(np.mean(degree_all_nodes))
    all_cc = [
        c for c in sorted(nx.connected_components(bee_graph), key=len, reverse=True)
    ]
    pdb.set_trace()
    seem = [Counter(pd.DataFrame(x, columns=["gid", "layer"]).layer) for x in all_cc]
    df2 = pd.DataFrame(seem)
    df2 = df2.fillna(0).astype(np.int32)
    df2["length_cc"] = df2.sum(axis=1)
    cc_layers = (df2 == 0).sum(axis=1)
    df2["num_zeroes"] = cc_layers
    print(f"There are {len(all_cc)} connected components in the graph")
    print(Counter(cc_layers.to_numpy()))
    # edge_diff = edgedf2.n1_length - edgedf2.n2_length
    # edge_diff = np.abs(edge_diff.to_numpy())
    # large_outliers = np.where(edge_diff > np.quantile(edge_diff,0.75))[0]

    # Adding annotation to bee_graph (prokka edges)
    pdb.set_trace()

    el = list(edges.data())
    edge_df = pd.DataFrame(el)
    edge_df.columns = ["tup1", "tup2", "info"]
    edge_df["node1"], edge_df["layer1"] = zip(*edge_df.tup1)
    edge_df["node2"], edge_df["layer2"] = zip(*edge_df.tup2)
    edgedf2 = edge_df.drop(["tup1", "tup2", "info"], axis=1)
    res2 = np.where(edgedf2.layer1 == edgedf2.layer2)[0]
    res3 = np.where(edgedf2.layer1 != edgedf2.layer2)[0]
    print(f"The number of intralayer edges {len(res2)}")
    # tmp = len(c1edge) +len(c2edge) + len(c3edge)
    print(f"The number of interlayer edges {len(res3) }")
    cc_info = []
    for cc in all_cc:
        tt = list(bee_graph.edges(cc))
        # print(f"The number of nodes in the connected componenet is {len(cc)} and the number of edges if {len(tt)}")
        edata = [bee_graph.get_edge_data(*x)["annot"] for x in tt]
        # print(set(edata))
        cc_info.append([len(cc), len(tt), edata])
    cc_infodf = pd.DataFrame(cc_info)

    def total_num_possedges(n: int) -> int:
        return (n * (n - 1)) / 2

    tt = list(bee_graph.edges(all_cc[1]))
    print(total_num_possedges(len(all_cc[1])))
    print(len(tt))

    pdb.set_trace()
    return edgedf2
