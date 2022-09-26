import networkx as nx
import sys
import pdb
import numpy as np
import scipy.stats as st
import pandas as pd
from collections import Counter, OrderedDict
from kedro.extras.datasets.networkx import JSONDataSet
from programs.py3plex.py3plex.core import multinet
from kedromcbee.extras.dataset.mlnetwork_dataset import MLNetworkDataSet
MULTINET_PATH = "/groups/mcbee/radja/kedro_MCBee/src/programs/py3plex"
if not MULTINET_PATH in sys.path:
    print("Adding path to sys")
    sys.path.append(MULTINET_PATH)
    print(sys.path)


def build_multilayer_network(cdhit_edges: pd.DataFrame, prokka_edges: pd.DataFrame, gff_prokka: pd.DataFrame) -> [MLNetworkDataSet, JSONDataSet]:
    gff_prokka = gff_prokka.set_index('gid')
    A = multinet.multi_layer_network(directed=False)
    A.add_edges(cdhit_edges.values.tolist(), input_type='list')
    A.add_edges(prokka_edges.values.tolist(),input_type='list')
    print(A.basic_stats())
    G = nx.Graph()
    one = list(zip(cdhit_edges.node1,cdhit_edges.layer1))
    two = list(zip(cdhit_edges.node2,cdhit_edges.layer2))
    cweights = cdhit_edges.weight.tolist()
    xx = gff_prokka[gff_prokka.length.str.contains(r'\|')].length.str.split(r'\|').apply(lambda x: np.mean(list(map(int,x))))
    gff_prokka.loc[xx.index,'length'] = xx
    gff_prokka["length"] = gff_prokka.length.astype(float)
    for row_num in range(len(cdhit_edges)):
        slice_prokka = gff_prokka.loc[[one[row_num][0],two[row_num][0]]]
        length_arr = slice_prokka.length.to_numpy()
        edge_length = np.abs(np.diff(length_arr))[0]
        check = slice_prokka.scaffold.str.contains(r'\|').any()
        if check:
            scaf_match = set.intersection(*slice_prokka.scaf_level.str.split('|').apply(set).to_list())
        else:
            set_scafs = set(slice_prokka.scaffold)
            if len(set_scafs) ==1:
                scaf_match = set_scafs
            else:
                scaf_match = ''
        G.add_edge(
            one[row_num],two[row_num],
            weight = cweights[row_num],
            annot = 'cdhit',
            length = edge_length,
            scaffold_match = scaf_match)

    one = list(zip(prokka_edges.node1,prokka_edges.layer1))
    two = list(zip(prokka_edges.node2,prokka_edges.layer2))
    pweights = prokka_edges.weight.tolist()
    for row_num in range(len(prokka_edges)):
        slice_prokka = gff_prokka.loc[[one[row_num][0],two[row_num][0]]]
        length_arr = slice_prokka.length.to_numpy()
        edge_length = np.abs(np.diff(length_arr))[0]
        check = slice_prokka.scaffold.str.contains(r'\|').any()
        if check:
            scaf_match = set.intersection(*slice_prokka.scaf_level.str.split('|').apply(set).to_list())
        else:
            set_scafs = set(slice_prokka.scaffold)
            if len(set_scafs) ==1:
                scaf_match = set_scafs
            else:
                scaf_match = ''
        G.add_edge(
            one[row_num],two[row_num],
            weight = pweights[row_num],
            annot = gff_prokka.loc[one[row_num][0]].annot,
            length = edge_length,
            scaffold_match = scaf_match)
    pdb.set_trace()
    return A,G

def analyze_networks(bee_mlnetwork: MLNetworkDataSet, bee_graph: JSONDataSet, gff_prokka: pd.DataFrame) -> pd.DataFrame:
    nodes = bee_mlnetwork.core_network.nodes
    edges = bee_mlnetwork.core_network.edges
    adj = bee_mlnetwork.core_network.edges._adjdict
    node_list = []
    for x in nodes:
        node_list.append(x)
    df = pd.DataFrame(node_list)
    df.columns = ["id","layer"]
    res = df[df.duplicated("id")] # If a merged node was in different layers. It is a duplicate
    bee_mlnetwork.basic_stats()
    print(f"The number of edges is {len(edges)}")
    print(f"The number of nodes is {len(nodes)}")
    print(f"The number of unique nodes is {len(set(df.id))}")
    print(f"The number of duplicate nodes across layers: {len(set(res.id))}")
    C1 = bee_mlnetwork.subnetwork(['high'],subset_by="layers")
    C2 = bee_mlnetwork.subnetwork(['low'],subset_by="layers")
    C3 = bee_mlnetwork.subnetwork(['control'],subset_by="layers")
    c1edge = list(C1.get_edges())
    c2edge = list(C2.get_edges())
    c3edge = list(C3.get_edges())
    print(f"Number of edges in high {len(c1edge)}")
    print(f"Number of edges in low {len(c2edge)}")
    print(f"Number of edges in control {len(c3edge)}")
    el = list(edges.data())
    edge_df = pd.DataFrame(el)
    edge_df.columns = ["tup1","tup2","info"]
    edge_df["node1"], edge_df["layer1"] = zip(*edge_df.tup1)
    edge_df["node2"], edge_df["layer2"] = zip(*edge_df.tup2)
    edgedf2 = edge_df.drop(["tup1","tup2","info"],axis=1)
    res2 = np.where(edgedf2.layer1 == edgedf2.layer2)[0]
    res3 = np.where(edgedf2.layer1 != edgedf2.layer2)[0]
    print(f"The number of intralayer edges {len(res2)}")
    tmp = len(c1edge) +len(c2edge) + len(c3edge)
    print(f"The number of interlayer edges {len(res3) }")
    #gff_prokka = gff_prokka.set_index('gid')
    #node_df = pd.DataFrame(nodes.data(),columns = ['infos','empty'])
    #node_df['gene'], node_df['level'] = zip(*node_df.infos)
    #node_df = node_df.drop(["infos","empty","level"],axis=1)
    #node_df['length'] = node_df.gene.map(gff_prokka['length'])
    #test = node_df[node_df.length.str.contains(r'\|')].length.str.split(r'\|').apply(lambda x: np.mean(list(map(int,x))))
    #node_df.loc[test.index, 'length'] = test
    #node_length_vec = node_df.length.to_numpy()
    #nlv = node_length_vec.astype(float)
    #print(np.var(nlv))
    #print(st.variation(nlv)) # std/mean
    #print(np.std(nlv))
    #print(np.mean(nlv))
    #print(np.min(nlv))
    #print(np.max(nlv))
    #node_df["length"] = node_df.length.astype(float)
    #node_df['stable_length'] = node_df.length**(1/2)
    #node_df['nb_trans_length'] = node_df.length.apply(lambda x: np.arcsinh(np.sqrt(x+(3/8))/(54904+(6/8))))
    #node_df = node_df.drop_duplicates()
    #new_ndf = node_df.set_index('gene')
    #edgedf2['n1_length'] = edgedf2.node1.map(new_ndf['length'])
    #edgedf2['n1_scaf'] = edgedf2.node1.map(gff_prokka['scaf_level'])
    #edgedf2['n2_length'] = edgedf2.node2.map(new_ndf['length'])
    #edgedf2['n2_scaf'] = edgedf2.node2.map(gff_prokka['scaf_level'])
    G = bee_mlnetwork.core_network
    G_degs = list(G.degree)
    degree_all_nodes = np.array([x[1] for x in G_degs])
    print(f"The minimum degree of all the nodes is {np.min(degree_all_nodes)}")
    print(np.max(degree_all_nodes))
    print(np.mean(degree_all_nodes))
    all_cc = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    seem = [Counter(pd.DataFrame(x,columns = ['gid','layer']).layer) for x in all_cc]
    df2 = pd.DataFrame(seem)
    df2 = df2.fillna(0).astype(np.int32)
    df2['length_cc'] = df2.sum(axis=1)
    cc_layers = (df2==0).sum(axis=1)
    df2['num_zeroes'] = cc_layers
    print(f"There are {len(all_cc)} connected components in the graph")
    print(Counter(cc_layers.to_numpy()))
    edge_diff = edgedf2.n1_length - edgedf2.n2_length
    edge_diff = np.abs(edge_diff.to_numpy())
    large_outliers = np.where(edge_diff > np.quantile(edge_diff,0.75))[0]

    cc_info = []
    for cc in all_cc:
        tt = list(bee_graph.edges(cc))
        #print(f"The number of nodes in the connected componenet is {len(cc)} and the number of edges if {len(tt)}")
        edata = [bee_graph.get_edge_data(*x)['annot'] for x in tt]
        #print(set(edata))
        cc_info.append([len(cc),len(tt),edata])
    cc_infodf = pd.DataFrame(cc_info)
    def total_num_possedges(n: int) -> int:
        return (n*(n-1))/2
    tt = list(bee_graph.edges(all_cc[1]))
    print(total_num_possedges(len(all_cc[1])))
    print(len(tt))

    pdb.set_trace()
    return edgedf2
