import pandas as pd
from Bio import SeqIO
import subprocess
import glob
import itertools
import pdb
import os
from kedro.extras.datasets.text import TextDataSet
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.io import PartitionedDataSet
import re

def _creating_edges_list_multilayer(elist,edge,prokka_bins):
    """ creating edges for multilayer networks using the list method
    [node1 layer1 node2 layer2 weight]

    What if a node comes from multiple bins. I need to check
    if both nodes have / then I am not accounting the types for both of them!
    """
    tmp = []
    posit = {}
    for node in edge:
        if '/' in node:
            exp_types = set()
            for dup_node in node.split('/'):
                unique_id = dup_node[:dup_node.find('_')]
                exp_types.add(prokka_bins[unique_id]) 
            exp_types = list(exp_types)
            tmp.append(node)
            posit[len(tmp)] = exp_types
            tmp.append(exp_types[0])
        else:
            unique_id = node[:node.find('_')]
            exp_type = prokka_bins[unique_id]
            tmp.append(node)
            tmp.append(exp_type)
    if 'low' in tmp: tmp.append(1)
    else: tmp.append(-1)
    elist.append(tmp[:])
    for x in posit:
        if len(posit[x]) > 1:
            set_iter = posit[x][1:]
            for exp in set_iter:
                tmp[x] = exp
                del tmp[-1]
                if 'low' in tmp: tmp.append(1)
                else: tmp.append(-1)
                elist.append(tmp[:])
    #else:
    #    exp_types = set()
    return elist

def filter_clusters(clusterdf: pd.DataFrame) -> pd.DataFrame:
    """remove clusters that are 99% remove sequences below 100
    """
    remove_rows = []
    merged_ids =[]
    temp_merge_ids = []
    for main in set(clusterdf.index.get_level_values(0)):
        res = clusterdf.loc[main]
        if float(res.length.values[0]) <= 100: # 100 nucleotide base pairs
            remove_rows.extend(res.index.tolist())
            continue
        above_99 = res[res.percent.astype(float) >= 99]
        if not above_99.empty:
            if len(res) - len(above_99) == 1:
                remove_rows.extend(res.index.tolist())
            else:
                temp_merge_ids.append(main)
                temp_merge_ids.extend(above_99.gene_id.tolist())
                joined_ids = '/'.join(temp_merge_ids)
                merged_ids.append((main,joined_ids))
                remove_rows.extend(above_99.index.tolist())
                clusterdf.rename(index={},inplace=True)
                temp_merge_ids = []

    cdf = clusterdf.drop(remove_rows,axis=0,level=1).reset_index('main')
    for ids in merged_ids:
        cdf.replace(ids[0],ids[1],inplace=True)
    return cdf

def clustered_hypo_prot_edges(clusterdf: pd.DataFrame, prokka_bins: pd.DataFrame) -> pd.DataFrame:
    """ Creating edges based off the clusers from CD-hit and layers from the unique prokka ID given to each run
    
    """
    cdhit_edges = []
    cdf = clusterdf.set_index('main')
    prokka_bins = prokka_bins.set_index('prokka_unique')['level'].T.to_dict()
    for main in set(cdf.index):
        res = cdf.loc[main]
        if len(res) > 2:
            edges = list(itertools.combinations(res.gene_id.tolist(),2))
        else:
            edges = [res.gene_id.tolist()]
        for edge in edges:
            cdhit_edges = _creating_edges_list_multilayer(cdhit_edges,edge,prokka_bins)
    #for val in cluster.values():
    #    if len(val) > 2:
    #        edges = list(itertools.combinations(val,2))
    #    else: edges = [val]

    #    for edge in edges:
    #        cdhit_edges = creating_edges_list_multilayer(cdhit_edges,edge,prokka_bin_id)

    return pd.DataFrame(cdhit_edges)

def prokka_annotation_edges(input_dict,prokka_bins):
    """ Creating edges based off the prokka annotations
        Cutoff based off of quartile range
    
    """
    prokka_edges = []
    prokka_bins = prokka_bins.set_index('prokka_unique').T.to_dict('level')
    for val in input_dict.values():
        if len(val) == 1:
            continue
        elif len(val) == 2:
            prokka_edges = _creating_edges_list_multilayer(prokka_edges,val,prokka_bins)
        else:
            edges = list(itertools.combinations(val,2))
            # Check the length of the edges and cutoff
            for edge in edges:
                prokka_edges = _creating_edges_list_multilayer(prokka_edges,edge,prokka_bins)
    return prokka_edges

#clusterdf = filter_clusters(clusterdf)
#cdhit_edges = clustered_hypo_prot_edges(clusterdf, prokka_bin_id)
#A = multinet.multi_layer_network(directed=False)
#A.add_edges(cdhit_edges, input_type='list')
#print(A.basic_stats())
#annot_df = gff_df[gff_df.ainfo != 'hypothetical protein']
#annot_groups = annot_df.groupby('annot')['gid'].apply(list)
#annot_groups = annot_groups.drop([''])
#pdb.set_trace()
#prokka_annotation_edges = prokka_annotation_edges(annot_groups.to_dict(),prokka_bin_id)
#A.add_edges(prokka_annotation_edges,input_type='list')
#ledge,unodes = A.get_unique_entity()
#core_net = A.core_network
#Amat = A.get_supra_adjacency_matrix()
#Anode_list = A.node_order_in_matrix
#network_analysis(A)

#test_data = pd.read_csv(os.path.join(datasets,'testing_3multilayer.csv'))
#test_datamap = test_data.values.tolist()

#dum = multinet.multi_layer_network(directed=False)
#dum.add_edges(test_datamap, input_type='list')
#print(dum.basic_stats())
