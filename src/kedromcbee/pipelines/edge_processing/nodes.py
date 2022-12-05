import itertools
import pdb
from kedro.extras.datasets.json import JSONDataSet
import pandas as pd
import numpy as np

class ReturnDict(dict):
    def __missing__(self,key):
        return key

def _creating_edges_list_multilayer(elist, edge, prokka_bins):
    """creating edges for multilayer networks using the list method
    [node1 layer1 node2 layer2 weight]

    What if a node comes from multiple bins. I need to check
    if both nodes have / then I am not accounting the types for both of them!
    """
    tmp = []
    posit = {}
    for node in edge:
        if "/" in node:
            exp_types = set()
            for dup_node in node.split("/"):
                unique_id = dup_node[: dup_node.find("_")]
                exp_types.add(prokka_bins[unique_id])
            exp_types = list(exp_types)
            tmp.append(node)
            posit[len(tmp)] = exp_types
            tmp.append(exp_types[0])
        else:
            unique_id = node[: node.find("_")]
            exp_type = prokka_bins[unique_id]
            tmp.append(node)
            tmp.append(exp_type)
    if "low" in tmp:
        tmp.append(1)
    else:
        tmp.append(-1)
    elist.append(tmp[:])
    for x in posit: # Only in case of /
        if len(posit[x]) > 1:
            set_iter = posit[x][1:]
            for exp in set_iter:
                tmp[x] = exp
                del tmp[-1]
                if "low" in tmp:
                    tmp.append(1)
                else:
                    tmp.append(-1)
                elist.append(tmp[:])
    # else:
    #    exp_types = set()
    return elist


def filter_clusters(clusterdf: pd.DataFrame,gff_prokka: pd.DataFrame) -> [pd.DataFrame, pd.DataFrame]:
    """remove clusters that are 99% remove sequences below 100"""
    single = clusterdf[clusterdf.gene_id =='']
    cluster= clusterdf[clusterdf.gene_id !='']
    print('Running')
    remove_rows = []
    remove_rows_prokka = set()
    merged_ids = {}
    temp_merge_ids = []
    for main in set(cluster.index):
        res = cluster.loc[main]
        above_99 = res[res.percent.astype(float) >= 99]
        #if not above_99.empty:
        # IS there at least one real edge. something within the cluster below 99
        if len(res) - len(above_99) == 1: #Duplicate node with no connections, just delete it
            remove_rows.append(main)
            remove_rows_prokka.update(set(res.gene_id))
        else:
            temp_merge_ids.append(main)
            temp_merge_ids.extend(above_99.gene_id.tolist())
            merged_ids[main] = temp_merge_ids
            temp_merge_ids = []
        #if main == 'CBFCJLLA_00755':
        #    pdb.set_trace()

    # I need to see when I'm using merged_ids. It should just be a dict
    # I am not actually replacing the index values
    cluster = cluster[cluster.percent.astype(float) < 99]
    cluster = cluster.drop(remove_rows, axis=0)#.reset_index("index")
    jmerged_ids = ReturnDict({key:'/'.join(val) for key,val in merged_ids.items()})
    for ids,genes in merged_ids.items():
        #if 'CBFCJLLA_00755' in ids:
        #    pdb.set_trace()
        #joined_ids = "/".join(merged_ids[ids])
        #cdf.replace(ids,joined_ids , inplace=True)
        for gene in genes:
            gff_prokka.loc[gff_prokka.gid == gene, "gid"] = jmerged_ids[ids]
    cluster.index = cluster.index.map(jmerged_ids)
    cluster.gene_id= cluster.gene_id.map(jmerged_ids)
    jtmp = list(jmerged_ids.values())
    gp = gff_prokka.set_index('gid')
    #gp = gp.drop(list(remove_rows_prokka),axis=0)
    rr = gp.loc[jtmp].groupby('gid')
    gp.loc[jtmp,'length'] = rr['length'].agg('mean')
    info_col = gp.columns[1:].tolist()
    set_rrcols = rr[info_col].agg(set)
    for cname in info_col:
        # I'm only adding sort to pass the test
        gp.loc[jtmp,cname] = set_rrcols[cname].agg(lambda x: '|'.join(sorted(map(str,x))))
    gp = gp[~gp.index.duplicated(keep='first')]
    #gp = gp.drop_duplicates()

    # gp[gp.scaf_level.str.contains(r'\|')]
    # Done with loop
    #xx = gff_prokka.groupby('gid').agg(set)
    #prokka_gff_colnames = xx.columns.tolist()
    #for cname in prokka_gff_colnames:
    #    xx[cname] = xx[cname].agg(lambda x: '|'.join(map(str,x)))
    #print(len(cluster))
    #xx = xx.drop(remove_rows,axis=0) #I need to add index information to gff_prokka
    return cluster, gp #xx.reset_index()#,','.join(single.index)

def clustered_hypo_prot_edges(
    clusterdf: pd.DataFrame, prokka_bins: JSONDataSet
) -> pd.DataFrame:
    """Creating edges based off the clusers from CD-hit and layers from the unique prokka ID given to each run"""
    cdhit_edges = []
    cdf = clusterdf
    #prokka_bins = prokka_bins.set_index("prokka_unique")["level"].T.to_dict()
    for main in set(cdf.index):
        res = cdf.loc[[main]]
        edgelist = res.gene_id.tolist()
        if len(edgelist) > 2: # and isinstance(res, pd.DataFrame):
            edges = list(itertools.combinations(edgelist, 2))
        else:
            edges = [edgelist]
        for edge in edges:
            cdhit_edges = _creating_edges_list_multilayer(
                cdhit_edges, edge, prokka_bins
            )
    return pd.DataFrame(
        cdhit_edges, columns=["node1", "layer1", "node2", "layer2", "weight"]
    )


def prokka_annotation_edges(
    gff_prokka: pd.DataFrame, prokka_bins: JSONDataSet
) -> pd.DataFrame:
    """Creating edges based off the prokka annotations
    Cutoff based off of quartile range

    """
    prokka_edges = []
    #prokka_bins = prokka_bins.set_index("prokka_unique")["level"].T.to_dict()
    annot_df = gff_prokka[gff_prokka.ainfo != 'hypothetical protein']
    #annot_df = annot_df[annot_df.length > 300]
    tmp = annot_df.annot.reset_index()
    annot_groups = tmp.groupby('annot')['gid'].apply(list)
    # There is no zero. Write a unit test for this
    for val in annot_groups.values:
        if len(val) == 1:
            continue
        elif len(val) == 2:
            prokka_edges = _creating_edges_list_multilayer(
                prokka_edges, val, prokka_bins
            )
        else:
            edges = list(itertools.combinations(val, 2))
            # Check the length of the edges and cutoff
            for edge in edges:
                prokka_edges = _creating_edges_list_multilayer(
                    prokka_edges, edge, prokka_bins
                )
    return pd.DataFrame(
        prokka_edges, columns=["node1", "layer1", "node2", "layer2", "weight"]
    )


def binned_edges(
    gff_prokka: pd.DataFrame, prokka_bins: JSONDataSet
) -> pd.DataFrame:
    """ Creating edges for everything within a bin. This will create a connected network between bins
    """
    bin_edges = []
    gff_index = set(gff_prokka.index)
    merged_gff  = [x for x in gff_index if '/' in x]
    tmp = list(gff_index - set(merged_gff))
    bins = pd.DataFrame(tmp,columns=['gid'])
    bins['uid'] = bins.gid.str.split('_').str[0]
    groups = bins.groupby('uid').agg(list).gid.to_dict()
    #groups = {key:list(g) for key,g in itertools.groupby(tmp,lambda a: a.split('_')[0])}
    res = {x:[sub[:sub.find('_')] for sub in x.split('/')] for x in merged_gff}
    for merged,vals in res.items():
        for key in vals:
            #if key == 'CBFCJLLA' or key == 'CPFNGDDN':
            #    pdb.set_trace()
            groups[key].append(merged)
    #tmp = annot_df.annot.reset_index()
    #annot_groups = tmp.groupby('annot')['gid'].apply(list)
    # There is no zero Write a unit test for this
    for val in groups.values():
        if len(val) == 1:
            continue
        elif len(val) == 2:
            bin_edges = _creating_edges_list_multilayer(
                bin_edges, val, prokka_bins
            )
        else:
            edges = list(itertools.combinations(val, 2))
            # Check the length of the edges and cutoff
            for edge in edges:
                bin_edges = _creating_edges_list_multilayer(
                    bin_edges, edge, prokka_bins
                )
    return pd.DataFrame(
        bin_edges, columns=["node1", "layer1", "node2", "layer2", "weight"]
    )
