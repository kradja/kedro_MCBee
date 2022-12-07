import itertools

import numpy as np
import pandas as pd
from Bio import SeqIO
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.extras.datasets.json import JSONDataSet
from kedro.io import PartitionedDataSet


class ReturnDict(dict):
    def __missing__(self, key):
        return key


def _creating_edges_list_multilayer(annot, elist, edge, prokka_gff, prokka_bins):
    """creating edges for multilayer networks using the list method
    [node1 layer1 node2 layer2 weight]

    What if a node comes from multiple bins. I need to check
    if both nodes have / then I am not accounting the types for both of them!
    """
    tmp = []
    edge_length = []
    edge_scaf = set()
    for node in edge:
        length, level, scaf_level = prokka_gff.loc[node][
            ["length", "level", "scaf_level"]
        ]
        tmp.append(f"{node}|{level}")
        edge_length.append(length)
        edge_scaf.add(scaf_level)
    if "low" in tmp[0] or "low" in tmp[1]:
        tmp.append(1)
    else:
        tmp.append(-1)
    tmp.append(np.mean(edge_length))
    if len(edge_scaf) == 1:
        tmp.append(edge_scaf.pop())
    else:
        tmp.append("")
    tmp.append(annot)
    elist.append(tmp[:])
    return elist


def prokka_bins_faa(
    partition_prokka_faa: PartitionedDataSet,
) -> [BioSequenceDataSet, JSONDataSet]:
    """
    Reading each of the prokka faa files split into each of the maxbin clusters and merging duplicate genes
    inputs:
    outputs:
        1. All unique protein sequences
        2. Merged genes
    """
    unique_merge_seq = {}
    merged_ids = set()
    for fasta_file in partition_prokka_faa:
        record_list = partition_prokka_faa[fasta_file]()
        for record in record_list:
            sequence = str(record.seq)
            if sequence in unique_merge_seq:
                if (
                    unique_merge_seq[sequence].id in merged_ids
                ):  # Checking if sequence has already been added
                    merged_ids.remove(unique_merge_seq[sequence].id)
                unique_merge_seq[sequence].id = (
                    unique_merge_seq[sequence].id + "/" + record.id
                )
                merged_ids.add(unique_merge_seq[sequence].id)
            else:
                unique_merge_seq[sequence] = record

    dict_merged_ids = {}
    for merged_id in merged_ids:
        for original_id in merged_id.split("/"):
            dict_merged_ids[original_id] = merged_id
    # df_merged_ids = pd.DataFrame(
    #    dict_merged_ids.items(), columns=["single", "merged"]
    # ).set_index("single")
    return list(unique_merge_seq.values()), dict_merged_ids


def prokka_bins_gff(
    partition_prokka_gff: PartitionedDataSet,
) -> [pd.DataFrame, JSONDataSet]:
    """Parsing prokka annotation information from the multiple runs"""
    list_gff_df = []
    for gff_file in partition_prokka_gff:
        list_gff_df.append(partition_prokka_gff[gff_file]())
    concat_gff = pd.concat(list_gff_df).reset_index(drop=True)
    prokka_bins = dict(concat_gff[["prokka_unique", "level"]].values)
    concat_gff["scaf_level"] = concat_gff["scaffold"] + ":" + concat_gff["level"]
    prokka_gff = concat_gff.drop(["prokka_unique", "scaffold"], axis=1)
    uni_prokka_gff = prokka_gff[prokka_gff.annot.str.contains("UniProtKB")]
    rest_prokka_gff = prokka_gff[~prokka_gff.annot.str.contains("UniProtKB")]
    """I should just merge the merged rows right here!
    wasn't there some difference between the gff and fasta files?
    Why do I not remember that difference
    """
    return uni_prokka_gff, rest_prokka_gff, prokka_bins


def hypo_prot_sequences(
    prokka_seq: BioSequenceDataSet, merged_ids: JSONDataSet, prokka_gff: pd.DataFrame
) -> pd.DataFrame:
    """There are more sequences than what's in the gff file
     Creates sequences of proteins annotated as hypothetical by prokka
    Would be a good idea to update prokka_gff with the merged_ids
    """
    merged_ids = ReturnDict(merged_ids)
    prokka_seq_dict = SeqIO.to_dict(prokka_seq)
    prokka_gff.gid = prokka_gff.gid.map(merged_ids)
    prokka_dups = prokka_gff[prokka_gff.gid.duplicated(keep=False)]
    prokka_gff = prokka_gff[~prokka_gff.gid.duplicated(keep=False)]
    tmp = prokka_dups.groupby("gid").agg(set)
    res = [[",".join(map(str, x)) for x in tmp[col]] for col in tmp.columns]
    df = pd.DataFrame(res).T
    df = df.set_index(tmp.index)
    df.columns = tmp.columns
    prokka_gff = prokka_gff.set_index("gid")
    fres = pd.concat([prokka_gff, df], axis=0)
    fres["length"] = [len(prokka_seq_dict[x]) for x in fres.index]
    return fres


def prokka_edges(prokka_gff: pd.DataFrame, prokka_bins: JSONDataSet) -> pd.DataFrame:
    prokka_edges = []
    prokka_gff = prokka_gff[prokka_gff.length > 300].copy()
    prokka_gff["uni_annot"] = prokka_gff.annot.str.split(":").str[1]
    annot_groups = (
        prokka_gff.uni_annot.reset_index().groupby("uni_annot").gid.apply(list)
    )
    for uni_annot in annot_groups.keys():
        val = annot_groups[uni_annot]
        if len(val) == 1:
            continue
        elif len(val) == 2:
            prokka_edges = _creating_edges_list_multilayer(
                uni_annot, prokka_edges, val, prokka_gff, prokka_bins
            )
        else:
            edges = list(itertools.combinations(val, 2))
            for edge in edges:
                prokka_edges = _creating_edges_list_multilayer(
                    uni_annot, prokka_edges, edge, prokka_gff, prokka_bins
                )
    return pd.DataFrame(
        prokka_edges,
        columns=["node1", "node2", "weight", "edge_length", "edge_scaf", "edge_annot"],
    )
