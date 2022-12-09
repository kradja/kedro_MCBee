import itertools
import pdb
from collections import Counter

import networkx as nx
import numpy as np
import obonet
import pandas as pd
from Bio import SeqIO
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.extras.datasets.json import JSONDataSet
from kedro.io import PartitionedDataSet
from SPARQLWrapper import JSON, SPARQLWrapper


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
    prokka_gff["gid"] = prokka_gff.gid.map(merged_ids)
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
    fres["annot"] = fres.annot.str.split(":").str[1]
    return fres


def prokka_edges(prokka_gff: pd.DataFrame, prokka_bins: JSONDataSet) -> pd.DataFrame:
    prokka_edges = []
    prokka_gff = prokka_gff[prokka_gff.length > 300].copy() # Change to 100 50% of the mean of Archae
    annot_groups = prokka_gff.annot.reset_index().groupby("annot").gid.apply(list)
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


def _sparql_uniprot_query(sparql, prot_list: list[str]):
    # Filter out the uniprot keywords
    prot_list = '("' + '") ("'.join(prot_list) + '")'
    sparql.setQuery(
        f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    SELECT
        ?ac
        ?goterm_id
    WHERE \u007b
        VALUES (?ac) \u007b{prot_list}\u007d
        BIND (IRI(CONCAT("http://purl.uniprot.org/uniprot/",?ac)) AS ?protein)
        ?protein a up:Protein .
        ?protein up:classifiedWith ?goTerm .
        BIND (REPLACE(STR(?goTerm), "^.*/([^/]*)$","$1") as ?goterm_id)
        OPTIONAL \u007b
        ?goTerm rdfs:subClassOf <http://purl.obolibrary.org/obo/GO_0003674> .
        ?goTerm rdfs:label ?moltype .
        \u007d
        FILTER(?moltype)

    \u007d ORDER BY DESC(?ac)
    """
    )
    return sparql.query().convert()


def _chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def go_annotations(prokka_gff: pd.DataFrame) -> pd.DataFrame:
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql/")
    sparql.setReturnFormat(JSON)
    uni_annots = set(prokka_gff.annot)
    xres = []
    for annot_bin in _chunks(list(uni_annots), 350):
        res = _sparql_uniprot_query(sparql, annot_bin)
        xx = res["head"]["vars"]
        rbinds = res["results"]["bindings"]
        xres.extend([[s["value"] for s in x.values()] for x in rbinds])
    df = pd.DataFrame(xres, columns=xx)
    tmp = df.groupby("ac").agg(set)
    res = [",".join(map(str, x)) for x in tmp["goterm_id"]]
    df2 = pd.DataFrame(res, index=tmp.index, columns=tmp.columns)
    dict_df = df2.goterm_id.to_dict()
    prokka_gff["go"] = prokka_gff.annot.map(dict_df)
    return prokka_gff


def _flat_func(x):
    return [item for sublist in x for item in sublist]


def go_ontology(go_prokka: pd.DataFrame) -> pd.DataFrame:
    # Read the taxrank ontology
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    graph = obonet.read_obo(url)
    id_to_name = {id_: data["name"] for id_, data in graph.nodes(data=True)}
    nx.descendants(graph, "GO:0000036")
    nx.ancestors(graph, "GO:0003677")

    # Analysis of go_prokka
    goh_prokka = go_prokka[~go_prokka.go.isnull()].copy()
    goh_prokka["go"] = goh_prokka.go.str.split(",")
    res = goh_prokka.groupby("annot").go.agg(list)
    # dictionary of uniprot mapped to go
    flat_res = {key: set(_flat_func(val)) for key, val in res.items()}
    df = pd.DataFrame.from_dict(flat_res, orient="index").reset_index()
    go_df = df.melt(id_vars=["index"], value_name="go_term").dropna()
    # go mapped to uniprot
    go_annots = go_df.groupby("go_term").index.agg(set)
    total_annots = [x for x in go_annots.values]
    shared_go = Counter(_flat_func(total_annots))
    go_shared = {}
    for key, val in shared_go.items():
        go_shared[val] = go_shared.get(val, []) + [key]
    shared_go_counts = Counter(list(shared_go.values()))
    count_lengths = Counter([len(x) for x in flat_res.values()])
    go_annots = go_annots.reset_index()
    go_annots["go_term"] = go_annots.go_term.str.replace("_", ":")
    go_annots["name"] = go_annots.go_term.map(id_to_name)
    sources = [
        x for x in graph.nodes() if graph.out_degree(x) == 1 and graph.in_degree(x) == 0
    ]
    sour = set(go_annots.go_term) & set(sources)
    pdb.set_trace()
    print("her")
    return go_prokka
