import itertools
import pdb
import time
import urllib
from collections import Counter
from typing import Union

import networkx as nx
import numpy as np
import obonet
import pandas as pd
from Bio import SeqIO
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.extras.datasets.json import JSONDataSet
from kedro.extras.datasets.networkx import JSONDataSet as nxJSONDataSet
from kedro.io import PartitionedDataSet
from SPARQLWrapper import JSON, SPARQLWrapper

pd.set_option("display.max_rows", 500)
pd.set_option("display.max_columns", 500)


class ReturnDict(dict):
    def __missing__(self, key):
        return key


def _parse_annotation_id(x):
    """Removing identifier for each identification.
    ex: UniProtKB:P25762; parsing
    """
    return x.str.split(":").str[1]


def _parse_protein_length(prot_sequences, unique_gids):
    prot_sequences = SeqIO.to_dict(prot_sequences)
    return {gid: len(prot_sequences[gid]) for gid in unique_gids}


def _merge_duplicate_info(df):
    """This function merges duplicate information in a dataframe and then merges based on
    a comma. It had some hardcoded column for index and groupby
    """
    unique_df = df[~df.gid.duplicated(keep=False)]
    unique_df = unique_df.set_index("gid")
    duplicated_df = df[df.gid.duplicated(keep=False)]
    tmp = duplicated_df.groupby("gid").agg(set)
    res = [[",".join(map(str, x)) for x in tmp[col]] for col in tmp.columns]
    joined_duplicates = pd.DataFrame(res).T.set_index(tmp.index)
    joined_duplicates.columns = tmp.columns
    final = pd.concat([unique_df, joined_duplicates], axis=0)
    return final


def _creating_edges_list_multilayer(annot, elist, edge, prokka_gff):
    """creating edges for multilayer networks using the list method
    [node1 layer1 node2 layer2 weight]

    What if a node comes from multiple bins. I need to check
    if both nodes have / then I am not accounting the types for both of them!
    """
    tmp = []
    edge_length = []
    edge_scaf = set()
    for node in edge:
        length, level, scaf_level = prokka_gff.loc[node][["length", "bin", "scaffold"]]
        # tmp.append(f"{node}|{level}")
        if "," in length:
            length = np.mean(np.array(length.split(",")).astype(np.int32))
        else:
            length = int(length)
        tmp.append(node)
        tmp.append(level)
        edge_length.append(length)
        edge_scaf.add(scaf_level)
    # if "low" in tmp[0] or "low" in tmp[1]:
    #    tmp.append(1)
    # else:
    #    tmp.append(-1)
    tmp.append(np.mean(edge_length))
    if len(edge_scaf) == 1:
        tmp.append(edge_scaf.pop())
    else:
        tmp.append("")
    tmp.append(annot)
    elist.append(tmp[:])
    return elist


def _merge_duplicate_seqrecords(seqrecords):
    """Removing protein sequences that are have the description and annotation of hypothetical protein
    If I merge two sequences that are the same, I need to merge the ids as well
    But if I merge those ids what will I do with the Nosema annotation for the gene?
    """
    unique_records = {}
    for rec in seqrecords:
        if "hypothetical protein" in rec.description:
            continue
        if rec.seq in unique_records:
            merged_rec = unique_records[rec.seq].id + "/" + rec.id
            unique_records[rec.seq].id = merged_rec
        else:
            unique_records[rec.seq] = rec
    return unique_records.values()


def _find_merged_ids(records):
    list_merged_id = [x.id for x in records if "/" in x.id]
    dict_merged_ids = {}
    for merged_id in list_merged_id:
        for original_id in merged_id.split("/"):
            dict_merged_ids[original_id] = merged_id
    return dict_merged_ids


def preprocess_prokka_sequences(
    partition_prokka_faa: PartitionedDataSet,
) -> Union[BioSequenceDataSet, JSONDataSet]:
    """Merging duplicate protein sequence from various partitions of prokka.

    partition_prokka_faa: A list of biopython sequence records from various partitions

    Returns:
    1. Biosequence dataset which loads fasta file
    2. JSON file of merged ids
    """
    start = time.time()
    combine_all = []
    for partition_load_func in partition_prokka_faa.values():
        partition_data = partition_load_func()
        combine_all.extend(partition_data)
    print(f"For loop {time.time() - start}")
    unique_records = _merge_duplicate_seqrecords(combine_all)
    merged_ids = _find_merged_ids(unique_records)
    return list(unique_records), merged_ids


def prokka_bins_gff(
    partition_prokka_gff: PartitionedDataSet,
) -> Union[pd.DataFrame, pd.DataFrame, JSONDataSet]:
    """Parsing prokka annotation information from the multiple runs
    I don't need to keep most of the information like gene or annotation info
    """
    list_gff_df = []
    for partition_load_func in partition_prokka_gff.values():
        list_gff_df.append(partition_load_func())
    concat_gff = pd.concat(list_gff_df, ignore_index=True)
    # concat_gff["scaf_level"] = concat_gff["scaffold"] + ":" + concat_gff["level"]
    prokka_bins = dict(concat_gff[["prokka_unique", "bin"]].values)
    prokka_gff = concat_gff.drop(["prokka_unique"], axis=1)

    uni_prokka_gff = prokka_gff[prokka_gff.annot.str.contains("UniProtKB")]
    rest_prokka_gff = prokka_gff[~prokka_gff.annot.str.contains("UniProtKB")]
    return uni_prokka_gff, rest_prokka_gff, prokka_bins


def preprocess_prokka_annotations(
    prot_sequences: BioSequenceDataSet,
    merged_gene_ids: JSONDataSet,
    annotations: pd.DataFrame,
) -> pd.DataFrame:
    """Preprocesses data for prokka annotations pulled from the gff files

    prot_sequences: A BioSequence dataset containing protein sequences
    merged_gene_ids: A JSON file containing a mapping of duplicated gene ids to a merged id
    annotations: A pandas dataframe containing information pulled from prokka's gff files

    Returns:
    A pandas dataframe containing processed prokka annotations, including gene id,
    """
    merged_gene_ids = ReturnDict(merged_gene_ids)
    annotations["gid"] = annotations.gid.map(merged_gene_ids)  # Duplicates
    # Why do a lot of the merged ids not show up in annnotations
    # pdb.set_trace()
    # all_bins = set(annotations.bin)
    # fixed_annotations = _merge_duplicate_info(annotations)
    annotations["length"] = annotations.gid.map(
        _parse_protein_length(prot_sequences, set(annotations.gid))
    )
    annotations = annotations[annotations.length > 100].copy()
    annotations["annot"] = _parse_annotation_id(annotations["annot"])
    # fixed_annotations["test"] = range(len(fixed_annotations))
    # fixed_annotations["test"] = fixed_annotations.test.astype(str)
    return annotations


def _sparql_uniprot_query(sparql, prot_list: list[str]):
    # Filter out the uniprot keywords
    prot_list = '("' + '") ("'.join(prot_list) + '")'
    sparql.setQuery(
        f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    SELECT DISTINCT
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
    res = []
    for i in range(4):
        try:
            res = sparql.query().convert()
            if res != []:
                break
        except urllib.error.HTTPError as exec:
            print(f"Sleeping for 30 seconds try {i+1}")
            print(f"HTTPError {exec}")
            time.sleep(30)
            continue
            # res = sparql.query().convert()
    return res


def _chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def go_annotations(prokka_gff: pd.DataFrame) -> pd.DataFrame:
    """I need to add asyncio here so that I can send all the requests all at once instead of waiting for the
    computation to be done and then requesting again
    This will improve the speed of performance
    """
    sparql = SPARQLWrapper("https://sparql.uniprot.org/sparql/")
    sparql.setReturnFormat(JSON)
    uni_annots = set(prokka_gff.annot)
    xres = []
    for annot_bin in _chunks(list(uni_annots), 150):
        res = _sparql_uniprot_query(sparql, annot_bin)
        if res == []:
            print("Broken, skipping this bin")
            print(annot_bin)
            continue
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


def prokka_edges(prokka_gff: pd.DataFrame) -> Union[pd.DataFrame, pd.DataFrame]:
    """I need to add the nosema intensity here. I need to implement multiprocessing here. It is taking wayyyyy too long
    I need to also take into account what's being not included in prokka_gff
    If two proteins are in the same bin with the same annotation then they are the same!
    """
    prokka_gff = prokka_gff.fillna("")
    prokka_gff["length"] = prokka_gff.length.astype(str)
    prokka_edges = []
    # prokka_gff["bin"] = ["b" + str(b) for b in prokka_gff.bin]
    # prokka_gff["gid"] = ["g" + str(g) for g in prokka_gff.gid]
    merged_prokka_gff_rows = []
    annot_groups = prokka_gff.groupby(["annot", "bin"]).gid.agg(list)
    # Merging genes whose annotation and bin are the same, thus they are effectively the same
    for annot_bin, gids in annot_groups.items():
        if len(gids) > 1:
            if "CHMOIKBA_00524" in gids:
                print(prokka_gff[prokka_gff.gid.str.contains('CHMOIKBA_00535')])
                # What's happening is that there is a gid merge in high.012 which also affects
                # low.015. This
                pdb.set_trace()
            # annot_groups[annot_bin] = "/".join(gids)
            #tmp = prokka_gff[prokka_gff.gid.isin(gids)].index
            tmp = prokka_gff[
                (prokka_gff.annot == annot_bin[0]) & (prokka_gff.bin == annot_bin[1])
            ].index
            prokka_gff.loc[tmp,'gid'] = "/".join(gids)
            # From preprocess 
            #merged_gene_ids = ReturnDict(merged_gene_ids)
            #annotations["gid"] = annotations.gid.map(merged_gene_ids)  # Duplicates
            # merged_prokka_gff_rows.extend(tmp.index.tolist())
            # res = ["/".join(set(tmp[col])) for col in tmp.columns]
            # prokka_gff.loc[len(prokka_gff)] = res
        #else:
        #    annot_groups[annot_bin] = gids[0]
    annot_groups = prokka_gff.groupby(['annot','bin']).gid.apply(lambda x: '/'.join(set(x)))
    # Removing the merged rows
    # prokka_gff = prokka_gff.drop(merged_prokka_gff_rows, axis=0)
    # for duplicate_rows in merged_prokka_gff_rows:
    #    tmp = prokka_gff.loc[duplicate_rows]
    #    res = [",".join(set(tmp[col])) for col in tmp.columns]
    #    # prokka_gff.drop(duplicate_rows, axis=0)
    #    prokka_gff.loc[len(prokka_gff)] = res

    id_bin = prokka_gff.bin.unique().astype('U')
    # bin_types_dict = {} # dict.fromkeys(bin_types,[])
    xx = prokka_gff.groupby("gid")[["scaffold", "go"]].agg(set)
    # go columns is always the same. I must have sorted it
    res = [[",".join(set(x)) for x in xx[col]] for col in xx.columns]
    node_features = pd.DataFrame(res).T.set_index(xx.index)
    node_features.columns = xx.columns
    node_list = list(node_features.index)
    node_list.extend(list(id_bin))
    node_id_conversion = dict(zip(node_list,range(len(node_list))))
    # Creating the gene-gene edges and bin-gene edges
    bin_types = set([x[:x.rfind('.')] for x in id_bin])
    for btypes in bin_types:
        # bin_types_dict[btypes] = 
        sbin = id_bin[np.where(np.char.find(id_bin,btypes) == 0)]
        # convert from bin to id
        convert_sbin = np.vectorize(node_id_conversion.__getitem__)(sbin)
        edges = list(itertools.combinations(convert_sbin,2))
        prokka_edges.extend(edges)
    annot_groups2 = annot_groups.map(node_id_conversion)#.astype(int)
    pdb.set_trace()
    annot_conversion = annot_groups.droplevel(0)
    annot_conversion.index = annot_conversion.index.map(node_id_conversion).astype(int)
    bin_genes = list(annot_conversion.items())
    prokka_edges.extend(bin_genes)
    xp = annot_groups.droplevel(1)
    values = xp.groupby('annot').agg(list)
    # annot_groups.index.get_level_values(1).map(node_id_conversion)
    #for uniprot_annots in set(annot_groups.index.get_level_values("annot")):
        #vals = annot_groups[uniprot_annots].to_list()
    for vals in values:
        # bin_gene = list(annot_groups[uniprot_annots].items())
        # prokka_edges.extend(bin_gene)
        if len(vals) == 1:
            continue
        elif len(vals) == 2:
            prokka_edges.append(tuple(vals))
        else:
            edges = list(itertools.combinations(vals, 2))
            prokka_edges.extend(edges)
    pdb.set_trace()
    return (
        pd.DataFrame(
            prokka_edges,
            columns=["node1", "node2"],  # ,
        ),
        prokka_gff, node_features
    )


def _flat_func(x):
    return [item for sublist in x for item in sublist]


def go_ontology(go_prokka: pd.DataFrame) -> nxJSONDataSet:
    # Read the taxrank ontology
    url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
    graph = obonet.read_obo(url)
    mapping = dict(
        zip(graph.nodes.keys(), [x.replace(":", "_") for x in graph.nodes.keys()])
    )
    graph = nx.relabel_nodes(graph, mapping)
    id_to_name = {id_: data["name"] for id_, data in graph.nodes(data=True)}
    # nx.descendants(graph, "GO:0000036")

    # Analysis of go_prokka
    goh_prokka = go_prokka[~go_prokka.go.isnull()].copy()
    goh_prokka["go"] = goh_prokka.go.str.split(",")
    res = goh_prokka.groupby("annot").go.agg(list)
    # dictionary of uniprot mapped to go
    flat_res = {key: set(_flat_func(val)) for key, val in res.items()}
    df = (
        pd.DataFrame.from_dict(flat_res, orient="index")
        .rename_axis("uniprots")
        .reset_index()
    )
    go_df = df.melt(id_vars=["uniprots"], value_name="go_term").dropna()
    # go mapped to uniprot
    go_annots = go_df.groupby("go_term").uniprots.agg(set)
    total_annots = [x for x in go_annots.values]
    shared_go = Counter(_flat_func(total_annots))
    go_shared = {}
    for key, val in shared_go.items():
        go_shared[val] = go_shared.get(val, []) + [key]
    shared_go_counts = Counter(list(shared_go.values()))
    count_lengths = Counter([len(x) for x in flat_res.values()])
    go_annots = go_annots.reset_index()
    go_annots["name"] = go_annots.go_term.map(id_to_name)
    go_annots.uniprots = go_annots.uniprots.map(list)
    sources = [
        x for x in graph.nodes() if graph.out_degree(x) == 1 and graph.in_degree(x) == 0
    ]
    sour = set(go_annots.go_term) & set(sources)
    print(shared_go_counts)
    print(count_lengths)
    print(sour)
    xx = go_annots.go_term.to_list()
    sub = graph.subgraph(xx)
    return (
        dict(zip(go_df.uniprots, go_df.go_term)),
        dict(zip(go_annots.go_term, go_annots.uniprots)),
        sub,
    )
