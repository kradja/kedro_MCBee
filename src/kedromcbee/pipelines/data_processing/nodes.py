import pandas as pd
import pdb
from Bio import SeqIO
from kedro.extras.datasets.biosequence import BioSequenceDataSet
from kedro.io import PartitionedDataSet


def prokka_bins_faa(
    partition_prokka_faa: PartitionedDataSet,
) -> [BioSequenceDataSet, pd.DataFrame]:
    unique_merge_seq = {}
    merged_ids = set()
    for fasta_file in partition_prokka_faa:
        record_list = partition_prokka_faa[fasta_file]()
        for record in record_list:
            sequence = str(record.seq)
            if sequence in unique_merge_seq:
                if unique_merge_seq[sequence].id in merged_ids:
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
    df_merged_ids = pd.DataFrame(
        dict_merged_ids.items(), columns=["single", "merged"]
    ).set_index("single")
    return list(unique_merge_seq.values()), df_merged_ids


def prokka_bins_gff(
    partition_prokka_gff: PartitionedDataSet,
) -> [pd.DataFrame, pd.DataFrame]:
    """ removing protein sequences under the length of 300 base pairs
    """
    list_gff_df = []
    for gff_file in partition_prokka_gff:
        list_gff_df.append(partition_prokka_gff[gff_file]())
    concat_gff = pd.concat(list_gff_df).reset_index(drop=True)
    prokka_bins = concat_gff[["prokka_unique", "level"]].drop_duplicates()
    # .set_index('prokka_unique')
    prokka_gff = concat_gff.drop(["prokka_unique", "level"], axis=1)
    prokka_gff["prokka_unique"] = prokka_gff.gid.apply(lambda x: x[:x.find("_")])
    prokka_gff = prokka_gff.merge(prokka_bins, on = "prokka_unique")
    prokka_gff['scaf_level'] = prokka_gff['scaffold'] + ':' +prokka_gff['level']
    return prokka_gff[prokka_gff.length > 300], prokka_bins


def hypo_prot_sequences(
    prokka_seq: BioSequenceDataSet, merged_ids: pd.DataFrame, prokka_gff: pd.DataFrame
) -> [BioSequenceDataSet, pd.DataFrame]:
    """Creates sequences of proteins annotated as hypothetical by prokka
       Would be a good idea to update prokka_gff with the merged_ids
    """

    hypo_prot = prokka_gff[prokka_gff.ainfo == "hypothetical protein"]
    hypo_prot_list = hypo_prot.gid.tolist()  # this list won't work since I merged some'
    hypo_prot_seq = {}
    prokka_seq_dict = SeqIO.to_dict(prokka_seq)
    for x in hypo_prot_list:
        if x in prokka_seq_dict:
            hypo_prot_seq[x] = prokka_seq_dict[x]
        elif x in merged_ids.index:
            hypo_prot_seq[merged_ids.loc[x].merged] = prokka_seq_dict[
                merged_ids.loc[x].merged
            ]
            prokka_gff.loc[prokka_gff.gid == x, "gid"] = merged_ids.loc[x].merged
            # merge scaffold with level?
            # merge
        else: # IF x is not in prokka_seq_dict or merged_ids
            print(x)
            continue
    xx = prokka_gff.groupby('gid').agg(set)
    prokka_gff_colnames = xx.columns.tolist()
    for cname in prokka_gff_colnames:
        xx[cname] = xx[cname].apply(lambda x: '|'.join(map(str,x)))
    return list(hypo_prot_seq.values()), xx.reset_index()
