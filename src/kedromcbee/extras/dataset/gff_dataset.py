# import pdb
import re
from pathlib import PurePosixPath

import fsspec
import pandas as pd
from kedro.io import AbstractDataSet
from kedro.io.core import get_filepath_str, get_protocol_and_path


class GffDataSet(AbstractDataSet):
    def __init__(self, filepath: str):
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)
        self._columns = [
            "scaffold",
            "gid",
            "annot",
            "ainfo",
            "gene",
            "prokka_unique",
            "level",
        ]

    def _load(self) -> pd.DataFrame:
        load_path = get_filepath_str(self._filepath, self._protocol)
        gff_cont = []
        prokka_bin_unique_id = ""
        id_pattern = re.compile(r"ID=(.+?);")
        product_pattern = re.compile("product=(.+?)(?=\n)")
        # uniprot_pattern = re.compile(r"UniProtKB:(.+?)[;|,]")
        # hmap_pattern = re.compile(r"HAMAP:(.+?)[;|,]")
        aa_seq_pattern = re.compile(r"AA sequence:(.+?)[;|,]")
        motif_pattern = re.compile(r"motif:(.+?)[;|,]")
        gene_pattern = re.compile(r"gene=(.+?);")
        with self._fs.open(load_path, mode="r") as f:
            level = self._filepath.stem
            level = level[: level.find(".")]
            for line in f:
                if "ID=" in line:
                    unique_id = id_pattern.search(line)[1]
                    product = product_pattern.search(line)[1]
                    if product == "hypothetical protein":
                        continue
                    # bin ID
                    if prokka_bin_unique_id == "":
                        prokka_bin_unique_id = unique_id[: unique_id.find("_")]
                    # product annotation
                    # if "UniProtKB:" in line:
                    #     prod_annot = uniprot_pattern.search(line)[1]
                    # elif "HAMAP:" in line:
                    #     prod_annot = hmap_pattern.search(line)[1]
                    if "AA sequence:" in line:
                        prod_annot = aa_seq_pattern.search(line)[1]
                    elif "motif:" in line:
                        prod_annot = motif_pattern.search(line)[1]
                    else:
                        prod_annot = ""
                    # Gene name
                    if "gene=" in line:
                        gene_id = gene_pattern.search(line)[1]
                    else:
                        gene_id = ""

                    scaf = line[: line.find("\t")]
                    gff_cont.append(
                        [
                            scaf,
                            unique_id,
                            prod_annot,
                            product,
                            gene_id,
                            prokka_bin_unique_id,
                            level,
                        ]
                    )
                    prod_annot = ""
                elif line.startswith(">"):
                    break
                else:
                    continue
        prokka_gff = pd.DataFrame(gff_cont, columns=self._columns)
        return prokka_gff

    def _save(self, input_data: pd.DataFrame, sep: str = ",") -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        input_data.to_csv(save_path, sep=sep, columns=self._columns)

    def _describe(self) -> dict[str, str]:
        return dict(filepath=self._filepath, protocol=self._protocol)
