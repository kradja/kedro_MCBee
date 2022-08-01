import re
from pathlib import PurePosixPath
from typing import Any, Dict, List
import pdb

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
        self._columns = ["scaffold", "length", "gid", "annot", "ainfo", "gene", "prokka_unique", "level"]

    def _load(self) -> pd.DataFrame:
        load_path = get_filepath_str(self._filepath, self._protocol)
        gff_cont = []
        prokka_bin_unique_id = ""
        with self._fs.open(load_path, mode="r") as f:
            level = self._filepath.stem
            level = level[: level.find(".")]
            for line in f:
                if "ID=" in line:
                    unique_id = re.findall(r"ID=(.+?);", line)[0]
                    if prokka_bin_unique_id == "":
                        prokka_bin_unique_id = unique_id[: unique_id.find("_")]
                    if "UniProtKB:" in line:
                        prod_annot = re.findall(r"UniProtKB:(.+?)[;|,]", line)
                    elif "HAMAP:" in line:
                        prod_annot = re.findall(r"HAMAP:(.+?)[;|,]", line)
                    else:
                        prod_annot = ""
                    if "gene=" in line:
                        gene_id = re.findall(r"gene=(.+?);", line)[0]
                    else:
                        gene_id = ""

                    scaf = line[: line.find("\t")]
                    product = re.findall(r"product=(.+?)\n", line)
                    sline = line.split("\t")
                    prot_length = int(sline[4]) - int(sline[3]) # Protein length end - start
                    gff_cont.append(
                        [
                            scaf,
                            prot_length,
                            unique_id,
                            "/".join(prod_annot),
                            "/".join(product),
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
        prokka_gff = pd.DataFrame(
            gff_cont, columns = self._columns
        )
        return prokka_gff

    def _save(
        self, input_data: pd.DataFrame, sep: str = ","
    ) -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        input_data.to_csv(save_path, sep=sep, columns=self._columns)

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath, protocol=self._protocol)
