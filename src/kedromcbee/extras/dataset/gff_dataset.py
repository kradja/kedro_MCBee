from pathlib import PurePosixPath
from typing import Dict, List, Any
import fsspec
import pandas as pd
import re
from kedro.io import AbstractDataSet
from kedro.io.core import get_filepath_str, get_protocol_and_path

class GffDataSet(AbstractDataSet):

    def __init__(self,filepath: str):
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> pd.DataFrame:
        load_path = get_filepath_str(self._filepath, self._protocol)
        gff_cont = []
        with self._fs.open(load_path, mode='r') as f:
            for line in f:
                if 'ID=' in line:
                    unique_id = re.findall(r'ID=(.+?);',line)[0]
                    if 'UniProtKB:' in line:
                        prod_annot = re.findall(r'UniProtKB:(.+?)[;|,]',line)
                    elif  'HAMAP:' in line:
                        prod_annot = re.findall(r'HAMAP:(.+?)[;|,]',line)
                    else:
                        prod_annot = ''
                    if 'gene=' in line:
                        gene_id = re.findall(r'gene=(.+?);',line)[0]
                    else: gene_id = ''

                    scaf = line[:line.find('\t')]
                    product = re.findall(r'product=(.+?)\n',line)
                    gff_cont.append([scaf,unique_id,'/'.join(prod_annot),'/'.join(product),gene_id])
                    prod_annot = ''
                elif line.startswith('>'):
                    break
                else:
                    continue
        prokka_gff = pd.DataFrame(gff_cont,columns=['scaffold','gid','annot','ainfo','gene'])
        return prokka_gff

    def _save(self, input_data: pd.DataFrame, sep: str=',', columns: List=None) -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        input_data.to_csv(save_path,sep=sep, columns=columns)

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath, protocol=self._protocol)
