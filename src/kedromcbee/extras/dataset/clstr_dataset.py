from pathlib import PurePosixPath
from typing import Dict, List, Any
import fsspec
import pandas as pd
import re
from kedro.io import AbstractDataSet
from kedro.io.core import get_filepath_str, get_protocol_and_path

class CDhitClusterDataSet(AbstractDataSet):

    def __init__(self,filepath: str):
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> pd.DataFrame:
        """ I only care about the clusters. The single clusters are hypothetical proteins that aren't similar to anything else!!!!
            So I don't care about them and I don't need to keep them around

            I only ran cdhit to try and find similarities between proteins that were labelled as hypothetical
        """
        load_path = get_filepath_str(self._filepath, self._protocol)
        cluster = {}
        tmp = []
        main = None
        with self._fs.open(load_path, mode='r') as f:
            for line in f:
                line = line.rstrip()
                if '>Cluster' in line:
                    if tmp:
                        cluster[main] = tmp
                        tmp = []
                    continue
                elif '0\t' in line:
                    gene_id = re.search(r'>(.+)(?=\.\.\.)',line).group(0)[1:]
                    prod_length = re.search(r'(\d+)(?=aa)',line).group(0)
                    percentage = line[line.rfind(' ')+1:]
                    if not percentage == '*':
                        percentage = percentage[:-1]
                    else:
                        main = gene_id
                else:
                    if '1\t' in line:
                        tmp.append([gene_id,prod_length,percentage])
                    gene_id = re.search(r'>(.+)(?=\.\.\.)',line).group(0)[1:]
                    prod_length = re.search(r'(\d+)(?=aa)',line).group(0)
                    percentage = line[line.rfind(' ')+1:]
                    if not percentage == '*':
                        percentage = percentage[:-1]
                    else:
                        main = gene_id
                    tmp.append([gene_id,prod_length,percentage])
        if tmp:
            cluster[main] = tmp
        cdf = pd.concat({k: pd.Series(v) for k,v in cluster.items()})
        cdf = cdf.droplevel(1).rename_axis('main').reset_index()
        cdf[['gene_id','length','percent']]= pd.DataFrame(cdf[0].tolist(),index = cdf.index)
        clusterdf = cdf.drop(0,axis=1).set_index(['main',cdf.index])
        clusterdf.percent = clusterdf.percent.replace('*',0)
        return clusterdf

    def _save(self, input_data: pd.DataFrame, sep: str=',', columns: List=None) -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        input_data.to_csv(save_path,sep=sep, columns=columns)

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath, protocol=self._protocol)
