import re
from pathlib import PurePosixPath
from typing import Any, Dict, List
import pdb
import fsspec
import pandas as pd
from kedro.io import AbstractDataSet
from kedro.io.core import get_filepath_str, get_protocol_and_path
from timeit import default_timer as timer


class CDhitClusterDataSet(AbstractDataSet):
    def __init__(self, filepath: str):
        protocol, path = get_protocol_and_path(filepath)
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> pd.DataFrame:
        """I only care about the clusters. The single clusters are hypothetical proteins that aren't similar to
        anything else!!!! So I don't care about them and I don't need to keep them around
        I only ran cdhit to try and find similarities between proteins that were labelled as hypothetical
        """
        load_path = get_filepath_str(self._filepath, self._protocol)
        cluster = {}
        single = set()
        tmp = []
        main = None
        with self._fs.open(load_path, mode="r") as f:
            for line in f:
                line = line.rstrip()
                if ">Cluster" in line:
                    if tmp:
                        cluster[main] = tmp
                        tmp = []
                    continue
                elif "0\t" in line:
                    gene_id = re.search(r">(.+)(?=\.\.\.)", line).group(0)[1:]
                    single.add(gene_id)
                    #prod_length = re.search(r"(\d+)(?=aa)", line).group(0)
                    percentage = line[line.rfind(" ") + 1 :]
                    if not percentage == "*":
                        percentage = percentage[:-1] #removing %
                    else:
                        main = gene_id
                else:
                    if "1\t" in line:
                        tmp.append([gene_id,  percentage])
                        single.remove(gene_id)
                    gene_id = re.search(r">(.+)(?=\.\.\.)", line).group(0)[1:]
                    #prod_length = re.search(r"(\d+)(?=aa)", line).group(0)
                    percentage = line[line.rfind(" ") + 1 :]
                    if not percentage == "*":
                        percentage = percentage[:-1]
                    else:
                        main = gene_id
                    tmp.append([gene_id,  percentage])
        if tmp:
            cluster[main] = tmp
        #start = timer()
        #cdf = pd.concat({k: pd.Series(v) for k, v in cluster.items()})
        #cdf = cdf.droplevel(1).rename_axis("main").reset_index()
        #cdf[["gene_id", "percent"]] = pd.DataFrame(
        #    cdf[0].tolist(), index=cdf.index
        #)
        #clusterdf = cdf.drop(0, axis=1).set_index(["main", cdf.index])
        #clusterdf.percent = clusterdf.percent.replace("*", 0)
        #end = timer()
        #print(end-start)
        start = timer()
        test = pd.DataFrame.from_dict(cluster,orient='index')
        cdf = pd.melt(test.reset_index(),id_vars='index').dropna()
        cdf = cdf.set_index('index').drop('variable',axis=1).value.apply(pd.Series)
        cdf.columns = ['gene_id','percent']
        cdf.percent = cdf.percent.replace("*", 0)
        end = timer()
        print(end-start)
        xtmp = cdf.columns.tolist()
        xtmp.insert(0,cdf.index.name)
        empty = [[x,'',''] for x in single]
        em = pd.DataFrame(empty,columns=xtmp).set_index('index')
        combined = pd.concat([cdf,em])
        return combined

    def _save(
        self, input_data: pd.DataFrame, sep: str = ",", columns: List = None
    ) -> None:
        save_path = get_filepath_str(self._filepath, self._protocol)
        input_data.to_csv(save_path, sep=sep, columns=columns)

    def _describe(self) -> Dict[str, Any]:
        return dict(filepath=self._filepath, protocol=self._protocol)
