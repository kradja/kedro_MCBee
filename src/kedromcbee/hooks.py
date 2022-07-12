from kedro.framework.hooks import hook_impl
from kedro.pipeline.node import Node
from kedro.config import ConfigLoader
from kedro.io import DataCatalog
import os
import subprocess
import pdb
from typing import Dict, Any, List

def say_hello(node: Node):
    """An extra behaviour for a node to say hello before running."""
    print(f"Hello from {node.name}")
    
class DataCatalogHooks:
    @hook_impl
    def after_node_run(self, node: Node, catalog: DataCatalog, outputs: Dict[str, Any]) -> None:
        if node.short_name == 'hypo_prot_seq':
            hypo_prot_file = catalog.datasets.hypo_prot._filepath
            seq_id = catalog.datasets.params__data_processing__sequence_id._data
            cdhit_file = catalog.datasets.params__data_processing__cdhit_file._data
            cdhit_file = cdhit_file + str(seq_id)
            if not os.path.isfile(cdhit_file):
                cmd = f'/groups/mcbee/radja/kedro_MCBee/programs/cdhit/cd-hit -i {hypo_prot_file} -c {seq_id} -aL {seq_id} -d 0 -o {cdhit_file}'
                _ = subprocess.call(cmd, shell=True)

class ProjectHooks:
    @hook_impl
    def before_node_run(self, node: Node):
        # adding extra behaviour to a single node
        if node.short_name == "prokka_bins_gff_node":
            say_hello(node)
