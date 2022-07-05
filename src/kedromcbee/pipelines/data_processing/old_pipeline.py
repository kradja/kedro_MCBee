# from kedro.pipeline import Pipeline, node
# from kedro.pipeline.modular_pipeline import pipeline
import pdb
import os
from kedro.runner import SequentialRunner
from kedro.io import DataCatalog, MemoryDataSet
from kedro.pipeline import node, pipeline
from kedro.extras.datasets.biosequence import BioSequenceDataSet

from .extract import walkthrough_prokka_bins, parse_gff, parse_write_unique_fasta, comparethings

#def data_catalog():
#    return data_catalog

def run(param, catalog):
    #io = DataCatalog({
    #            'prokka_dir': MemoryDataSet()
    #        })
    inter = param['dp']['inter_data']
    data_catalog = DataCatalog({
        'prokka_dir': MemoryDataSet(data=param['dp']['prokka_dir']),
        'gff_files': MemoryDataSet(),
        'fasta_files': MemoryDataSet(),
        'inter_data': MemoryDataSet(data= inter),
        'umerged': MemoryDataSet(),#data=catalog['umerged']['filepath']),
        'merged_ids': MemoryDataSet(),
        'hypo_prot': MemoryDataSet()
        })

    wpb = node(func=walkthrough_prokka_bins,inputs='prokka_dir',outputs=['gff_files','fasta_files'])
    pgf = node(func=parse_gff,inputs='gff_files',outputs=['gff_cont','unique_prokka_id'])
    pff = node(func=parse_write_unique_fasta,inputs=['fasta_files','inter_data'],outputs=['umerged','merged_ids'])
    com = node(func=comparethings,inputs=['umerged','merged_ids','gff_cont','inter_data'],outputs='hypo_prot')
    dp_pipeline = pipeline([wpb,pgf,pff,com])
    runner = SequentialRunner()
    return runner.run(dp_pipeline,data_catalog)
