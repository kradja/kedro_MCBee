from kedromcbee.pipelines.edge_processing.nodes import filter_clusters
import pdb
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest

@pytest.fixture(scope="class")
def input_cluster_data():
    cluster_data = pd.DataFrame([
        ["g1","g1",0],
        ["g1","g3",92],
        ["g5","g4",99],
        ["g5","g5",0],
        ["g6","g7",99],
        ["g6","g8",90],
        ["g6","g6",0],
    ],columns=['index','gene_id','percent'])
    return cluster_data

@pytest.fixture(scope="class")
def input_gff_data():
    gff_data = pd.DataFrame([
        ["g1",1046,'','hypothetical protein','','scaf_1:control'],
        ["g2",346,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
        ["g3",946,'','hypothetical protein','','scaf_15:high'],
        ["g4",846,'','hypothetical protein','','scaf_7:low'],
        ["g5",746,'','hypothetical protein','','scaf_8:control'],
        ["g6",546,'','hypothetical protein','','scaf_9:control'],
        ["g7",886,'','hypothetical protein','','scaf_9:low'],
        ["g8",896,'','hypothetical protein','','scaf_2:low'],
    ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
    return gff_data

#@pytest.mark.usefixtures("cluster_data")
class Test_filter_cluster:
    def test_normal_filter_cluster(self,input_cluster_data,input_gff_data):
        normal_cluster = input_cluster_data[0:2].copy()
        normal_gff = input_gff_data[0:3].copy()

        normal_cluster = normal_cluster.set_index('index')
        output1,output2 = filter_clusters(normal_cluster,normal_gff)
        cluster_final = pd.DataFrame([
            ["g1","g1",0],
            ["g1","g3",92],
        ],columns=['index','gene_id','percent'])
        cluster_final = cluster_final.set_index('index')
        # I want to merge the rows even though I'm not using them
        # I want to delete the rows. Delete rows that I'm not interested in
        gff_final = pd.DataFrame([
            ["g1",1046.0,'','hypothetical protein','','scaf_1:control'],
            ["g2",346.0,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
            ["g3",946.0,'','hypothetical protein','','scaf_15:high'],
        ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
        assert_frame_equal(output1,cluster_final,check_dtype=False)
        gff_final = gff_final.set_index('gid')
        assert_frame_equal(output2,gff_final,check_dtype=False)

    def test_duplicate_filter_clusters(self,input_cluster_data,input_gff_data):
        duplicate_cluster = input_cluster_data[2:4].copy()
        duplicate_gff = input_gff_data.copy()

        duplicate_cluster = duplicate_cluster.set_index('index')
        output1,output2 = filter_clusters(duplicate_cluster,duplicate_gff)

        cluster_final = pd.DataFrame([
        ],columns=['index','gene_id','percent'])
        cluster_final = cluster_final.set_index('index')
        gff_final = pd.DataFrame([
            ["g1",1046.0,'','hypothetical protein','','scaf_1:control'],
            ["g2",346.0,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
            ["g3",946.0,'','hypothetical protein','','scaf_15:high'],
            ["g6",546.0,'','hypothetical protein','','scaf_9:control'],
            ["g7",886.0,'','hypothetical protein','','scaf_9:low'],
            ["g8",896.0,'','hypothetical protein','','scaf_2:low'],
        ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
        gff_final = gff_final.set_index('gid')
        assert_frame_equal(output1,cluster_final,check_dtype=False)
        assert_frame_equal(output2,gff_final,check_dtype=False)

    def test_merge_filter_clusters(self,input_cluster_data,input_gff_data):
        merge_cluster = input_cluster_data[4:7].copy()
        merge_gff = input_gff_data.copy()

        merge_cluster = merge_cluster.set_index('index')
        output1,output2 = filter_clusters(merge_cluster,merge_gff)
        cluster_final = pd.DataFrame([
            ["g6/g7","g8",90],
            ["g6/g7","g6/g7",0],
        ],columns=['index','gene_id','percent'])
        cluster_final = cluster_final.set_index('index')
        gff_final = pd.DataFrame([
            ["g1",1046.0,'','hypothetical protein','','scaf_1:control'],
            ["g2",346.0,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
            ["g3",946.0,'','hypothetical protein','','scaf_15:high'],
            ["g4",846.0,'','hypothetical protein','','scaf_7:low'],
            ["g5",746.0,'','hypothetical protein','','scaf_8:control'],
            ["g6/g7",716.0,'','hypothetical protein','','scaf_9:control|scaf_9:low'],
            ["g8",896.0,'','hypothetical protein','','scaf_2:low'],
        ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
        gff_final = gff_final.set_index('gid')
        assert_frame_equal(output1,cluster_final,check_dtype=False)
        assert_frame_equal(output2,gff_final,check_dtype=False)

    #def test_filter_clusters():
    #    # Split test into multiple functions. One function for each case
    #    cluster_data = pd.DataFrame([
    #        ["g1","g1",0],
    #        ["g1","g3",92],
    #        ["g5","g4",99],
    #        ["g5","g5",0],
    #        ["g6","g7",99],
    #        ["g6","g8",90],
    #        ["g6","g6",0],
    #    ],columns=['index','gene_id','percent'])
    #    cluster_data = cluster_data.set_index('index')
    #    gff_data = pd.DataFrame([
    #        ["g1",1046,'','hypothetical protein','','scaf_1:control'],
    #        ["g2",346,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
    #        ["g3",946,'','hypothetical protein','','scaf_15:high'],
    #        ["g4",846,'','hypothetical protein','','scaf_7:low'],
    #        ["g5",746,'','hypothetical protein','','scaf_8:control'],
    #        ["g6",546,'','hypothetical protein','','scaf_9:control'],
    #        ["g7",886,'','hypothetical protein','','scaf_9:low'],
    #        ["g8",896,'','hypothetical protein','','scaf_2:low'],
    #    ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
    #    print('hey')
    #    output1,output2 = filter_clusters(cluster_data,gff_data)
    #    assert 0
    #    cluster_final = pd.DataFrame([
    #        ["g1","g1",0],
    #        ["g1","g3",92],
    #        ["g6|g7","g6|g7",0],
    #        ["g6|g7","g8",90],
    #    ])
    #    # I want to merge the rows even though I'm not using them
    #    # I want to delete the rows. Delete rows that I'm not interested in
    #    gff_final = pd.DataFrame([
    #        ["g1",1046,'','hypothetical protein','','scaf_1:control'],
    #        ["g2",346,'MF_1','DNA_binding','xerC_1','scaf_14:low'],
    #        ["g3",946,'','hypothetical protein','','scaf_15:high'],
    #        ["g6|g7",'546|886','','hypothetical protein','','scaf_9:control|scaf_9:low'],
    #        ["g8",896,'','hypothetical protein','','scaf_2:low']
    #    ],columns = ['gid','length','annot','ainfo','gene','scaf_level'])
    #    assert output1.equals(cluster_final)
    #    assert output2.equals(gff_final)
