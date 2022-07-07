from kedro.pipeline import Pipeline, node
from kedro.pipeline.modular_pipeline import pipeline

from .nodes import walkthrough_prokka_bins, parse_gff, merge_fasta_files,unique_seq_merged_ids, 
                   preprocess_companies, preprocess_shuttles, create_model_input_table, hypo_prot_sequences


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=walkthrough_prokka_bins,
                inputs=["params:prokka_dir",'params:gff_dir','params:fasta_dir'],
                outputs=['gff_files','fasta_files'],
                name="walkthrough_prokka_bins_node",
            ),
            node(
                func=parse_gff,
                inputs="gff_files",
                outputs=['prokka_gff','prokka_bins'],
                name="parse_gff_node",
            ),
            node(
                func=merge_fasta_files,
                inputs=['fasta_files','params:merged_fasta'],
                outputs=None,
                ),
            node(
                func=unique_seq_merged_ids,
                inputs='merged_fasta',
                outputs=['prokka_seq','merged_ids'],
                name="unique_seq_merged_ids_node",
            ),
            node(
                func=hypo_prot_sequences,
                inputs=['prokka_seq','merged_ids','gff_cont','inter_data'],
                outputs='hypo_prot')
            node(
                func=preprocess_companies,
                inputs="companies",
                outputs="preprocessed_companies",
                name="preprocess_companies_node",
            ),
            node(
                func=preprocess_shuttles,
                inputs="shuttles",
                outputs="preprocessed_shuttles",
                name="preprocess_shuttles_node",
            ),
            node(
                func=create_model_input_table,
                inputs=["preprocessed_shuttles", "preprocessed_companies", "reviews"],
                outputs="model_input_table",
                name="create_model_input_table_node",
            ),
        ],
        namespace="data_processing",
        inputs=["companies", "shuttles", "reviews"],
        outputs="model_input_table",
    )
