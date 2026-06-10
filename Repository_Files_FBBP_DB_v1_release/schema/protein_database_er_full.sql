-- FBBP database schema generic ER-friendly DDL

CREATE TABLE sources (
  `source_id` TEXT NOT NULL,
  `identifier` TEXT,
  `title` TEXT,
  `authors` TEXT,
  `journal_or_office` TEXT,
  `publication_date` TEXT,
  `source_file_name` TEXT,
  `source_type` TEXT,
  PRIMARY KEY (`source_id`)
);

CREATE TABLE proteins (
  `protein_id` INT PRIMARY KEY,
  `structure_unique_sequence_id` TEXT,
  `canonical_name` TEXT,
  `organism` TEXT,
  `source_sequence` TEXT,
  `description` TEXT,
  `scaffold_category` TEXT,
  `final_tested_sequence` TEXT,
  `sequence_length` INT,
  UNIQUE (`structure_unique_sequence_id`)
);

CREATE TABLE protein_identifiers (
  `identifier_id` TEXT NOT NULL,
  `protein_id` INT,
  `id_type` TEXT,
  `id_value` TEXT,
  `status` TEXT,
  PRIMARY KEY (`identifier_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`)
);

CREATE TABLE domains (
  `domain_id` TEXT NOT NULL,
  `protein_id` INT,
  `source_id` TEXT,
  `structure_unique_sequence_id` TEXT,
  `domain_name` TEXT,
  `sequence` TEXT,
  `is_engineered` BOOLEAN,
  `construct_details` TEXT,
  `scaffold_type` TEXT,
  PRIMARY KEY (`domain_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`),
  FOREIGN KEY (`source_id`) REFERENCES sources (`source_id`)
);

CREATE TABLE targets_conceptual (
  `target_concept_id` TEXT NOT NULL,
  `gene_name` TEXT,
  `protein_name` TEXT,
  PRIMARY KEY (`target_concept_id`)
);

CREATE TABLE target_species_variants (
  `target_variant_id` TEXT NOT NULL,
  `target_concept_id` TEXT,
  `species_name` TEXT,
  `uniprot_id` TEXT,
  `gene_name_species` TEXT,
  PRIMARY KEY (`target_variant_id`),
  FOREIGN KEY (`target_concept_id`) REFERENCES targets_conceptual (`target_concept_id`)
);

CREATE TABLE interactions (
  `interaction_id` TEXT NOT NULL,
  `domain_id` TEXT,
  `target_variant_id` TEXT,
  `source_id` TEXT,
  `interaction_class` TEXT,
  `is_inhibitory` BOOLEAN,
  `protein_row_id` INT,
  `structure_unique_sequence_id` TEXT,
  PRIMARY KEY (`interaction_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`),
  FOREIGN KEY (`target_variant_id`) REFERENCES target_species_variants (`target_variant_id`),
  FOREIGN KEY (`source_id`) REFERENCES sources (`source_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE affinity_data (
  `affinity_id` TEXT NOT NULL,
  `interaction_id` TEXT,
  `determination_method` TEXT,
  `value` TEXT,
  PRIMARY KEY (`affinity_id`),
  FOREIGN KEY (`interaction_id`) REFERENCES interactions (`interaction_id`)
);

CREATE TABLE digestive_assays (
  `assay_id` TEXT NOT NULL,
  `domain_id` TEXT,
  `enzyme_name` TEXT,
  `result_value` FLOAT,
  PRIMARY KEY (`assay_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`)
);

CREATE TABLE functional_annotations (
  `annotation_id` TEXT NOT NULL,
  `protein_id` INT,
  `annotation_type` TEXT,
  `annotation_value` TEXT,
  PRIMARY KEY (`annotation_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`)
);

CREATE TABLE structural_source_info (
  `structure_id` TEXT NOT NULL,
  `domain_id` TEXT,
  `structure_unique_sequence_id` TEXT,
  `pdb_id` TEXT,
  `method` TEXT,
  `resolution` FLOAT,
  `source_structure_identifier` TEXT,
  `source_structure_type` TEXT,
  `source_structure_note` TEXT,
  PRIMARY KEY (`structure_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`),
  FOREIGN KEY (`structure_unique_sequence_id`) REFERENCES proteins (`structure_unique_sequence_id`)
);

CREATE TABLE active_structure_context (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `structure_best_available_source` TEXT,
  `structure_source` TEXT,
  `structure_type` TEXT,
  `structure_locator` TEXT,
  `chain_id` TEXT,
  `resnum_start` INT,
  `resnum_end` INT,
  `chain_length` INT,
  `final_chain_assignment` TEXT,
  `complex_chain_annotation_completed` BOOLEAN,
  `pdb_path` TEXT,
  `pdb_code` TEXT,
  `chain_count` BOOLEAN,
  `structure_match_type` TEXT,
  `structure_match_detail` TEXT,
  `structure_match_identity` FLOAT,
  `structure_match_coverage` FLOAT,
  `structure_match_chain_coverage` FLOAT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`) REFERENCES proteins (`protein_id`)
);

CREATE TABLE bcell_epitope_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `BP3_Input_ID` TEXT,
  `BP3_Residue_Count` INTEGER,
  `BP3_Mean_Score` FLOAT,
  `BP3_Max_Score` FLOAT,
  `BP3_Mean_Linear` FLOAT,
  `BP3_Max_Linear` FLOAT,
  `BP3_Top20_Ratio` FLOAT,
  `BP3_BCell_Level` TEXT,
  `BP3_Confidence_Score_0_100` FLOAT,
  `BP3_Confidence_Level` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE mhci_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `mhci_total_rows` INTEGER,
  `mhci_unique_peptides` INTEGER,
  `mhci_hits_percentile_le_1pct` INTEGER,
  `mhci_hits_percentile_le_2pct` INTEGER,
  `mhci_hits_ba_ic50_lt50nM` INTEGER,
  `mhci_hits_ba_ic50_lt500nM` INTEGER,
  `mhci_hits_ba_ic50_lt5000nM` INTEGER,
  `mhci_immunogenicity_positive_count_score_gt0` INTEGER,
  `mhci_immunogenicity_positive_ratio` FLOAT,
  `mhci_best_peptide` TEXT,
  `mhci_best_allele` TEXT,
  `mhci_best_median_percentile` FLOAT,
  `mhci_best_ba_ic50_nM` FLOAT,
  `mhci_best_el_percentile` FLOAT,
  `mhci_best_ba_percentile` FLOAT,
  `mhci_best_immunogenicity_score` FLOAT,
  `mhci_binding_level` TEXT,
  `mhci_immunogenicity_level` TEXT,
  `mhci_confidence_score_0_100` FLOAT,
  `mhci_confidence_level` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE mhcii_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `mhcii_total_rows` INTEGER,
  `mhcii_unique_peptides` INTEGER,
  `mhcii_hits_percentile_le_10pct` INTEGER,
  `mhcii_hits_percentile_le_5pct` INTEGER,
  `mhcii_hits_ba_ic50_lt50nM` INTEGER,
  `mhcii_hits_ba_ic50_lt500nM` INTEGER,
  `mhcii_best_peptide` TEXT,
  `mhcii_best_allele` TEXT,
  `mhcii_best_median_percentile` FLOAT,
  `mhcii_best_ba_ic50_nM` FLOAT,
  `mhcii_best_el_percentile` FLOAT,
  `mhcii_best_ba_percentile` FLOAT,
  `mhcii_best_cd4_score` FLOAT,
  `mhcii_best_combined_score` FLOAT,
  `mhcii_binding_level` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE immunogenicity_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `overall_final_judgement` TEXT,
  `overall_confidence_score_0_100` FLOAT,
  `overall_confidence_level` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE ecoli_expression_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `EcoliExpr_screen_input_id` TEXT,
  `EcoliExpr_screen_row_index` INTEGER,
  `EcoliExpr_screen_track` TEXT,
  `EcoliExpr_soluprot_score` FLOAT,
  `EcoliExpr_netsolp_solubility` FLOAT,
  `EcoliExpr_netsolp_usability` FLOAT,
  `EcoliExpr_rp3net_score` FLOAT,
  `EcoliExpr_ensemble_score` FLOAT,
  `EcoliExpr_engineering_required` BOOLEAN,
  `EcoliExpr_out_of_domain` BOOLEAN,
  `EcoliExpr_recommended_ecoli_route` TEXT,
  `EcoliExpr_candidate_bucket` TEXT,
  `EcoliExpr_screening_note` TEXT,
  `EcoliExpr_track_reason` TEXT,
  `EcoliExpr_score_components_available` INTEGER,
  `EcoliExpr_priority_rank_within_track` FLOAT,
  `EcoliExpr_priority_rank_within_scaffold` INTEGER,
  `EcoliExpr_ensemble_score_weight_sum` FLOAT,
  `EcoliExpr_netsolp_tool_score` FLOAT,
  `EcoliExpr_confidence_score_0_100` INTEGER,
  `EcoliExpr_confidence_tier` TEXT,
  `EcoliExpr_tool_score_spread` FLOAT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE oral_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `Oral_ID` TEXT,
  `Oral_SGF_%` INTEGER,
  `Oral_SIF_%` FLOAT,
  `Oral_Pepsin_Frags` INTEGER,
  `Oral_Intestine_Frags` INTEGER,
  `Oral_Max_Frag_Len` INTEGER,
  `Oral_DeepDigest_Risk` INTEGER,
  `Oral_Uptake_Rating` TEXT,
  `Oral_Uptake_Evidence` TEXT,
  `Oral_CPP_Best_Prob` FLOAT,
  `Oral_CPP_Best_Seq` TEXT,
  `Oral_CPP_Positive_Unique_Count` INTEGER,
  `Oral_CPP_HighEff_Unique_Count` INTEGER,
  `Oral_PEPT1_Candidate_Count` INTEGER,
  `Oral_PEPT1_Dipeptide_Count` INTEGER,
  `Oral_PEPT1_Tripeptide_Count` INTEGER,
  `Oral_Rating_Strict` TEXT,
  `Oral_Rating` TEXT,
  `Oral_Confidence` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE solubility_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `QCBundleServerOfficial_protein_name` TEXT,
  `QCBundleServerOfficial_fasta_header` TEXT,
  `QCBundleServerOfficial_scaffold_category` TEXT,
  `QCBundleServerOfficial_protein_label` TEXT,
  `QCBundleServerOfficial_sequence_length` INTEGER,
  `QCBundleServerOfficial_pdb_file` TEXT,
  `QCBundleServerOfficial_pdb_path` TEXT,
  `QCBundleServerOfficial_assessment_branch` TEXT,
  `QCBundleServerOfficial_gravy_score` TEXT,
  `QCBundleServerOfficial_instability_index` FLOAT,
  `QCBundleServerOfficial_soluprot_score` FLOAT,
  `QCBundleServerOfficial_camsol_score` FLOAT,
  `QCBundleServerOfficial_deepsol_score` FLOAT,
  `QCBundleServerOfficial_a3d_avg_score` TEXT,
  `QCBundleServerOfficial_a3d_max_score` FLOAT,
  `QCBundleServerOfficial_final_assessment` TEXT,
  `QCBundleServerOfficial_protparam_interpretation` TEXT,
  `QCBundleServerOfficial_soluprot_call` TEXT,
  `QCBundleServerOfficial_soluprot_interpretation` TEXT,
  `QCBundleServerOfficial_deepsol_call` TEXT,
  `QCBundleServerOfficial_deepsol_interpretation` TEXT,
  `QCBundleServerOfficial_camsol_call` TEXT,
  `QCBundleServerOfficial_camsol_interpretation` TEXT,
  `QCBundleServerOfficial_a3d_call` TEXT,
  `QCBundleServerOfficial_a3d_interpretation` TEXT,
  `QCBundleServerOfficial_sequence_consensus_code` TEXT,
  `QCBundleServerOfficial_sequence_consensus` TEXT,
  `QCBundleServerOfficial_structure_interpretation_code` TEXT,
  `QCBundleServerOfficial_structure_interpretation` TEXT,
  `QCBundleServerOfficial_final_assessment_code` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE protparam_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `ProtParam_MW_Da` FLOAT,
  `ProtParam_pI` FLOAT,
  `ProtParam_GRAVY` FLOAT,
  `ProtParam_Instability_Index` FLOAT,
  `ProtParam_Aromaticity` FLOAT,
  `ProtParam_Ext_Coeff_Reduced` INTEGER,
  `ProtParam_Ext_Coeff_Oxidized` INTEGER,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE ted_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `TED_n_residues` FLOAT,
  `TED_n_high` FLOAT,
  `TED_n_medium` FLOAT,
  `TED_n_low` FLOAT,
  `TED_n_consensus` FLOAT,
  `TED_high_domains` TEXT,
  `TED_medium_domains` TEXT,
  `TED_low_domains` TEXT,
  `TED_searchable_domain_count` FLOAT,
  `TED_searchable_high_domain_count` FLOAT,
  `TED_searchable_medium_domain_count` FLOAT,
  `TED_domains_with_cath_hit_count` FLOAT,
  `TED_domains_without_cath_hit_count` FLOAT,
  `TED_domain_hit_coverage` FLOAT,
  `TED_mean_novelty_l2_topk` FLOAT,
  `TED_max_novelty_l2_topk` FLOAT,
  `TED_best_top1_target` TEXT,
  `TED_best_top1_cosine` FLOAT,
  `TED_best_cath_domain_id` TEXT,
  `TED_best_cath_class_code` FLOAT,
  `TED_best_cath_class_name` TEXT,
  `TED_best_cath_architecture_code` FLOAT,
  `TED_best_cath_architecture_name` TEXT,
  `TED_best_cath_topology_code` TEXT,
  `TED_best_cath_topology_name` TEXT,
  `TED_best_cath_superfamily_code` TEXT,
  `TED_best_cath_superfamily_name` TEXT,
  `TED_top_hit_targets_all` TEXT,
  `TED_top_hit_superfamilies_all` TEXT,
  `TED_top_hit_topologies_all` TEXT,
  `TED_best_domain_id` TEXT,
  `TED_best_domain_confidence` TEXT,
  `TED_best_domain_ranges` TEXT,
  `TED_best_domain_novelty_l2_topk` FLOAT,
  `TED_best_domain_top1_target` TEXT,
  `TED_best_domain_top1_cosine` FLOAT,
  `TED_best_domain_cath_superfamily_code` TEXT,
  `TED_best_domain_cath_superfamily_name` TEXT,
  `TED_domain_summary_compact` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE protrek_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `protrek_confidence` TEXT,
  `protrek_confidence_reason` TEXT,
  `protrek_official_seq2seq_band` TEXT,
  `protrek_official_seq2seq_score` FLOAT,
  `protrek_official_function_band` TEXT,
  `protrek_official_function_score` FLOAT,
  `protrek_structure_support` TEXT,
  `protrek_structure_support_detail` TEXT,
  `protrek_seq_top1_hit` TEXT,
  `protrek_seq_function_top1` TEXT,
  `protrek_str_top1_hit` TEXT,
  `protrek_str_top1_score` FLOAT,
  `protrek_str_function_top1` TEXT,
  `protrek_str_function_score` FLOAT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE foldseek_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `foldseek_status` TEXT,
  `foldseek_target` TEXT,
  `foldseek_evalue` TEXT,
  `foldseek_bits` FLOAT,
  `foldseek_alntmscore` FLOAT,
  `foldseek_qtmscore` FLOAT,
  `foldseek_ttmscore` FLOAT,
  `foldseek_lddt` FLOAT,
  `foldseek_prob` FLOAT,
  `foldseek_fident` FLOAT,
  `foldseek_alnlen` FLOAT,
  `foldseek_qcov` FLOAT,
  `foldseek_tcov` FLOAT,
  `foldseek_avg_tm` FLOAT,
  `foldseek_confidence` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE plmsearch_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `PLMSearch_query_original` TEXT,
  `PLMSearch_query_normalized` TEXT,
  `PLMSearch_hit_count` INTEGER,
  `PLMSearch_top1_target` TEXT,
  `PLMSearch_top1_score` FLOAT,
  `PLMSearch_top1_confidence_level` TEXT,
  `PLMSearch_top2_target` TEXT,
  `PLMSearch_top3_target` TEXT,
  `PLMSearch_top4_target` TEXT,
  `PLMSearch_top5_target` TEXT,
  `PLMSearch_top6_target` TEXT,
  `PLMSearch_top7_target` TEXT,
  `PLMSearch_top8_target` TEXT,
  `PLMSearch_top9_target` TEXT,
  `PLMSearch_top10_target` TEXT,
  `PLMSearch_top11_target` TEXT,
  `PLMSearch_top12_target` TEXT,
  `PLMSearch_top13_target` TEXT,
  `PLMSearch_top14_target` TEXT,
  `PLMSearch_top15_target` TEXT,
  `PLMSearch_top16_target` TEXT,
  `PLMSearch_top17_target` TEXT,
  `PLMSearch_top18_target` TEXT,
  `PLMSearch_top19_target` TEXT,
  `PLMSearch_top20_target` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE loop_annotations (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `loop_slot` INTEGER,
  `loop_id` TEXT NOT NULL,
  `loop_chain_id` TEXT,
  `loop_start` INTEGER,
  `loop_end` INTEGER,
  `loop_length` INTEGER,
  `loop_flank_pattern` TEXT,
  `loop_status` TEXT,
  `loop_label` TEXT,
  `loop_manual_decision` TEXT,
  `loop_confidence_tier` TEXT,
  `loop_candidate_scope` TEXT,
  `loop_note` TEXT,
  `loop_is_active` BOOLEAN,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE loop_flexibility_public_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `loop_id` TEXT NOT NULL,
  `loop_slot` INTEGER,
  `loop_chain_id` TEXT,
  `loop_start` FLOAT,
  `loop_end` FLOAT,
  `loop_length` FLOAT,
  `loop_status` TEXT,
  `loop_manual_decision` TEXT,
  `loop_confidence_tier` TEXT,
  `loop_candidate_scope` TEXT,
  `selected_structure_source` TEXT,
  `avg_rsa` FLOAT,
  `avg_asa` TEXT,
  `loop_constraint_score_total_norm` FLOAT,
  `itsflex_score` TEXT,
  `itsflex_class` TEXT,
  `flexscore_ensemble_diversity` FLOAT,
  `flexibility_score_static` TEXT,
  `flexibility_score_dynamic` FLOAT,
  `flexibility_overall_score` FLOAT,
  `flexibility_evidence_count` INTEGER,
  `flexibility_consensus_label` TEXT,
  `flexibility_conflict_flag` BOOLEAN,
  `flexibility_result_version` TEXT,
  `flexibility_result_timestamp` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`) REFERENCES loop_annotations (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`)
);

CREATE TABLE protein_flexibility_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `flexibility_loop_count_total` INTEGER,
  `flexibility_loop_count_primary` INTEGER,
  `flexibility_loop_count_with_dssp` INTEGER,
  `flexibility_loop_count_with_getcontact` INTEGER,
  `flexibility_loop_count_with_itsflex` INTEGER,
  `flexibility_loop_count_with_colabfold` INTEGER,
  `flexibility_avg_rsa_mean` FLOAT,
  `flexibility_avg_rsa_median` FLOAT,
  `flexibility_constraint_score_mean` FLOAT,
  `flexibility_constraint_score_median` FLOAT,
  `flexibility_itsflex_score_mean` TEXT,
  `flexibility_itsflex_score_median` TEXT,
  `flexibility_dynamic_score_mean` FLOAT,
  `flexibility_dynamic_score_median` FLOAT,
  `flexibility_consensus_rigid_count` INTEGER,
  `flexibility_consensus_intermediate_count` INTEGER,
  `flexibility_consensus_flexible_count` INTEGER,
  `flexibility_conflict_loop_count` INTEGER,
  `flexibility_summary_version` TEXT,
  `flexibility_summary_timestamp` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE trusted_source_qc (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `qc_review_mode` TEXT,
  `trusted_source_group` TEXT,
  `trusted_source_sampled_for_qc` BOOLEAN,
  `trusted_source_qc_verdict` TEXT,
  `trusted_source_qc_issue_type` TEXT,
  `trusted_source_qc_sampling_stratum` TEXT,
  `trusted_source_qc_sampling_rate` FLOAT,
  `trusted_source_origin_row_index` INTEGER,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE manual_review_audit (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `review_timestamp` TEXT,
  `manual_overall_verdict` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE structure_refinement_audit (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` TEXT NOT NULL,
  `pre_refine_Final_Tested_Sequence` TEXT,
  `pre_refine_Sequence_Length` INTEGER,
  `pre_refine_Sequence_PDB_ID` TEXT,
  `pre_refine_local_pdb_exact_structure_paths` TEXT,
  `pre_refine_local_pdb_fasta_paths` TEXT,
  `pre_refine_structure_locator` TEXT,
  `pre_refine_pdb_path` TEXT,
  `pre_refine_chain_id` TEXT,
  `pre_refine_chain_count` BOOLEAN,
  `pre_refine_resnum_start` INTEGER,
  `pre_refine_resnum_end` INTEGER,
  `structure_refinement_applied` BOOLEAN,
  `structure_refinement_status` TEXT,
  `structure_refinement_notes` TEXT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
);

CREATE TABLE pipeline_versions (
  `pipeline_name` TEXT NOT NULL,
  `source_table` TEXT NOT NULL,
  `version` TEXT,
  `timestamp` TEXT,
  PRIMARY KEY (`pipeline_name`, `source_table`)
);

CREATE TABLE result_metadata (
  `key` TEXT NOT NULL,
  `value` TEXT,
  PRIMARY KEY (`key`)
);

-- Indexes
CREATE INDEX idx_protein_identifiers_protein_id ON protein_identifiers (`protein_id`);
CREATE INDEX idx_domains_protein_id ON domains (`protein_id`);
CREATE INDEX idx_domains_source_id ON domains (`source_id`);
CREATE INDEX idx_target_species_variants_target_concept_id ON target_species_variants (`target_concept_id`);
CREATE INDEX idx_interactions_domain_id ON interactions (`domain_id`);
CREATE INDEX idx_interactions_target_variant_id ON interactions (`target_variant_id`);
CREATE INDEX idx_interactions_source_id ON interactions (`source_id`);
CREATE INDEX idx_interactions_protein_row_id_structure_unique_sequence_id ON interactions (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_affinity_data_interaction_id ON affinity_data (`interaction_id`);
CREATE INDEX idx_digestive_assays_domain_id ON digestive_assays (`domain_id`);
CREATE INDEX idx_functional_annotations_protein_id ON functional_annotations (`protein_id`);
CREATE INDEX idx_structural_source_info_domain_id ON structural_source_info (`domain_id`);
CREATE INDEX idx_active_structure_context_protein_row_id ON active_structure_context (`protein_row_id`);
CREATE INDEX idx_bcell_epitope_results_protein_row_id_structure_unique_sequence_id ON bcell_epitope_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_mhci_results_protein_row_id_structure_unique_sequence_id ON mhci_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_mhcii_results_protein_row_id_structure_unique_sequence_id ON mhcii_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_immunogenicity_summary_protein_row_id_structure_unique_sequence_id ON immunogenicity_summary (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_ecoli_expression_results_protein_row_id_structure_unique_sequence_id ON ecoli_expression_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_oral_results_protein_row_id_structure_unique_sequence_id ON oral_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_solubility_results_protein_row_id_structure_unique_sequence_id ON solubility_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_protparam_results_protein_row_id_structure_unique_sequence_id ON protparam_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_ted_results_protein_row_id_structure_unique_sequence_id ON ted_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_protrek_results_protein_row_id_structure_unique_sequence_id ON protrek_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_foldseek_results_protein_row_id_structure_unique_sequence_id ON foldseek_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_plmsearch_results_protein_row_id_structure_unique_sequence_id ON plmsearch_results (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_loop_annotations_protein_row_id_structure_unique_sequence_id ON loop_annotations (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_loop_flexibility_public_summary_protein_row_id_structure_unique_sequence_id_loop_id ON loop_flexibility_public_summary (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`);
CREATE INDEX idx_loop_flexibility_public_summary_protein_structure ON loop_flexibility_public_summary (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_protein_flexibility_summary_protein_row_id_structure_unique_sequence_id ON protein_flexibility_summary (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_trusted_source_qc_protein_row_id_structure_unique_sequence_id ON trusted_source_qc (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_manual_review_audit_protein_row_id_structure_unique_sequence_id ON manual_review_audit (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_structure_refinement_audit_protein_row_id_structure_unique_sequence_id ON structure_refinement_audit (`protein_row_id`, `structure_unique_sequence_id`);
