-- FBBP database schema final
SET FOREIGN_KEY_CHECKS = 0;

CREATE TABLE sources (
  `source_id` VARCHAR(64) NOT NULL,
  `identifier` VARCHAR(38),
  `title` VARCHAR(195),
  `authors` VARCHAR(383),
  `journal_or_office` VARCHAR(161),
  `publication_date` VARCHAR(32),
  `source_file_name` VARCHAR(172),
  `source_type` VARCHAR(32),
  PRIMARY KEY (`source_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE proteins (
  `protein_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64),
  `canonical_name` VARCHAR(255),
  `organism` VARCHAR(128),
  `source_sequence` TEXT,
  `description` VARCHAR(182),
  `scaffold_category` VARCHAR(32),
  `final_tested_sequence` TEXT,
  `sequence_length` INT,
  PRIMARY KEY (`protein_id`),
  UNIQUE KEY `uk_proteins_structure_unique_sequence_id` (`structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protein_identifiers (
  `identifier_id` VARCHAR(64) NOT NULL,
  `protein_id` INT,
  `id_type` VARCHAR(32),
  `id_value` VARCHAR(32),
  `status` VARCHAR(32),
  PRIMARY KEY (`identifier_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE domains (
  `domain_id` VARCHAR(64) NOT NULL,
  `protein_id` INT,
  `source_id` VARCHAR(64),
  `structure_unique_sequence_id` VARCHAR(64),
  `domain_name` VARCHAR(85),
  `sequence` TEXT,
  `is_engineered` TINYINT(1),
  `construct_details` VARCHAR(64),
  `scaffold_type` VARCHAR(32),
  PRIMARY KEY (`domain_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`),
  FOREIGN KEY (`source_id`) REFERENCES sources (`source_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE targets_conceptual (
  `target_concept_id` VARCHAR(64) NOT NULL,
  `gene_name` VARCHAR(144),
  `protein_name` VARCHAR(34),
  PRIMARY KEY (`target_concept_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE target_species_variants (
  `target_variant_id` VARCHAR(64) NOT NULL,
  `target_concept_id` VARCHAR(64),
  `species_name` VARCHAR(32),
  `uniprot_id` VARCHAR(64),
  `gene_name_species` VARCHAR(144),
  PRIMARY KEY (`target_variant_id`),
  FOREIGN KEY (`target_concept_id`) REFERENCES targets_conceptual (`target_concept_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE interactions (
  `interaction_id` VARCHAR(64) NOT NULL,
  `domain_id` VARCHAR(64),
  `target_variant_id` VARCHAR(64),
  `source_id` VARCHAR(64),
  `interaction_class` VARCHAR(47),
  `is_inhibitory` TINYINT(1),
  `protein_row_id` INT,
  `structure_unique_sequence_id` VARCHAR(64),
  PRIMARY KEY (`interaction_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`),
  FOREIGN KEY (`target_variant_id`) REFERENCES target_species_variants (`target_variant_id`),
  FOREIGN KEY (`source_id`) REFERENCES sources (`source_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE affinity_data (
  `affinity_id` VARCHAR(64) NOT NULL,
  `interaction_id` VARCHAR(64),
  `determination_method` VARCHAR(32),
  `value` VARCHAR(32),
  PRIMARY KEY (`affinity_id`),
  FOREIGN KEY (`interaction_id`) REFERENCES interactions (`interaction_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE digestive_assays (
  `assay_id` VARCHAR(64) NOT NULL,
  `domain_id` VARCHAR(64),
  `enzyme_name` VARCHAR(44),
  `result_value` DOUBLE,
  PRIMARY KEY (`assay_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE functional_annotations (
  `annotation_id` VARCHAR(64) NOT NULL,
  `protein_id` INT,
  `annotation_type` VARCHAR(32),
  `annotation_value` TEXT,
  PRIMARY KEY (`annotation_id`),
  FOREIGN KEY (`protein_id`) REFERENCES proteins (`protein_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE structural_source_info (
  `structure_id` VARCHAR(64) NOT NULL,
  `domain_id` VARCHAR(64),
  `structure_unique_sequence_id` VARCHAR(64),
  `pdb_id` VARCHAR(64),
  `method` VARCHAR(32),
  `resolution` DOUBLE,
  `source_structure_identifier` VARCHAR(32),
  `source_structure_type` VARCHAR(32),
  `source_structure_note` VARCHAR(172),
  PRIMARY KEY (`structure_id`),
  FOREIGN KEY (`domain_id`) REFERENCES domains (`domain_id`),
  FOREIGN KEY (`structure_unique_sequence_id`) REFERENCES proteins (`structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE active_structure_context (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `structure_best_available_source` VARCHAR(32),
  `structure_source` VARCHAR(32),
  `structure_type` VARCHAR(39),
  `structure_locator` VARCHAR(32),
  `chain_id` VARCHAR(64),
  `resnum_start` INT,
  `resnum_end` INT,
  `chain_length` INT,
  `final_chain_assignment` VARCHAR(32),
  `complex_chain_annotation_completed` TINYINT(1),
  `pdb_path` VARCHAR(250),
  `pdb_code` VARCHAR(32),
  `chain_count` TINYINT(1),
  `structure_match_type` VARCHAR(32),
  `structure_match_detail` VARCHAR(32),
  `structure_match_identity` DOUBLE,
  `structure_match_coverage` DOUBLE,
  `structure_match_chain_coverage` DOUBLE,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`) REFERENCES proteins (`protein_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE bcell_epitope_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `BP3_Input_ID` VARCHAR(64),
  `BP3_Residue_Count` INT,
  `BP3_Mean_Score` DOUBLE,
  `BP3_Max_Score` DOUBLE,
  `BP3_Mean_Linear` DOUBLE,
  `BP3_Max_Linear` DOUBLE,
  `BP3_Top20_Ratio` DOUBLE,
  `BP3_BCell_Level` VARCHAR(32),
  `BP3_Confidence_Score_0_100` DOUBLE,
  `BP3_Confidence_Level` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE mhci_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `mhci_total_rows` INT,
  `mhci_unique_peptides` INT,
  `mhci_hits_percentile_le_1pct` INT,
  `mhci_hits_percentile_le_2pct` INT,
  `mhci_hits_ba_ic50_lt50nM` INT,
  `mhci_hits_ba_ic50_lt500nM` INT,
  `mhci_hits_ba_ic50_lt5000nM` INT,
  `mhci_immunogenicity_positive_count_score_gt0` INT,
  `mhci_immunogenicity_positive_ratio` DOUBLE,
  `mhci_best_peptide` VARCHAR(32),
  `mhci_best_allele` VARCHAR(32),
  `mhci_best_median_percentile` DOUBLE,
  `mhci_best_ba_ic50_nM` DOUBLE,
  `mhci_best_el_percentile` DOUBLE,
  `mhci_best_ba_percentile` DOUBLE,
  `mhci_best_immunogenicity_score` DOUBLE,
  `mhci_binding_level` VARCHAR(32),
  `mhci_immunogenicity_level` VARCHAR(32),
  `mhci_confidence_score_0_100` DOUBLE,
  `mhci_confidence_level` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE mhcii_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `mhcii_total_rows` INT,
  `mhcii_unique_peptides` INT,
  `mhcii_hits_percentile_le_10pct` INT,
  `mhcii_hits_percentile_le_5pct` INT,
  `mhcii_hits_ba_ic50_lt50nM` INT,
  `mhcii_hits_ba_ic50_lt500nM` INT,
  `mhcii_best_peptide` VARCHAR(32),
  `mhcii_best_allele` VARCHAR(32),
  `mhcii_best_median_percentile` DOUBLE,
  `mhcii_best_ba_ic50_nM` DOUBLE,
  `mhcii_best_el_percentile` DOUBLE,
  `mhcii_best_ba_percentile` DOUBLE,
  `mhcii_best_cd4_score` DOUBLE,
  `mhcii_best_combined_score` DOUBLE,
  `mhcii_binding_level` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE immunogenicity_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `overall_final_judgement` VARCHAR(32),
  `overall_confidence_score_0_100` DOUBLE,
  `overall_confidence_level` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE ecoli_expression_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `EcoliExpr_screen_input_id` VARCHAR(64),
  `EcoliExpr_screen_row_index` INT,
  `EcoliExpr_screen_track` VARCHAR(32),
  `EcoliExpr_soluprot_score` DOUBLE,
  `EcoliExpr_netsolp_solubility` DOUBLE,
  `EcoliExpr_netsolp_usability` DOUBLE,
  `EcoliExpr_rp3net_score` DOUBLE,
  `EcoliExpr_ensemble_score` DOUBLE,
  `EcoliExpr_engineering_required` TINYINT(1),
  `EcoliExpr_out_of_domain` TINYINT(1),
  `EcoliExpr_recommended_ecoli_route` VARCHAR(36),
  `EcoliExpr_candidate_bucket` VARCHAR(32),
  `EcoliExpr_screening_note` VARCHAR(224),
  `EcoliExpr_track_reason` VARCHAR(43),
  `EcoliExpr_score_components_available` INT,
  `EcoliExpr_priority_rank_within_track` DOUBLE,
  `EcoliExpr_priority_rank_within_scaffold` INT,
  `EcoliExpr_ensemble_score_weight_sum` DOUBLE,
  `EcoliExpr_netsolp_tool_score` DOUBLE,
  `EcoliExpr_confidence_score_0_100` INT,
  `EcoliExpr_confidence_tier` VARCHAR(32),
  `EcoliExpr_tool_score_spread` DOUBLE,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE oral_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `Oral_ID` VARCHAR(64),
  `Oral_SGF_%` INT,
  `Oral_SIF_%` DOUBLE,
  `Oral_Pepsin_Frags` INT,
  `Oral_Intestine_Frags` INT,
  `Oral_Max_Frag_Len` INT,
  `Oral_DeepDigest_Risk` INT,
  `Oral_Uptake_Rating` VARCHAR(32),
  `Oral_Uptake_Evidence` VARCHAR(64),
  `Oral_CPP_Best_Prob` DOUBLE,
  `Oral_CPP_Best_Seq` VARCHAR(32),
  `Oral_CPP_Positive_Unique_Count` INT,
  `Oral_CPP_HighEff_Unique_Count` INT,
  `Oral_PEPT1_Candidate_Count` INT,
  `Oral_PEPT1_Dipeptide_Count` INT,
  `Oral_PEPT1_Tripeptide_Count` INT,
  `Oral_Rating_Strict` VARCHAR(32),
  `Oral_Rating` VARCHAR(32),
  `Oral_Confidence` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE solubility_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `QCBundleServerOfficial_protein_name` VARCHAR(38),
  `QCBundleServerOfficial_fasta_header` TEXT,
  `QCBundleServerOfficial_scaffold_category` VARCHAR(32),
  `QCBundleServerOfficial_protein_label` VARCHAR(255),
  `QCBundleServerOfficial_sequence_length` INT,
  `QCBundleServerOfficial_pdb_file` VARCHAR(65),
  `QCBundleServerOfficial_pdb_path` VARCHAR(123),
  `QCBundleServerOfficial_assessment_branch` VARCHAR(32),
  `QCBundleServerOfficial_gravy_score` VARCHAR(32),
  `QCBundleServerOfficial_instability_index` DOUBLE,
  `QCBundleServerOfficial_soluprot_score` DOUBLE,
  `QCBundleServerOfficial_camsol_score` DOUBLE,
  `QCBundleServerOfficial_deepsol_score` DOUBLE,
  `QCBundleServerOfficial_a3d_avg_score` VARCHAR(32),
  `QCBundleServerOfficial_a3d_max_score` DOUBLE,
  `QCBundleServerOfficial_final_assessment` VARCHAR(176),
  `QCBundleServerOfficial_protparam_interpretation` VARCHAR(154),
  `QCBundleServerOfficial_soluprot_call` VARCHAR(32),
  `QCBundleServerOfficial_soluprot_interpretation` VARCHAR(98),
  `QCBundleServerOfficial_deepsol_call` VARCHAR(32),
  `QCBundleServerOfficial_deepsol_interpretation` VARCHAR(77),
  `QCBundleServerOfficial_camsol_call` VARCHAR(32),
  `QCBundleServerOfficial_camsol_interpretation` VARCHAR(130),
  `QCBundleServerOfficial_a3d_call` VARCHAR(32),
  `QCBundleServerOfficial_a3d_interpretation` VARCHAR(88),
  `QCBundleServerOfficial_sequence_consensus_code` VARCHAR(64),
  `QCBundleServerOfficial_sequence_consensus` VARCHAR(116),
  `QCBundleServerOfficial_structure_interpretation_code` VARCHAR(64),
  `QCBundleServerOfficial_structure_interpretation` VARCHAR(88),
  `QCBundleServerOfficial_final_assessment_code` VARCHAR(64),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protparam_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `ProtParam_MW_Da` DOUBLE,
  `ProtParam_pI` DOUBLE,
  `ProtParam_GRAVY` DOUBLE,
  `ProtParam_Instability_Index` DOUBLE,
  `ProtParam_Aromaticity` DOUBLE,
  `ProtParam_Ext_Coeff_Reduced` INT,
  `ProtParam_Ext_Coeff_Oxidized` INT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE ted_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `TED_n_residues` DOUBLE,
  `TED_n_high` DOUBLE,
  `TED_n_medium` DOUBLE,
  `TED_n_low` DOUBLE,
  `TED_n_consensus` DOUBLE,
  `TED_high_domains` VARCHAR(49),
  `TED_medium_domains` VARCHAR(86),
  `TED_low_domains` VARCHAR(242),
  `TED_searchable_domain_count` DOUBLE,
  `TED_searchable_high_domain_count` DOUBLE,
  `TED_searchable_medium_domain_count` DOUBLE,
  `TED_domains_with_cath_hit_count` DOUBLE,
  `TED_domains_without_cath_hit_count` DOUBLE,
  `TED_domain_hit_coverage` DOUBLE,
  `TED_mean_novelty_l2_topk` DOUBLE,
  `TED_max_novelty_l2_topk` DOUBLE,
  `TED_best_top1_target` VARCHAR(32),
  `TED_best_top1_cosine` DOUBLE,
  `TED_best_cath_domain_id` VARCHAR(64),
  `TED_best_cath_class_code` DOUBLE,
  `TED_best_cath_class_name` VARCHAR(32),
  `TED_best_cath_architecture_code` DOUBLE,
  `TED_best_cath_architecture_name` VARCHAR(32),
  `TED_best_cath_topology_code` VARCHAR(32),
  `TED_best_cath_topology_name` VARCHAR(80),
  `TED_best_cath_superfamily_code` VARCHAR(32),
  `TED_best_cath_superfamily_name` VARCHAR(80),
  `TED_top_hit_targets_all` VARCHAR(64),
  `TED_top_hit_superfamilies_all` VARCHAR(36),
  `TED_top_hit_topologies_all` VARCHAR(32),
  `TED_best_domain_id` VARCHAR(67),
  `TED_best_domain_confidence` VARCHAR(32),
  `TED_best_domain_ranges` VARCHAR(49),
  `TED_best_domain_novelty_l2_topk` DOUBLE,
  `TED_best_domain_top1_target` VARCHAR(32),
  `TED_best_domain_top1_cosine` DOUBLE,
  `TED_best_domain_cath_superfamily_code` VARCHAR(32),
  `TED_best_domain_cath_superfamily_name` VARCHAR(80),
  `TED_domain_summary_compact` VARCHAR(419),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protrek_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `protrek_confidence` VARCHAR(32),
  `protrek_confidence_reason` VARCHAR(49),
  `protrek_official_seq2seq_band` VARCHAR(32),
  `protrek_official_seq2seq_score` DOUBLE,
  `protrek_official_function_band` VARCHAR(32),
  `protrek_official_function_score` DOUBLE,
  `protrek_structure_support` VARCHAR(32),
  `protrek_structure_support_detail` VARCHAR(41),
  `protrek_seq_top1_hit` VARCHAR(32),
  `protrek_seq_function_top1` TEXT,
  `protrek_str_top1_hit` VARCHAR(32),
  `protrek_str_top1_score` DOUBLE,
  `protrek_str_function_top1` TEXT,
  `protrek_str_function_score` DOUBLE,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE foldseek_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `foldseek_status` VARCHAR(32),
  `foldseek_target` VARCHAR(32),
  `foldseek_evalue` VARCHAR(32),
  `foldseek_bits` DOUBLE,
  `foldseek_alntmscore` DOUBLE,
  `foldseek_qtmscore` DOUBLE,
  `foldseek_ttmscore` DOUBLE,
  `foldseek_lddt` DOUBLE,
  `foldseek_prob` DOUBLE,
  `foldseek_fident` DOUBLE,
  `foldseek_alnlen` DOUBLE,
  `foldseek_qcov` DOUBLE,
  `foldseek_tcov` DOUBLE,
  `foldseek_avg_tm` DOUBLE,
  `foldseek_confidence` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE plmsearch_results (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `PLMSearch_query_original` VARCHAR(760),
  `PLMSearch_query_normalized` VARCHAR(756),
  `PLMSearch_hit_count` INT,
  `PLMSearch_top1_target` VARCHAR(32),
  `PLMSearch_top1_score` DOUBLE,
  `PLMSearch_top1_confidence_level` VARCHAR(32),
  `PLMSearch_top2_target` VARCHAR(32),
  `PLMSearch_top3_target` VARCHAR(32),
  `PLMSearch_top4_target` VARCHAR(32),
  `PLMSearch_top5_target` VARCHAR(32),
  `PLMSearch_top6_target` VARCHAR(32),
  `PLMSearch_top7_target` VARCHAR(32),
  `PLMSearch_top8_target` VARCHAR(32),
  `PLMSearch_top9_target` VARCHAR(32),
  `PLMSearch_top10_target` VARCHAR(32),
  `PLMSearch_top11_target` VARCHAR(32),
  `PLMSearch_top12_target` VARCHAR(32),
  `PLMSearch_top13_target` VARCHAR(32),
  `PLMSearch_top14_target` VARCHAR(32),
  `PLMSearch_top15_target` VARCHAR(32),
  `PLMSearch_top16_target` VARCHAR(32),
  `PLMSearch_top17_target` VARCHAR(32),
  `PLMSearch_top18_target` VARCHAR(32),
  `PLMSearch_top19_target` VARCHAR(32),
  `PLMSearch_top20_target` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE loop_annotations (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `loop_slot` INT,
  `loop_id` VARCHAR(64) NOT NULL,
  `loop_chain_id` VARCHAR(64),
  `loop_start` INT,
  `loop_end` INT,
  `loop_length` INT,
  `loop_flank_pattern` VARCHAR(32),
  `loop_status` VARCHAR(32),
  `loop_label` VARCHAR(35),
  `loop_manual_decision` VARCHAR(32),
  `loop_confidence_tier` VARCHAR(32),
  `loop_candidate_scope` VARCHAR(32),
  `loop_note` VARCHAR(64),
  `loop_is_active` TINYINT(1),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE loop_flexibility_public_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `loop_id` VARCHAR(64) NOT NULL,
  `loop_slot` INT,
  `loop_chain_id` VARCHAR(64),
  `loop_start` DOUBLE,
  `loop_end` DOUBLE,
  `loop_length` DOUBLE,
  `loop_status` VARCHAR(32),
  `loop_manual_decision` VARCHAR(32),
  `loop_confidence_tier` VARCHAR(32),
  `loop_candidate_scope` VARCHAR(32),
  `selected_structure_source` VARCHAR(32),
  `avg_rsa` DOUBLE,
  `avg_asa` TEXT,
  `loop_constraint_score_total_norm` DOUBLE,
  `itsflex_score` VARCHAR(32),
  `itsflex_class` VARCHAR(32),
  `flexscore_ensemble_diversity` DOUBLE,
  `flexibility_score_static` VARCHAR(32),
  `flexibility_score_dynamic` DOUBLE,
  `flexibility_overall_score` DOUBLE,
  `flexibility_evidence_count` INT,
  `flexibility_consensus_label` VARCHAR(32),
  `flexibility_conflict_flag` TINYINT(1),
  `flexibility_result_version` VARCHAR(32),
  `flexibility_result_timestamp` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`) REFERENCES loop_annotations (`protein_row_id`, `structure_unique_sequence_id`, `loop_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE protein_flexibility_summary (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `flexibility_loop_count_total` INT,
  `flexibility_loop_count_primary` INT,
  `flexibility_loop_count_with_dssp` INT,
  `flexibility_loop_count_with_getcontact` INT,
  `flexibility_loop_count_with_itsflex` INT,
  `flexibility_loop_count_with_colabfold` INT,
  `flexibility_avg_rsa_mean` DOUBLE,
  `flexibility_avg_rsa_median` DOUBLE,
  `flexibility_constraint_score_mean` DOUBLE,
  `flexibility_constraint_score_median` DOUBLE,
  `flexibility_itsflex_score_mean` VARCHAR(32),
  `flexibility_itsflex_score_median` VARCHAR(32),
  `flexibility_dynamic_score_mean` DOUBLE,
  `flexibility_dynamic_score_median` DOUBLE,
  `flexibility_consensus_rigid_count` INT,
  `flexibility_consensus_intermediate_count` INT,
  `flexibility_consensus_flexible_count` INT,
  `flexibility_conflict_loop_count` INT,
  `flexibility_summary_version` VARCHAR(64),
  `flexibility_summary_timestamp` VARCHAR(64),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE trusted_source_qc (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `qc_review_mode` VARCHAR(36),
  `trusted_source_group` VARCHAR(32),
  `trusted_source_sampled_for_qc` TINYINT(1),
  `trusted_source_qc_verdict` VARCHAR(32),
  `trusted_source_qc_issue_type` VARCHAR(32),
  `trusted_source_qc_sampling_stratum` VARCHAR(37),
  `trusted_source_qc_sampling_rate` DOUBLE,
  `trusted_source_origin_row_index` INT,
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE manual_review_audit (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `review_timestamp` VARCHAR(32),
  `manual_overall_verdict` VARCHAR(32),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE structure_refinement_audit (
  `protein_row_id` INT NOT NULL,
  `structure_unique_sequence_id` VARCHAR(64) NOT NULL,
  `pre_refine_Final_Tested_Sequence` TEXT,
  `pre_refine_Sequence_Length` INT,
  `pre_refine_Sequence_PDB_ID` VARCHAR(64),
  `pre_refine_local_pdb_exact_structure_paths` VARCHAR(367),
  `pre_refine_local_pdb_fasta_paths` VARCHAR(325),
  `pre_refine_structure_locator` VARCHAR(32),
  `pre_refine_pdb_path` VARCHAR(121),
  `pre_refine_chain_id` VARCHAR(64),
  `pre_refine_chain_count` TINYINT(1),
  `pre_refine_resnum_start` INT,
  `pre_refine_resnum_end` INT,
  `structure_refinement_applied` TINYINT(1),
  `structure_refinement_status` VARCHAR(33),
  `structure_refinement_notes` VARCHAR(118),
  PRIMARY KEY (`protein_row_id`, `structure_unique_sequence_id`),
  FOREIGN KEY (`protein_row_id`, `structure_unique_sequence_id`) REFERENCES active_structure_context (`protein_row_id`, `structure_unique_sequence_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE pipeline_versions (
  `pipeline_name` VARCHAR(64) NOT NULL,
  `source_table` VARCHAR(64) NOT NULL,
  `version` VARCHAR(32),
  `timestamp` VARCHAR(32),
  PRIMARY KEY (`pipeline_name`, `source_table`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE result_metadata (
  `key` VARCHAR(64) NOT NULL,
  `value` VARCHAR(55),
  PRIMARY KEY (`key`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

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
SET FOREIGN_KEY_CHECKS = 1;
