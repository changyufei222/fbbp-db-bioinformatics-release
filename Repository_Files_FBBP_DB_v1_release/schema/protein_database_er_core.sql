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

-- Indexes
CREATE INDEX idx_protein_identifiers_protein_id ON protein_identifiers (`protein_id`);
CREATE INDEX idx_domains_protein_id ON domains (`protein_id`);
CREATE INDEX idx_domains_source_id ON domains (`source_id`);
CREATE INDEX idx_target_species_variants_target_concept_id ON target_species_variants (`target_concept_id`);
CREATE INDEX idx_structural_source_info_domain_id ON structural_source_info (`domain_id`);
CREATE INDEX idx_active_structure_context_protein_row_id ON active_structure_context (`protein_row_id`);
CREATE INDEX idx_interactions_domain_id ON interactions (`domain_id`);
CREATE INDEX idx_interactions_target_variant_id ON interactions (`target_variant_id`);
CREATE INDEX idx_interactions_source_id ON interactions (`source_id`);
CREATE INDEX idx_interactions_protein_row_id_structure_unique_sequence_id ON interactions (`protein_row_id`, `structure_unique_sequence_id`);
CREATE INDEX idx_affinity_data_interaction_id ON affinity_data (`interaction_id`);
CREATE INDEX idx_digestive_assays_domain_id ON digestive_assays (`domain_id`);
CREATE INDEX idx_functional_annotations_protein_id ON functional_annotations (`protein_id`);
