# Strict Residue-Level Binding Report

This strict result only includes loops with explicit `target_chain_id`, existing experimental multichain structures, and successful residue-level contact extraction.

- Eligible loops with residue-level labels: `235`
- Eligible structures: `100`
- Binding residues: `316`
- Non-binding residues: `1516`
- Residues with DSSP rel_asa matched: `1832`
- ΔASA-supported residues (`ΔASA_proxy > 1.0 Å²`): `294`
- Core-interface residues (`ΔASA_proxy >= 5.0 Å²`): `220`
- Core-interface strict residues (`ΔASA_proxy >= 10.0 Å²`): `199`

## Output Files
- `<local_path_removed>/binding_residue_labels_strict_v1.csv`
- `<local_path_removed>/binding_loop_residue_summary_strict_v1.csv`
- `<local_path_removed>/binding_vs_nonbinding_residue_stats_strict_v1.csv`

## Overall Strict Statistics
- `distance_defined_label / rel_asa`: binding_mean=`0.3475264915588282`, nonbinding_mean=`0.44269720674727725`, delta_mean=`-0.09517071518844905`, p=`1.4593154702790064e-07`
- `distance_defined_label / delta_asa_proxy`: binding_mean=`27.559620281072597`, nonbinding_mean=`0.06244653134487734`, delta_mean=`27.497173749727718`, p=`0.0`
- `delta_asa_supported_label / rel_asa`: binding_mean=`0.34822329519531486`, nonbinding_mean=`0.44120265798051983`, delta_mean=`-0.09297936278520497`, p=`7.883599165334942e-07`
- `core_interface_label / rel_asa`: binding_mean=`0.3534033934943052`, nonbinding_mean=`0.43622741327091497`, delta_mean=`-0.08282401977660975`, p=`0.00014856727148126675`
- `core_interface_strict_label / rel_asa`: binding_mean=`0.348273968410643`, nonbinding_mean=`0.435787395620174`, delta_mean=`-0.08751342720953098`, p=`0.00012420676309637688`