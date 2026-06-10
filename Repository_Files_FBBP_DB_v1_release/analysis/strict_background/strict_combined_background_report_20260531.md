# Strict combined background analysis

Strict combined background analysis used only experimental PDB/X-ray structures with resolution <= 3.0 A, collapsed loop sequences into greedy same-length sequence clusters at >= 90% identity, removed exact FBBP-loop sequences from the AIR background, and then estimated the AIR expectation by matching the AIR clusters to the FBBP cluster length distribution.

After this combined restriction, the FBBP strict subset retained 177 records and 150 sequence clusters; the AIR background retained 20537 records and 20375 clusters. The FBBP cluster-level ITsFlexible-label flexible fraction was 0.820, compared with a length-matched strict AIR expectation of 0.245 (bootstrap 95% interval 0.180-0.310; one-sided resampling P=0.0010). The length-stratified Mantel-Haenszel OR was 13.52, and the length-adjusted logistic OR was 17.48 (11.44-26.71).

## Machine-readable summary

```json
{
  "analysis": "strict_combined_background_pdb_resolution_dedup_length_matched",
  "resolution_threshold_angstrom": 3.0,
  "sequence_identity_dedup_threshold": 0.9,
  "exact_fbbp_air_sequence_overlap_removed_from_air": 1,
  "n_fbbp_records_after_pdb_xray_resolution_filter": 177,
  "n_fbbp_sequence_clusters_90id": 150,
  "n_air_records_after_pdb_xray_resolution_filter": 20537,
  "n_air_sequence_clusters_90id": 20375,
  "fbbp_cluster_mixed_label_count": 0,
  "air_cluster_mixed_label_count": 33,
  "n_fbbp_clusters_length_matched": 150,
  "n_fbbp_clusters_missing_air_length": 0,
  "observed_fbbp_flexible_fraction": 0.82,
  "strict_air_length_matched_expected_fraction": 0.2449288888888889,
  "observed_minus_expected": 0.5750711111111111,
  "bootstrap_ci_lower": 0.18,
  "bootstrap_ci_upper": 0.31008333333333327,
  "bootstrap_p_one_sided_air_ge_fbbp": 0.000999000999000999,
  "n_bootstrap": 1000,
  "mh_strata_used": 14,
  "mantel_haenszel_or": 13.517253494734332,
  "logit_n_obs": 20525,
  "length_adjusted_dataset_or": 17.479918918647606,
  "length_adjusted_or_ci_lower": 11.438576865099913,
  "length_adjusted_or_ci_upper": 26.71202624294518,
  "length_adjusted_p_value": 6.402300207041636e-40,
  "logit_converged": 1
}
```