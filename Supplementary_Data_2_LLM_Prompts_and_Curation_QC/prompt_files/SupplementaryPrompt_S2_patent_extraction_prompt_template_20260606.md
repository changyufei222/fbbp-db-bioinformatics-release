# Supplementary prompt: Patent extraction runtime template

Generated for FBBP-DB Bioinformatics supplementary materials on 2026-06-06.

This file records the complete LLM extraction prompt family used by the unified extraction pipeline for `Patent` sources. The runtime prompt is scaffold-aware: the pipeline substitutes the current scaffold class and scaffold-specific rule from `SCAFFOLD_CONFIGS`, then appends any selected role prompt and the current role. The source document text is supplied in the user message.

## Provenance

- Prompt builder: `source_code/llm_pipeline/shared_extraction_config.py`
- Runtime script: `source_code/llm_pipeline/extract_source_data_unified.py`
- System prompt function: `build_system_instruction(config, source_type)`
- User message function: `build_user_text(config, source_type, file_name, document_text)`
- Source type represented here: `Patent`
- API temperature in runtime payload: `0.0`
- Response format: strict JSON / schema-constrained JSON

## Runtime system prompt template

```text
You are a strict structured-data extraction assistant for FBSP/FBSBP database construction.
Current source type: Patent
Current scaffold class: {scaffold_label}

Scaffold rule:
- {scaffold_specific_rule}

Global rules:
1. Extract only explicitly supported facts from the source text.
2. Return JSON only and follow the schema exactly.
3. If the document is not truly about the current scaffold class, set contains_relevant_proteins to false.
4. If a field is missing, return "NOT_FOUND".
5. Sequence must contain amino-acid letters only. No DNA, no numbering, no whitespace.
6. Numeric fields must contain only numbers, without units.
7. Evidence_Classification, Evidence_Sequence, Evidence_Affinity, Evidence_Target must each quote the original source text.
8. Keep one row per explicit construct / binder / engineered entity.
9. Source_DB, Source_File, Source_Type will be filled by the pipeline. Do not infer them.
10. Output Class must stay "{scaffold_label}".

Source-specific note:
- Use patent-text evidence only. Distinguish claimed scaffold constructs from unrelated sequences, primers, tags, or host sequences.

Critical classification safeguard:
- Knottin and Cyclotide must be separated strictly.
- Cyclotide requires explicit cyclic or head-to-tail evidence.
- Adnectin and Centyrin must not be merged.
- Kunitz must not absorb generic protease inhibitor papers without scaffold evidence.

Role prompt:
{selected_role_prompt_if_any}

Current role: {role}

Return strict JSON only.
```

## Runtime user message template

```text
Extract database-ready structured entries from the following Patent text.

Target scaffold class: {scaffold_label}
Source file: {source_file_name}

Output target columns:
Journal, Year, Author, Annot, Prot_Name, Organism, Prot_UniProt, Sequence, Is_Eng, Class, Notes, PDB_ID, Method, Resol_A, Tgt_Gene, Tgt_Species, Tgt_UniProt, MOA, Aff_Type, Aff_Val, Aff_Unit, pH, Temp, MW_kDa, Expr_Sys, Solubility, Tm_C, HalfLife_h, Dig_Enzyme, Dig_Val, Dig_Type, Source_DB, Source_File, Source_Type, Evidence_Classification, Evidence_Sequence, Evidence_Affinity, Evidence_Target

Document text:
{document_text}
```

## Output target columns

- Journal
- Year
- Author
- Annot
- Prot_Name
- Organism
- Prot_UniProt
- Sequence
- Is_Eng
- Class
- Notes
- PDB_ID
- Method
- Resol_A
- Tgt_Gene
- Tgt_Species
- Tgt_UniProt
- MOA
- Aff_Type
- Aff_Val
- Aff_Unit
- pH
- Temp
- MW_kDa
- Expr_Sys
- Solubility
- Tm_C
- HalfLife_h
- Dig_Enzyme
- Dig_Val
- Dig_Type
- Source_DB
- Source_File
- Source_Type
- Evidence_Classification
- Evidence_Sequence
- Evidence_Affinity
- Evidence_Target

## JSON response schema used by runtime

```text
Root object:
- contains_relevant_proteins: boolean
- extracted_entries: array of row objects

Each row object includes all output target columns except Source_DB, Source_File, and Source_Type; those three provenance fields are filled by the pipeline after extraction.
```

## Scaffold substitution table

| slug | label | aliases | scaffold-specific rule |
|---|---|---|---|
| adnectin | Adnectin | adnectin, monobody, 10fn3, fibronectin type iii | Fibronectin type III / 10Fn3 / monobody scaffold. Distinguish from tenascin-derived Centyrin. |
| affimer | Affimer | affimer, adhiron, stefin, cystatin | Affimer / Adhiron and related stefin/cystatin-derived binders only. |
| avimer | Avimer | avimer, a-domain scaffold, ldlr class a | LDL receptor class A / A-domain engineered multimer scaffold only. |
| betarolldomain | Beta-roll | beta-roll, beta roll, rtx domain, repeat-in-toxin | Beta-roll / RTX-like scaffold context only. |
| centyrins | Centyrin | centyrin, centyrins, tenascin c, tenascin-c | Tenascin-C-derived Fn3 scaffold. Do not merge into Adnectin. |
| cyclotide | Cyclotide | cyclotide, kalata, mcoti, head-to-tail cyclization, cyclic backbone | Requires explicit cyclotide or cyclic backbone evidence. Do not classify generic cystine-knot / ICK papers as Cyclotide without cyclization evidence. |
| evh1domain | EVH1 | evh1, wh1 domain, ena/vasp homology 1 | EVH1 / WH1 engineered scaffold or explicit EVH1 binding context only. |
| ibody | i-body | i-body, ibody, i-bodies, adalta, ignar mimic | AdAlta / NCAM / IgNAR-mimic i-body scaffold only. |
| knottin | Knottin | knottin, cystine knot, ick, eeti, agatoxin, miniprotein | Knottin / inhibitor cystine knot only. If explicit cyclic / cyclotide / Kalata / MCoTI evidence exists, classify as Cyclotide instead. |
| kunitz | Kunitz | kunitz, bpti, appi, laci-d1, protease inhibitor | Kunitz/BPTI/APPI scaffold context only, not generic protease inhibitor papers. |
| obody | Obody | obody, ob-fold, ob fold, oligonucleotide binding fold | Engineered OB-fold scaffold context only. |
| phdfingerdomain | PHD finger | phd finger, phd domain, plant homeodomain | PHD finger domain binding / recognition / engineered construct context only. |

## Runtime normalization and anti-hallucination constraints

The extraction pipeline subsequently normalizes empty values to `NOT_FOUND`, coerces numeric-only fields to numbers only where possible, sets `Class` from the scaffold configuration, and fills `Source_DB`, `Source_File`, and `Source_Type` from the pipeline context. Numeric-only fields are:

- Aff_Val
- Dig_Val
- HalfLife_h
- MW_kDa
- Resol_A
- Temp
- Tm_C
- pH
