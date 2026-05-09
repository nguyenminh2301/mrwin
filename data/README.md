# Data

This directory is intentionally empty. The simulation engine generates all data synthetically; no external data files are needed to run the replication or publication modes.

## UK Biobank data access (for companion paper P3)

The illustrative application in Section 7 of the manuscript uses UK Biobank data (Application Number: *[to be filled]*). UK Biobank data are available to bona fide researchers through a formal application process:

1. **Application:** https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access
2. **Data Fields used:**
   - Field 22006: White-British ancestry indicator
   - Field 22020: Used in genetic principal components
   - SBP PRS: 535 SNPs from Evangelou et al. (2018), *Nat Genet* 50:1412–1425
3. **Phenotype definitions:**
   - Priority 1: All-cause mortality (Data Field 40000, 40007)
   - Priority 2: Kidney replacement therapy (ICD-10 codes Z94.0, Z99.2 in HES; OPCS-4 codes M01–M02)
   - Priority 3: Sustained >40% eGFR decline (CKD-EPI from serum creatinine, Data Field 30700)
4. **GWAS summary statistics:** Evangelou et al. (2018) are publicly available from the GWAS Catalog (Study Accession: GCST006624)

## LD clumping parameters

```
--clump-p1 5e-8
--clump-r2 0.001
--clump-kb 10000
--clump-field P
```

## Checksums

The `REPLICATION_CHECKSUMS.txt` file in the repository root contains MD5 checksums for all JSON output files from the audited replication run. Use these to verify byte-identical reproduction.
