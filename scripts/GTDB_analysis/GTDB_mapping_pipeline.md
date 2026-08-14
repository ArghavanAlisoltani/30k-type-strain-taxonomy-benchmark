## Recommended workflow: metadata mapping, not GTDB-Tk

GTDB R232 was released in April 2026 and provides `bac120_metadata_r232.tsv.gz` (275.7 MB) and `ar53_metadata_r232.tsv.gz` (6.7 MB). ([GTDB Data][2]) This means the entire mapping job is computationally small.

I would do it as follows.

1. **Create a clean project directory and lightweight Python environment.**

```bash
cd 
mkdir -p gtdb_r232_mapping/{input,data/auxillary_files,scripts,results,logs}
cd gtdb_r232_mapping

# If using mamba
mamba create -n gtdb-r232-map -c conda-forge \
    python=3.12 pandas pyarrow openpyxl

conda activate gtdb-r232-map
```

You do **not** install GTDB-Tk for this step.

Your input should be the exact 30,495 accession.version identifiers used in the paper, for example:

```text
GCA_000005845.2
GCA_000006765.1
GCA_000006925.2
...
```

Put them in:

```text
nano input/type_strain_GCA_30495.txt
```

Confirm:

```bash
wc -l input/type_strain_GCA_30495.txt
sort -u input/type_strain_GCA_30495.txt | wc -l
```

Both should ideally report `30495`.

2. **Download the official GTDB R232 metadata.**

```bash
cd data

BASE="https://data.gtdb.ecogenomic.org/releases/release232/232.0"

wget -c ${BASE}/bac120_metadata_r232.tsv.gz
wget -c ${BASE}/ar53_metadata_r232.tsv.gz
wget -c ${BASE}/MD5SUM.txt

wget -c ${BASE}/auxillary_files/qc_failed_r232.tsv \
    -P auxillary_files/

wget -c ${BASE}/auxillary_files/sp_clusters_r232.tsv \
    -P auxillary_files/

wget -c ${BASE}/auxillary_files/metadata_field_desc.tsv \
    -P auxillary_files/

cd ..
```

The first two metadata files are the key files. The `qc_failed_r232.tsv` file is also important because an accession absent from the normal GTDB metadata may have been evaluated but excluded by GTDB quality control. GTDB explicitly provides this QC-failure file, and `sp_clusters_r232.tsv` records each GTDB species representative, its cluster members, and the ANI radius used for the species cluster. 

3. **Verify the downloads.**

GTDB publishes MD5 checksums for R232, including both metadata files and the auxiliary files. 

```bash
cd data

grep -E \
'bac120_metadata_r232.tsv.gz|ar53_metadata_r232.tsv.gz|qc_failed_r232.tsv|sp_clusters_r232.tsv' \
MD5SUM.txt | md5sum -c -

cd ..
```

You want `OK` for every downloaded file.

4. **Inspect the metadata columns before mapping.**

```bash
python - <<'PY'
import pandas as pd

for f in [
    "data/bac120_metadata_r232.tsv.gz",
    "data/ar53_metadata_r232.tsv.gz"
]:
    x = pd.read_csv(f, sep="\t", compression="gzip", nrows=0)
    print("\n", f)
    for c in x.columns:
        if any(k in c.lower() for k in
               ["accession", "taxonomy", "representative",
                "assembly_level", "type_material"]):
            print(c)
PY
```

The fields we particularly want are:

```text
accession
ncbi_genbank_assembly_accession
gtdb_taxonomy
ncbi_taxonomy
gtdb_genome_representative
gtdb_representative
ncbi_assembly_level
ncbi_type_material_designation
```

The important design point is that you should map your **GCA accession against the NCBI GenBank assembly accession field**, rather than simply matching it to the GTDB `accession` field. GTDB genomes can be represented internally by GenBank- or RefSeq-derived identifiers, so matching through the GenBank assembly accession avoids unnecessarily losing GCA records.

5. **Run the actual mapping.**

Save this as:

```text
scripts/01_map_GCA_to_GTDB_R232.py
```

```python
#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re

ROOT = Path(".")
INPUT = ROOT / "input/type_strain_GCA_30495.txt"

META_FILES = [
    ROOT / "data/bac120_metadata_r232.tsv.gz",
    ROOT / "data/ar53_metadata_r232.tsv.gz",
]

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

# ---------------------------------------------------------
# Read the 30,495 GCA accessions
# ---------------------------------------------------------

query = pd.read_csv(
    INPUT,
    header=None,
    names=["GCA"],
    dtype=str
)

query["GCA"] = query["GCA"].str.strip()

# Extract canonical GCA accession.version
query["GCA"] = query["GCA"].str.extract(
    r"(GCA_\d+\.\d+)",
    expand=False
)

if query["GCA"].isna().any():
    raise ValueError("Some input lines do not contain valid GCA accession.version IDs")

if query["GCA"].duplicated().any():
    print("WARNING: duplicated input GCAs detected")

print("Input genomes:", len(query))
print("Unique GCAs:", query["GCA"].nunique())

# ---------------------------------------------------------
# Read only useful GTDB metadata columns
# ---------------------------------------------------------

wanted = {
    "accession",
    "ncbi_genbank_assembly_accession",
    "gtdb_taxonomy",
    "ncbi_taxonomy",
    "gtdb_genome_representative",
    "gtdb_representative",
    "ncbi_assembly_level",
    "ncbi_type_material_designation",
}

tables = []

for f in META_FILES:

    header = pd.read_csv(
        f,
        sep="\t",
        compression="gzip",
        nrows=0
    )

    available = [x for x in wanted if x in header.columns]

    print(f"\nReading {f}")
    print("Columns retained:", available)

    x = pd.read_csv(
        f,
        sep="\t",
        compression="gzip",
        usecols=available,
        dtype=str,
        low_memory=False
    )

    tables.append(x)

gtdb = pd.concat(tables, ignore_index=True)

print("\nGTDB metadata rows:", len(gtdb))

# ---------------------------------------------------------
# Standardize GenBank accession
# ---------------------------------------------------------

gtdb["GCA"] = gtdb["ncbi_genbank_assembly_accession"].str.extract(
    r"(GCA_\d+\.\d+)",
    expand=False
)

# Check unexpected duplicates
duplicates = gtdb[
    gtdb["GCA"].notna() &
    gtdb["GCA"].duplicated(keep=False)
].copy()

duplicates.to_csv(
    OUT / "GTDB_R232_duplicate_GCA_records.tsv",
    sep="\t",
    index=False
)

print("Duplicate GCA records:", duplicates["GCA"].nunique())

# For ordinary mapping, retain first record.
# If duplicates occur among our 30,495 IDs, inspect them manually.
gtdb_unique = gtdb.drop_duplicates("GCA", keep="first")

# ---------------------------------------------------------
# Join
# ---------------------------------------------------------

mapped = query.merge(
    gtdb_unique,
    on="GCA",
    how="left"
)

mapped["GTDB_R232_status"] = mapped["gtdb_taxonomy"].notna().map(
    {True: "MAPPED", False: "UNMAPPED"}
)

# ---------------------------------------------------------
# Split GTDB taxonomy into standard ranks
# ---------------------------------------------------------

prefix_to_rank = {
    "d__": "GTDB_domain",
    "p__": "GTDB_phylum",
    "c__": "GTDB_class",
    "o__": "GTDB_order",
    "f__": "GTDB_family",
    "g__": "GTDB_genus",
    "s__": "GTDB_species",
}

for col in prefix_to_rank.values():
    mapped[col] = pd.NA

for idx, tax in mapped["gtdb_taxonomy"].items():

    if pd.isna(tax):
        continue

    for item in str(tax).split(";"):

        for prefix, col in prefix_to_rank.items():

            if item.startswith(prefix):
                mapped.at[idx, col] = item[len(prefix):]

# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

mapped.to_csv(
    OUT / "type_strains_30495_GTDB_R232_mapping.tsv",
    sep="\t",
    index=False
)

mapped.to_parquet(
    OUT / "type_strains_30495_GTDB_R232_mapping.parquet",
    index=False
)

mapped[mapped["GTDB_R232_status"] == "MAPPED"].to_csv(
    OUT / "GTDB_R232_mapped.tsv",
    sep="\t",
    index=False
)

mapped[mapped["GTDB_R232_status"] == "UNMAPPED"].to_csv(
    OUT / "GTDB_R232_unmapped.tsv",
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

n = len(mapped)
nm = (mapped["GTDB_R232_status"] == "MAPPED").sum()
nu = n - nm

summary = f"""GTDB R232 mapping summary

Input genomes: {n}
Mapped to GTDB R232: {nm}
Mapped percentage: {100*nm/n:.3f}%
Unmapped: {nu}
Unmapped percentage: {100*nu/n:.3f}%

Bacterial + archaeal GTDB metadata release: R232
"""

print("\n" + summary)

(OUT / "GTDB_R232_mapping_summary.txt").write_text(summary)
```

Run:

```bash
python scripts/01_map_GCA_to_GTDB_R232.py \
    | tee logs/01_map_GTDB_R232.log
```

6. **Check the result immediately.**

```bash
cat results/GTDB_R232_mapping_summary.txt

head results/GTDB_R232_unmapped.tsv
```

The first key revision result will be:

> Of 30,495 NCBI type-material genomes, **X (Y%) were represented in GTDB R232**.

Do **not** assume that all 30,495 will map. R232's files were produced in mid-April 2026, whereas your NCBI dataset was frozen on May 4, 2026, so some assemblies can legitimately postdate the GTDB snapshot. ([GTDB Data][2])

7. **Investigate every unmapped accession.**

For each unmapped GCA, classify it into:

```text
GTDB QC failed
not represented in R232
possible accession/version mismatch
other mapping issue
```

Start with `qc_failed_r232.tsv`. GTDB provides this specifically as the list of assemblies that failed its internal QC criteria. 

Do **not** silently discard the unmapped genomes. Their number and reason for absence should appear in your Methods/Supplementary table.

8. **Join GTDB onto your existing NCBI taxonomy—not the other way around.**

Your final file should look approximately like:

```text
GCA

NCBI_domain
NCBI_kingdom
NCBI_phylum
NCBI_class
NCBI_order
NCBI_family
NCBI_genus
NCBI_species

GTDB_domain
GTDB_phylum
GTDB_class
GTDB_order
GTDB_family
GTDB_genus
GTDB_species

GTDB_R232_status
gtdb_genome_representative
gtdb_representative
```

**Use your May 4 NCBI taxonomy columns as the NCBI side of the comparison.** Do not replace them with the `ncbi_taxonomy` field bundled inside GTDB metadata, because the whole purpose is to compare your frozen NCBI reference against GTDB R232.

Also note that standard GTDB taxonomy uses the major ranks domain → phylum → class → order → family → genus → species; you should therefore make the main NCBI/GTDB comparison at those common ranks rather than trying to force your NCBI kingdom categories into GTDB. GTDB's reference taxonomy and bacterial tree are explicitly organized using these phylogenomic ranks. 

## Time and infrastructure

For the **recommended metadata mapping**, this is easy:

**Compute:** 1–4 CPU cores are sufficient.
**RAM:** 8 GB is probably enough; 16 GB is comfortable.
**Disk:** I would reserve ~5–10 GB for metadata, decompression/temporary space, outputs and later analyses.
**GPU:** none.
**Genome FASTA files:** not required.

The two main compressed metadata files total only about **282 MB**, although the bacterial file represents nearly 879,000 bacterial genomes. ([GTDB Data][2]) I would budget roughly **10–30 minutes end-to-end**, mostly depending on download speed; once downloaded, the actual 30,495-accession join should take only minutes on an ordinary modern computer. This is an estimate based on the released file sizes, not a GTDB benchmark.

### What if an accession is not already in GTDB?

Only then would I consider **GTDB-Tk**, and initially only for the unmatched genomes.

That is a very different computational problem. Current GTDB-Tk documentation lists approximately **100 GB RAM for Archaea, 140 GB RAM for Bacteria, ~100 GB of reference-data storage, and about 90 minutes per 1,000 genomes using 64 CPUs**; R232 requires GTDB-Tk ≥2.7.0. ([EcoGenomics][3]) GTDB currently documents Bioconda installation as:

```bash
mamba create -n gtdbtk-2.7.2 \
    -c conda-forge -c bioconda \
    gtdbtk=2.7.2

conda activate gtdbtk-2.7.2

download-db.sh
```

and `GTDBTK_DATA_PATH` must point to the unarchived R232 reference package. ([EcoGenomics][3])

If you unnecessarily classified **all 30,495 genomes** from sequence, the published rate extrapolates to roughly **46 hours at 64 CPUs**, before allowing for I/O and workflow overhead; realistically I would budget about **2–3 days on a high-memory HPC node**. ([EcoGenomics][3])

**But I would not do that now.** First perform the metadata mapping. I expect that to resolve the overwhelming majority of the 30,495 genomes at negligible computational cost. Then we can look at the exact unmapped set and decide whether GTDB-Tk is necessary for any of them.


[1]: https://gtdb.ecogenomic.org/?utm_source=chatgpt.com "Genome Taxonomy Database: GTDB"
[2]: https://data.gtdb.ecogenomic.org/releases/release232/232.0/ "GTDB Data - /releases/release232/232.0/"
[3]: https://ecogenomics.github.io/GTDBTk/installing/index.html "Installing GTDB-Tk — GTDB-Tk 2.7.2 documentation"

