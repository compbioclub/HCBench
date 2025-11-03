

# HCBench

##  🧬 HCBench Parsers 

Welcome to the **HCBench Parser Suite** — a unified framework for standardizing single-cell copy number alteration (CNA) outputs from multiple tools.

Each parser converts heterogeneous output formats from different algorithms (CHISEL, Alleloscope, CNRein, SEACON, and SIGNALS) into a single consistent format suitable for downstream benchmarking and visualization.

------

### 📖 Overview

| Parser                           | Input Type                      | Input Format       | Output                                  |
| -------------------------------- | ------------------------------- | ------------------ | --------------------------------------- |
| [CHISEL](parsers/chisel.md)      | `calls.tsv` + `mapping.tsv`     | Tab-delimited text | Standardized CNA matrix + clone mapping |
| [Alleloscope](parsers/seacon.md) | `.rds` files                    | R serialized data  | Standardized CNA matrix + clusters      |
| [CNRein](parsers/cnrein.md)                       | `CNReinPrediction.csv`          | CSV                | Standardized CNA matrix                 |
| [SEACON](parsers/seacon.md)                       | `calls.tsv`                     | Tab-delimited text | Standardized CNA matrix                 |
| [SIGNALS](parsers/signals.md)                      | `hscn.rds` (exported to `.tsv`) | Tab-delimited text | Standardized CNA matrix                 |

All parsers generate outputs in the same canonical structure, making it possible to directly compare results across tools in the **HCBench** benchmarking pipeline.

------

### ⚙️ Unified Output Format

After parsing, every tool produces a standardized set of files in its designated output folder:

```
haplotype_combined.csv
haplotype_1.csv
haplotype_2.csv
minor.csv
major.csv
minor_major.csv
clusters.csv    # optional (for CHISEL / Alleloscope)
```

#### Main File

- **`haplotype_combined.csv`** — the core CNA matrix (regions × cells)
   Each entry represents the haplotype-specific copy number as `"hap1|hap2"`.

#### Optional Files

- `haplotype_1.csv`, `haplotype_2.csv` — per-haplotype copy number matrices.
- `minor.csv` stores copy number values of the *minor* allele for each genomic segment and cell.
- `major.csv` stores copy number values of the *major* allele.
- `minor_major.csv` contains both numeric columns (`minor`, `major`) side by side for convenience in downstream analysis.

------

### 🚀 Quick Start Example

```python
from hcbench.parsers.chisel import ChiselParser
from hcbench.parsers.alleloscope import AlleloscopeParser
from hcbench.parsers.cnrein import CNReinParser
from hcbench.parsers.seacon import SeaconParser
from hcbench.parsers.signals import SignalsParser
```

#### Example: Parsing CHISEL Output

```python
chisel_input = "/demo_output/chisel/calls/calls.tsv"
chisel_output = "/output/chisel/"

parser = ChiselParser(chisel_input, chisel_output)
parser.run()
```

#### Example: Parsing Alleloscope Output

```python
alleloscope_parser = AlleloscopeParser(
    genotypes_rds_path="/demo_output/alleloscope/genotypes_all_Sample.rds",
    seg_table_rds_path="/demo_output/alleloscope/seg_table_all_Sample.rds",
    output_path="/output/alleloscope/"
)
alleloscope_parser.run()
```

------

### 🌐 Documentation Index

- [CHISEL Parser](parsers/chisel.md)
- [Alleloscope Parser](parsers/alleloscope.md)
- [CNRein Parser](parsers/cnrein.md)
- [SEACON Parser](parsers/seacon.md)
- [SIGNALS Parser](parsers/signals.md)

