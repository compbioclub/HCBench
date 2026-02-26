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