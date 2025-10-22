# HCBench (Haplotype- and Clone-specific copy number Benchmarking)

### 🧬 hcbench Parser — CNA Caller Benchmark Parser

**hcbench Parser** is a tool designed to standardize the output files from CNA callers.
Currently, it supports the following tools:

* **CHISEL**
* **SEACON**
* **SIGNALS**

These parsers unify various caller output formats into a consistent **wide-format matrix**, where each row represents a genomic region and each column represents a single cell.

---

### 📥 Input File Requirements

#### 1. CHISEL

* File format: TSV (`calls.tsv`)

* Required columns:

  ```
  #CHR, START, END, CELL, CORRECTED_HAP_CN
  ```

* The `START` value is incremented by 1 to follow the standard 1-based genomic coordinate system.

#### 2. SEACON

* File format: TSV (`calls.tsv`)

* Required columns:

  ```
  chrom, start, end, cell, CN
  ```

* The `CN` column may contain comma-separated copy numbers, which will be replaced with vertical bars (`|`).

#### 3. SIGNALS

* File format: TSV (`hscn_data.tsv`)

* Required columns:

  ```
  chr, start, end, cell_id, state_AS_phased
  ```

---

### 📤 Output Format

The output file is a CSV:

* Row index: `region`, e.g. `chr1:1000001-1500000`
* Each column corresponds to a single-cell ID.
* Each cell value represents the corresponding copy number or state value.

| region                 | clone1_cell1 | clone1_cell10 | clone1_cell2 | clone1_cell3 |
| ---------------------- | ------------ | ------------- | ------------ | ------------ |
| chr1:5000001-10000000  | 1\|2         | 1\|2          | 1\|2         | 1\|2         |
| chr1:15000001-20000000 | 1\|2         | 1\|2          | 1\|2         | 1\|2         |
| chr1:20000001-25000000 | 1\|2         | 1\|2          | 1\|2         | 1\|2         |

---

### 🚀 Usage Example

```python
from hcbench.parsers.chisel import ChiselParser
from hcbench.parsers.seacon import SeaconParser
from hcbench.parsers.signals import SignalsParser

# Parse CHISEL output
chisel_input = "./demo_output/chisel/calls/calls.tsv"
chisel_output = "./output/chisel_haplotype_combined.csv"
chisel_parser = ChiselParser(chisel_input, chisel_output)
chisel_parser.run()

# Parse SEACON output
seacon_input = "./demo_output/seacon/calls.tsv"
seacon_output = "./output/seacon_haplotype_combined.csv"
seacon_parser = SeaconParser(seacon_input, seacon_output)
seacon_parser.run()

# Parse SIGNALS output
signals_input = "./demo_output/signals/output/hscn_data.tsv"
signals_output = "./output/signals_haplotype_combined.csv"
signals_parser = SignalsParser(signals_input, signals_output)
signals_parser.run()
```
---

### 🧩 Dependencies

```
pandas>=2.0.0
```

---


- gt-bench
  + ....
  + 
- real-bench
  + ...
