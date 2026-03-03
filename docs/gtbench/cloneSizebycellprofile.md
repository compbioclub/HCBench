## Function

```python
cloneSizebycellprofile(
    self,
    gt_cna_file: str,
    tool_cna_files: List[str],
    tool_names: List[str],
    outfile: str = "clone_size_by_cell_profile.csv",
) -> pd.DataFrame
```

This function evaluates clone-size consistency at the **cell profile level**.
It first derives clone sizes from the GT CNA matrix (cells sharing identical bin-level CNA profiles are grouped as one profile-clone), then compares each tool's inferred profile-clone sizes against GT.

The function writes:

- a detailed per-cell comparison table
- a mean summary grouped by GT `cluster_size`

------

## Parameters

| Name             | Type        | Description                                                  |
| ---------------- | ----------- | ------------------------------------------------------------ |
| `gt_cna_file`    | `str`       | Path to GT CNA profile CSV.                                  |
| `tool_cna_files` | `List[str]` | List of tool CNA profile CSV files. Must align with `tool_names` in order. |
| `tool_names`     | `List[str]` | Tool names used as output column prefixes (`{Tool}_pred_size`). |
| `outfile`        | `str`       | Detailed output filename. Default: `"clone_size_by_cell_profile.csv"`. |

------

## Input File Format

`gt_cna_file` and each file in `tool_cna_files` are expected to be CNA matrices:

- rows: genomic regions (first column as index when reading)
- columns: cell IDs
- values: CNA states (e.g., `"1|1"`)

Example:

```csv
region,cell_001,cell_002,cell_003
chr1:1-100000,1|1,1|1,2|1
chr1:100001-200000,1|1,1|1,2|1
```

Implementation details from source:

- missing values are filled with `"1|1"`
- for tool matrices, all-empty columns are dropped before processing

------

## Output

Two files are written to `self.output_dir`:

1. `os.path.join(self.output_dir, outfile)`
2. `os.path.join(self.output_dir, f"mean_{outfile}")`

### 1) Detailed Table (`outfile`)

Index: GT cell ID

| Column             | Meaning                                                      |
| ------------------ | ------------------------------------------------------------ |
| `cluster_size`     | GT profile-clone size for that cell                          |
| `{Tool}_pred_size` | Predicted profile-clone size for the same cell ID from each tool (NaN if cell is missing in tool result) |

### 2) Mean Table (`mean_{outfile}`)

Grouped by GT `cluster_size`, numeric prediction columns are averaged:

| Column             | Meaning                                                      |
| ------------------ | ------------------------------------------------------------ |
| `cluster_size`     | GT profile-clone size                                        |
| `{Tool}_pred_size` | Mean predicted size for GT clones of that size               |

------

## Return Type

The function signature annotates `-> pd.DataFrame`, but current implementation does not explicitly return a value (runtime return is `None`).

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

bench.cloneSizebycellprofile(
    gt_cna_file="/path/to/gt/haplotype_combined.csv",
    tool_cna_files=[
        "/path/to/signals/haplotype_combined.csv",
        "/path/to/seacon/haplotype_combined.csv",
    ],
    tool_names=["signals", "seacon"],
    outfile="clone_size_by_cell_profile.csv",
)
```
