## Function

```Python
unique_cna_profile_counts(
    self,
    tool_cna_files: List[str],
    tool_names: List[str],
    outfile: str = "unique_cna_profile_counts.csv"
) -> None
```

This function counts the **number of unique cells** contained in each tool’s CNA result file and exports a summary table to `self.output_dir`.

## Parameters

| **Name**         | **Type**    | **Description**                                              |
| ---------------- | ----------- | ------------------------------------------------------------ |
| `tool_cna_files` | `List[str]` | List of CNA result file paths produced by different tools. Must be aligned with `tool_names` by order. |
| `tool_names`     | `List[str]` | List of tool names. Used as the row index in the output table. Must have the same length as `tool_cna_files`. |
| `outfile`        | `str`       | The name of the output CSV file. Defaults to `"unique_cna_profile_counts.csv"`. |

## Output

The function writes a single CSV file:

- **Path**: `os.path.join(self.output_dir, outfile)`
- **Columns**:
  - `unique_cells_profile`: the number of unique cells for each tool
- **Index**:
  - `Tool`: tool name (from `tool_names`)

Example output (`unique_cna_profile_counts.csv`):

```
Tool,unique_cells_profile
CHISEL,1099
Alleloscope,1050
SIGNALS,980
```

## Example

```text
from hcbench.gtbench import gtbench

runner = gtbench.GTBench(output_dir="out/gt_output/", bin_size=100000 , genome="hg38")

tool_cna_files = [
    "out/chisel/haplotype_combined.csv",
    "out/signals/haplotype_combined.csv",
    "out/alleloscope/haplotype_combined.csv",
]
tool_names = ["CHISEL", "SIGNALS", "Alleloscope"]

runner.unique_cna_profile_counts(
    tool_cna_files=tool_cna_files,
    tool_names=tool_names,
    outfile="unique_cna_profile_counts.csv"
)
```

