## Function

```python
clusterConsistency(
    tool_clone_files: List[str],
    tool_names: List[str],
) -> pd.DataFrame
```

For each tool, the function loads  `clusters.csv`file and calculates:

- **ARI**: Adjusted Rand Index between `cell_id` and `clone_id`
- **AMI**: Adjusted Mutual Information between `cell_id` and `clone_id`

It aggregates results across tools into a single table, writes it to:

- `{self.output_dir}/clustering_result.csv`

and returns the result as a DataFrame.

------

## Parameters

| Name               | Type        | Description                                                  |
| ------------------ | ----------- | ------------------------------------------------------------ |
| `tool_clone_files` | `List[str]` | List of file paths, one per tool. Each file should be a CSV containing at least `cell_id` and `clone_id` columns. |
| `tool_names`       | `List[str]` | List of tool names aligned with `tool_clone_files` (same length and order). |

------

## Input File Format

Each `tool_clone_files[i]` is expected to be a CSV with at least:

- `cell_id`: string cell identifier, expected format like `"<cluster>_<rest>"` so that `cell_id.split("_")[0]` yields the reference cluster label
- `clone_id`: clone assignment label from the tool

Example:

```csv
cell_id,clone_id
A_cell0001,1
A_cell0002,1
B_cell0003,2
```

From this, the function derives:

- `clusters1 = ["A", "A", "B"]`
- `clusters2 = ["1", "1", "2"]`

------

## Return Type

- `pd.DataFrame`

------

## Returns

A DataFrame with one row per tool:

| Column | Meaning                                                      |
| ------ | ------------------------------------------------------------ |
| `Tool` | Tool name from `tool_names`                                  |
| `ARI`  | Adjusted Rand Index between derived clusters and `clone_id`  |
| `AMI`  | Adjusted Mutual Information between derived clusters and `clone_id` |

------

## Output

- Writes a CSV summary to:
  `os.path.join(self.output_dir, "clustering_result.csv")`

------

## Example

```python
from hcbench.gtbench import gtbench

# Suppose self.output_dir = "out/gt_output"

gtbench_runner = gtbench.GTBench(
    output_dir=f"out/gt_output/")

tool_clone_files = [
    "out/chisel/clusters.csv",
    "out/alleloscope/clusters.csv",
    "out/signals/clusters.csv",
]
tool_names = ["CHISEL", "Alleloscope", "SIGNALS"]

# Call from within your class instance that has self.output_dir (out/gt_output/) defined
df = gtbench_runner.clusterConsistency(
    tool_clone_files=tool_clone_files,
    tool_names=tool_names,
)

print(df)
# Expected output (values are illustrative):
#           Tool     ARI     AMI
# 0       CHISEL  0.230   0.301
# 1  Alleloscope  0.5321  0.4880
# 2      SIGNALS  0.8123  0.7451
```

After running, you will also find:

```
out/gt_output/clustering_result.csv
```

containing the same summary table.

