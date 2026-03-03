## Function

```python
dolazactree(
    self,
    tool_cna_files: List[str],
    tool_names: List[str],
    outfile: Optional[str] = "parsimony_score.csv"
) -> pd.DataFrame
```

This function evaluates **phylogenetic parsimony scores** for inferred CNA profiles using the external tool **LAZAC**.

For each tool:

1. Converts the CNA matrix into LAZAC-compatible input format.
2. Runs LAZAC with nearest-neighbor interchange (NNI) tree search.
3. Extracts the final parsimony score from the LAZAC output.
4. Collects scores across tools into a summary table.

The function writes a single CSV file containing the parsimony scores of all tools.

------

## Parameters

| Name             | Type            | Description                                                  |
| ---------------- | --------------- | ------------------------------------------------------------ |
| `tool_cna_files` | `List[str]`     | List of CNA profile CSV files from different tools. Must align with `tool_names` in order. |
| `tool_names`     | `List[str]`     | Tool names used for labeling outputs and result table.       |
| `outfile`        | `Optional[str]` | Output filename for the summary table. Default: `"parsimony_score.csv"`. |

------

## Input File Format

Each file in `tool_cna_files` is expected to be a CNA matrix:

- rows: genomic regions
- columns: cell IDs
- values: CNA states (e.g., `"1|1"`, `"2|1"`)

Example:

```csv
region,cell_001,cell_002,cell_003
chr1:1-100000,1|1,1|1,2|1
chr1:100001-200000,1|1,1|1,2|1
```

Before running LAZAC:

- The CNA matrix is converted into a LAZAC-compatible input file via `create_lazac_input`.
- The converted file is saved as:

```
{output_dir}/lazac_{tool_name}/{tool_name}_cn_profile.csv
```

------

## External Dependency

This function requires:

- `lazac` executable available in the system PATH.
- The command executed is:

```bash
lazac nni <cn_profile_file> -a 2 -o <output_prefix>
```

Where:

- `nni` performs nearest-neighbor interchange search.
- `-a 2` specifies diploid assumption.
- `-o` sets output prefix.

------

## Output

For each tool, a dedicated directory is created:

```
{self.output_dir}/lazac_{tool_name}/
```

LAZAC generates output files including:

- `{tool_name}_info.json` (contains final parsimony score)
- additional tree-related outputs

### Summary Table (`outfile`)

Saved to:

```
os.path.join(self.output_dir, outfile)
```

Structure:

| Column      | Meaning                                           |
| ----------- | ------------------------------------------------- |
| `Tool`      | Tool name                                         |
| `Parsimony` | Final parsimony score extracted from LAZAC output |

Example:

```csv
Tool,Parsimony
CHISEL,14331085
Alleloscope,34
SIGNALS,5577962
```

Lower parsimony scores indicate simpler evolutionary histories under the parsimony model.

------

## Return Value

Returns a `pd.DataFrame` containing:

| Column      | Meaning         |
| ----------- | --------------- |
| `Tool`      | Tool name       |
| `Parsimony` | Parsimony score |

------

## Example

```python
from hcbench.realbench.realbench import RealBench

bench = RealBench(output_dir="out/real_output")

bench.dolazactree(
    tool_cna_files=[
        "/path/to/chisel/cna.csv",
        "/path/to/signals/cna.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    outfile="parsimony_score.csv"
)
```

------

