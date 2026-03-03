## Function

```python
hccnstable(
    self,
    tool_hap1_cna_files: List[str],
    tool_hap2_cna_files: List[str],
    tool_names: List[str],
    changes_file: str,
    tree_file: str,
    outfile: str = "evolution_cn_stability_acc.csv",
    profile_bin_size=100000
) -> pd.DataFrame
```

This function evaluates **copy-number stability accuracy** along a given phylogenetic tree.

Given:

- a **tree** (Newick format),
- a **segment change table** (`changes_file`),
- and each tool’s haplotype-specific CNA profiles (hap1/hap2),

it builds a set of “stability checks” derived from the tree topology and then measures, for each tool, how often predicted CN states satisfy the expected stability constraints. Results are reported as **ACC** (mean of boolean check outcomes), stratified by change `Type`.

------

## Parameters

| Name                  | Type        | Description                                                  |
| --------------------- | ----------- | ------------------------------------------------------------ |
| `tool_hap1_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 1. Must align with `tool_names`. |
| `tool_hap2_cna_files` | `List[str]` | List of tool CNA profile CSVs for haplotype 2. Must align with `tool_names`. |
| `tool_names`          | `List[str]` | Tool names used in result rows.                              |
| `changes_file`        | `str`       | Path to the segment/change table used for stability checks (must include `Segment`, `Type`, and other fields used by internal helpers). |
| `tree_file`           | `str`       | Path to the phylogenetic tree file in **Newick** format.     |
| `outfile`             | `str`       | Output filename for the summary table. Default: `"evolution_cn_stability_acc.csv"`. |
| `profile_bin_size`    | `int`       | Bin size used to split/standardize both segments and predictions. Default: `100000` (100kb). |

------

## Input File Format

### `tree_file` (Newick tree)

Read using Biopython:

```python
tree = Phylo.read(tree_file, "newick")
```

So the tree must be a valid Newick string/file (tip names must match whatever identifiers your stability logic expects).

### `changes_file` (segment/change table)

Loaded via:

```python
change_df = read_and_drop_empty(changes_file)
```

Expected to contain at least:

- `Segment`: genomic interval key used for binning and joining
- `Type`: category label used for stratified reporting

Optional:

- `Haplotype`: if present, will be normalized via `self._normalize_hap_label(...)`

Preprocessing steps:

1. Add a checklist derived from the tree:

```python
change_df = self._add_check_list(tree, change_df)
```

1. Re-bin segments:

```python
change_df = split_all_regions(change_df.set_index("Segment"), profile_bin_size)
change_df = change_df.reset_index().rename(columns={"index": "Segment"})
```

1. Normalize haplotype labels if provided.

------

## Tool CNA Profile Format

Each file in `tool_hap1_cna_files` / `tool_hap2_cna_files` is expected to be a CNA matrix with:

- `region` column
- remaining columns as cell IDs
- values: haplotype-specific CN states

They are loaded and re-binned:

```python
p_h1 = split_all_regions(p_h1.set_index("region"), profile_bin_size).reset_index()
p_h2 = split_all_regions(p_h2.set_index("region"), profile_bin_size).reset_index()
```

------

## Evaluation Logic

For each tool:

1. Copy the processed change table:

```python
proc = change_df.copy()
```

1. Populate row-wise stability check outcomes using predictions:

```python
self._process_change_rows(proc, p_h1, p_h2)
```

This helper is expected to add/overwrite a column:

- `result`: a per-row boolean or numeric indicator (e.g., 1.0/0.0) showing whether the stability condition is satisfied for that row.

1. Aggregate accuracy per `Type`:

```python
acc = proc[proc["Type"] == t]["result"].mean()
```

1. Append one summary row per `(Tool × Type)`:

- `Tool`
- `Type`
- `ACC` (mean stability satisfaction rate)

------

## Output

A single CSV summary is written to `self.output_dir`:

```
os.path.join(self.output_dir, outfile)
```

Default:

```
{self.output_dir}/evolution_cn_stability_acc.csv
```

------

## Output Table Schema

One row per **(Tool × Type)**:

| Column | Meaning                                                  |
| ------ | -------------------------------------------------------- |
| `Tool` | Tool name                                                |
| `Type` | Change/stability category from `changes_file`            |
| `ACC`  | Mean of `result` for that type (stability accuracy rate) |

Notes:

- If `result` is boolean, pandas will treat `True/False` as `1/0` when computing `.mean()`.
- If a `Type` has no rows, it will not appear in the output (because types are taken from `proc['Type'].unique()`).

------

## Return Value

Returns a `pd.DataFrame` with columns:

- `Tool`
- `Type`
- `ACC`

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.hccnstable(
    tool_hap1_cna_files=[
        "/path/to/chisel/hap1.csv",
        "/path/to/signals/hap1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/chisel/hap2.csv",
        "/path/to/signals/hap2.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    changes_file="/path/to/gt/changes_stability.csv",
    tree_file="/path/to/gt/tree.newick",
    profile_bin_size=100000,
    outfile="evolution_cn_stability_acc.csv",
)

print(df)
```

