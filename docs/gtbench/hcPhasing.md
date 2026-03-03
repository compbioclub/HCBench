## Function

```python
hcPhasing(
    self,
    tool_hap1_cna_files: List[str],
    tool_hap2_cna_files: List[str],
    tool_names: List[str],
    ground_truth_hap1_file: str,
    ground_truth_hap2_file: str,
    outprefix="hcPhasing",
    profile_bin_size=100000,
    mode="heterozygous-only",
    is_clone=False,
) -> pd.DataFrame
```

This function evaluates **haplotype phasing accuracy** of inferred haplotype-specific CNAs (hap1/hap2) against ground-truth haplotype CNAs.

It converts both GT and tool haplotype CN profiles into a **binary phasing representation**, then computes:

- **mismatch error** (per cell or per clone)
- **switch error** (per cell or per clone)

It supports evaluation at **cell-level** (default) or **clone-level** (`is_clone=True`), and allows different evaluation modes (e.g., heterozygous-only).

------

## Parameters

| Name                     | Type        | Description                                                  |
| ------------------------ | ----------- | ------------------------------------------------------------ |
| `tool_hap1_cna_files`    | `List[str]` | List of tool CNA profile CSVs for inferred haplotype 1. Must align with `tool_names`. |
| `tool_hap2_cna_files`    | `List[str]` | List of tool CNA profile CSVs for inferred haplotype 2. Must align with `tool_names`. |
| `tool_names`             | `List[str]` | Tool names used in output rows.                              |
| `ground_truth_hap1_file` | `str`       | Path to ground-truth haplotype-1 CNA matrix CSV.             |
| `ground_truth_hap2_file` | `str`       | Path to ground-truth haplotype-2 CNA matrix CSV.             |
| `outprefix`              | `str`       | Prefix for the output CSV filename. Default: `"hcPhasing"`.  |
| `profile_bin_size`       | `int`       | Bin size used to split/standardize regions before alignment. Default: `100000` (100kb). |
| `mode`                   | `str`       | Evaluation mode passed to mismatch/switch error functions. Default: `"heterozygous-only"`. |
| `is_clone`               | `bool`      | If `True`, compute errors per clone label; otherwise per cell. Default: `False`. |

------

## Input File Format

### Ground truth files (`ground_truth_hap1_file`, `ground_truth_hap2_file`)

Loaded via `read_and_drop_empty(...)` and indexed by `region`:

```python
g1_r = read_and_drop_empty(ground_truth_hap1_file)
g2_r = read_and_drop_empty(ground_truth_hap2_file)
g1_r.set_index("region", inplace=True)
g2_r.set_index("region", inplace=True)
```

Expected structure:

- `region` column
- remaining columns: cell IDs (or clone labels if evaluating clone-level elsewhere in your pipeline)
- values: haplotype-specific CNA states (numeric or string convertible, depending on `_phase_to_binary`)

### Tool prediction files (`tool_hap1_cna_files`, `tool_hap2_cna_files`)

Also expected to have:

- `region` column
- cell/clone columns
- haplotype-specific CNA values

------

## Processing Steps

For each tool:

1. Read predicted haplotype profiles:

```python
t1 = read_and_drop_empty(f_h1)
t2 = read_and_drop_empty(f_h2)
```

1. Re-bin regions to `profile_bin_size`:

```python
t1 = split_all_regions(t1.set_index("region"), profile_bin_size)
t2 = split_all_regions(t2.set_index("region"), profile_bin_size)
```

1. Align predictions with GT on common regions and samples:

```python
g1, t1 = align(g1_r, t1)
g2, t2 = align(g2_r, t2)
```

1. Convert haplotype CNAs to a binary phasing representation:

```python
g1_bin, _ = self._phase_to_binary(g1, g2)
t1_bin, _ = self._phase_to_binary(t1, t2)
```

`_phase_to_binary` is expected to encode phasing comparisons bin-wise (e.g., indicating whether hap1 > hap2, hap1 < hap2, or equal), producing matrices suitable for mismatch/switch evaluation.

------

## Error Metrics

The function chooses evaluation routines based on `is_clone`:

- if `is_clone=False` (default):
  - `simu_cell_mismatch_error`
  - `simu_cell_switch_error`
- if `is_clone=True`:
  - `simu_clone_mismatch_error`
  - `simu_clone_switch_error`

Both mismatch and switch error functions are called with:

- predicted binary matrix `t1_bin`
- GT binary matrix `g1_bin`
- `mode` (e.g., `"heterozygous-only"`)

```python
mismatch_error_result = eval_mismatch_fn(t1_bin, g1_bin, mode)
switch_error_result   = eval_switch_fn(t1_bin, g1_bin, mode)
```

The two result tables are merged on:

- `cell` (cell-level), or
- `clone_label` (clone-level)

```python
label = "clone_label" if is_clone else "cell"
result = pd.merge(mismatch_error_result, switch_error_result, on=label, how="outer")
result["tool_name"] = name
```

------

## Output

A single CSV is written to `self.output_dir`:

```
os.path.join(self.output_dir, f"{outprefix}_{mode}.csv")
```

Default example:

```
{self.output_dir}/hcPhasing_heterozygous-only.csv
```

The output is **long-format**, with one row per cell/clone (depending on `is_clone`) per tool.

------

## Output Table Schema

Exact columns depend on what your mismatch/switch routines return, but typically include:

| Column                   | Meaning                                       |
| ------------------------ | --------------------------------------------- |
| `cell` or `clone_label`  | Identifier used for per-sample evaluation     |
| mismatch-related columns | e.g., mismatch count / mismatch ratio         |
| switch-related columns   | e.g., switch error count / switch error ratio |
| `tool_name`              | Tool name                                     |

------

## Return Value

Returns a `pd.DataFrame` containing the merged mismatch + switch error results for all tools.

------

## Example

```python
from hcbench.gtbench.gtbench import GTBench

bench = GTBench(output_dir="out/gt_output")

df = bench.hcPhasing(
    tool_hap1_cna_files=[
        "/path/to/chisel/hap1.csv",
        "/path/to/signals/hap1.csv",
    ],
    tool_hap2_cna_files=[
        "/path/to/chisel/hap2.csv",
        "/path/to/signals/hap2.csv",
    ],
    tool_names=["CHISEL", "SIGNALS"],
    ground_truth_hap1_file="/path/to/gt/hap1.csv",
    ground_truth_hap2_file="/path/to/gt/hap2.csv",
    outprefix="hcPhasing",
    profile_bin_size=100000,
    mode="heterozygous-only",
    is_clone=False,
)

print(df.head())
```

