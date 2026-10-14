# HCDSIM ground truth (GT)

Use `GTParser` to convert an HCDSIM ground-truth CNA profile into region-by-cell matrices, then annotate and count CNA events.

## Prepare the HCDSIM profile

The input is a CSV with `Chromosome`, `Start`, and `End` columns followed by one column per cell. Each row is a genomic region; each cell entry is a haplotype copy-number pair separated by `|`. For example:

```csv
Chromosome,Start,End,normal_cell1,normal_cell2
chr1,1,100000,1|1,1|1
chr1,100001,200000,1|1,1|1
chr1,200001,300000,2|1,1|1
```

Use `chr`-prefixed numeric chromosomes such as `chr1`; the current region sorter expects this form. Every data row should have the same number of fields as the header. Match the input coordinates to the `genome` you select.

## Parse the profile

`bin_size` is a required positive integer in base pairs, and `genome` is required as `"hg19"` or `"hg38"`. Create the output directory, then run the parser:

```python
from pathlib import Path
from hcbench.parsers.gt import GTParser

gt_file = "input/cell_cna_profile.csv"
output_dir = Path("output/gt")
output_dir.mkdir(parents=True, exist_ok=True)

gt_parser = GTParser(
    input_path=gt_file,
    output_path=str(output_dir),
    bin_size=100000,
    genome="hg38",
)
result = gt_parser.run()
```

Initialization prints the bin size, genome, and whether regions will be split. With the default `expand_bins=False`, regions are not split. Set `expand_bins=True` in the constructor to split regions into pieces of at most `bin_size` bp. `run()` writes `haplotype_combined.csv`, `haplotype_1.csv`, `haplotype_2.csv`, `minor.csv`, `major.csv`, and `minor_major.csv` in `output/gt/`; `result` is `None`. The combined CSV has a `region` column, such as `chr1:1-100000`, followed by cell columns.

## Annotate CNA events

```python
gt_parser.annotate_segments(
    str(output_dir / "haplotype_combined.csv"),
    detect_cnv_type=True,
    out_file="ground_truth_classified.csv",
)
```

This writes `output/gt/ground_truth_classified.csv` and returns `None`. With `detect_cnv_type=True`, the CSV includes `cell`, `Chrom`, `Start`, `End`, `value`, `Segment`, `size`, and `type`. The default annotation threshold is `0.8`. If you omit `detect_cnv_type=True`, the default is `False` and the output has no `type` column. Use `chromosomes=["chr1", "chr2"]` to annotate selected chromosomes; the parser also writes a CSV for each selected chromosome.

With the default size thresholds, focal events are shorter than 3 Mb, medium events are 3–10 Mb, and broad events are longer than 10 Mb. The configured `bin_size` limits which categories are analyzed:

| `bin_size` | Categories included |
| --- | --- |
| Up to 3 Mb | focal, medium, broad |
| Above 3 Mb and up to 10 Mb | medium, broad |
| Above 10 Mb | broad |

Categories outside that range are omitted from the annotation, rather than shown with a count of zero. You can set `small_size` and `mid_size` on `annotate_segments()` if your analysis needs different size thresholds.

