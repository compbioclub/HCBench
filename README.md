# HCBench

Welcome to the official documentation of HCBench.

## 🗺️ Overview

<img src="docs/assets/HCDSIM_Figure1_3.png" style="zoom: 50%;" />

`hcbench` is organized around three main components:

1. [**HCBench-parser**](https://compbioclub.github.io/HCBench/parsers/): standardizing CNA results from different callers
2. [**HCBench-sim**](https://compbioclub.github.io/HCBench/gtbench/): evaluating detection accuracy against simulated ground truth datasets from HCDSIM
3. [**HCBench-real**](https://compbioclub.github.io/HCBench/realbench/): assessing pairwise consistency between callers using large-scale real-world datasets

## 📦 Installation

```
pip install git+https://github.com/compbioclub/HCBench.git@main
```

Python 3.10 or newer is required. Core dependencies are NumPy, pandas
(`>=2.0,<3`), SciPy, scikit-learn, and Biopython; they are declared in `setup.py`.
The declared lower bounds are based on the APIs used, not a tested version matrix.

Optional features can be installed separately:

```sh
# Alleloscope RDS input
pip install "hcbench[rds] @ git+https://github.com/compbioclub/HCBench.git@main"
# Focal gene and pathway enrichment
pip install "hcbench[enrichment] @ git+https://github.com/compbioclub/HCBench.git@main"
```

VCF conversion requires the external `bcftools` executable. Follow the
[BCFtools installation instructions](https://samtools.github.io/bcftools/howtos/install.html)
and put `bcftools` on `PATH`. CNA phylogeny analysis requires the external
[Lazac copy-number tool](https://github.com/raphael-group/lazac-copy-number);
follow its build/install instructions and put `lazac` on `PATH`. HCBench uses
`bcftools query` and `lazac nni`; neither executable is installed by pip extras.
Record their versions (or source commits) for reproducibility. Their compatibility
and the Python dependency matrix still require runtime validation.



## 📌 Citation

If you use HCBench in your research, please cite the following paper:

```

```

