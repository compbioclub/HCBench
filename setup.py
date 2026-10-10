from setuptools import setup, find_packages

setup(
    name="hcbench",
    version="0.1.0",
    description="",
    author="",
    packages=find_packages(include=["hcbench", "hcbench.*"]),
    include_package_data=True, 
    package_data={
        "hcbench": [
            "ref/hg19_centromeres.csv",
            "ref/hg38_centromeres.csv",
            "lib/pathways/MSigDB_Hallmark_2020.gmt",
            "lib/pathways/Reactome_2022.gmt",
            "lib/pathways/KEGG_2021_Human.gmt",
        ],
    },
    python_requires=">=3.10",
    install_requires=[
        "numpy>=1.22",
        "pandas>=2.0,<3",
        "scipy>=1.8",
        "scikit-learn>=1.0",
        "biopython>=1.79",
    ],
    extras_require={
        "rds": ["pyreadr>=0.4.6"],
        "enrichment": ["gseapy>=1.1.9"],
    },
)
