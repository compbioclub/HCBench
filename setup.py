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
        ],
    },
    install_requires=["pandas>=2.0.0"]
)
