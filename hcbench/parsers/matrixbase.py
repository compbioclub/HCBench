import numpy as np
import pandas as pd
from .cna_parser_base import CNAParser
from .utils import read_table_auto, split_all_regions

class MatricBaseParser(CNAParser):

    start_offset: int = 0
    add_chr_prefix = False

    def before_pivot(self):
        cna_df = read_table_auto(self.input_path)
        cna_df.set_index("bin", inplace=True)
        cna_df.index.name = "region"

        print(cna_df.head())

        return cna_df

    

