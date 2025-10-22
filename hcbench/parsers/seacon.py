from .cna_parser_base import CNAParser

class SeaconParser(CNAParser):
    chrom_col = "chrom"
    start_col = "start"
    end_col = "end"
    cell_col = "cell"
    value_col = "CN"
    start_plus_one = False
    add_chr_prefix = False

    def preprocess_value(self, value):
        if isinstance(value, str):
            return value.replace(",", "|")
        return value
