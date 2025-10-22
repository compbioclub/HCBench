from .cna_parser_base import CNAParser

class SignalsParser(CNAParser):
    chrom_col = "chr"
    start_col = "start"
    end_col = "end"
    cell_col = "cell_id"
    value_col = "state_AS_phased"
    start_plus_one = False
    add_chr_prefix = True
