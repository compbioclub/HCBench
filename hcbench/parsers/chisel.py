from .cna_parser_base import CNAParser

class ChiselParser(CNAParser):
    chrom_col = "#CHR"
    start_col = "START"
    end_col = "END"
    cell_col = "CELL"
    value_col = "CORRECTED_HAP_CN"
    start_plus_one = True
    add_chr_prefix = False
