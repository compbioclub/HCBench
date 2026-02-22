from abc import ABC, abstractmethod
import os
import pandas as pd

class BaseParser(ABC):

    def __init__(self, input_path: str, output_path: str,change_hap: bool = False):
        self.input_path = input_path
        self.output_path = output_path
        self.change_hap = change_hap

    @abstractmethod
    def run(self):
        pass

    def _read(self, sep="\t") -> pd.DataFrame:
        return pd.read_csv(self.input_path, sep=sep)

    def _check_output_path(self):
        output_dir = os.path.dirname(self.output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"[hcbench] Created output directory: {output_dir}")
        else:
            print(f"[hcbench] Output directory verified: {output_dir}")