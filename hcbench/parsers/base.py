from abc import ABC, abstractmethod
import pandas as pd

class BaseParser(ABC):

    def __init__(self, input_path: str, output_path: str):
        self.input_path = input_path
        self.output_path = output_path

    @abstractmethod
    def run(self):
        pass

    def _read(self, sep="\t") -> pd.DataFrame:
        return pd.read_csv(self.input_path, sep=sep)

    def _save(self, df: pd.DataFrame):
        df.to_csv(self.output_path)
        print(f"[hcbench] {self.__class__.__name__} → {self.output_path}")
