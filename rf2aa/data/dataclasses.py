from omegaconf import DictConfig
from dataclasses import dataclass
from typing import Protocol
import mmap

@dataclass
class FFindexDB:
    index: list
    data: mmap

class ModelRunner(Protocol):
    config: DictConfig
    ffdb: FFindexDB
    deterministic: bool
