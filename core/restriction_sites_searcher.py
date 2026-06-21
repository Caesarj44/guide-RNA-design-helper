import csv
import re
from pathlib import Path
from functools import lru_cache

# IUPAC ambiguity codes → regex character classes
IUPAC_MAP = {
    'A': 'A',
    'T': 'T',
    'G': 'G',
    'C': 'C',
    'R': '[AG]',   # puRine
    'Y': '[CT]',   # pYrimidine
    'S': '[GC]',   # Strong
    'W': '[AT]',   # Weak
    'K': '[GT]',   # Keto
    'M': '[AC]',   # aMino
    'B': '[CGT]',  # not A
    'D': '[AGT]',  # not C
    'H': '[ACT]',  # not G
    'V': '[ACG]',  # not T
    'N': '[ACGT]', # aNy
}

