import csv
import re
from pathlib import Path
from functools import lru_cache


IUPAC_MAP = {
    'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
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

# Паттерн нотации вырожденного сайта разреза: (16/14) или (-5/-1)
_CLEAVAGE_NOTATION = re.compile(r'\(\s*-?\d+\s*/\s*-?\d+\s*\)')


def _iupac_to_regex(seq: str) -> str:
    return ''.join(IUPAC_MAP.get(c, c) for c in seq.upper())


def _parse_recognition_seq(raw: str) -> str | None:
    seq = raw.strip()
    if seq.startswith('('):
        return None
    seq = _CLEAVAGE_NOTATION.sub('', seq)
    seq = seq.replace('^', '')
    seq = seq.strip()
    if not seq:
        return None
    valid = set(IUPAC_MAP.keys())
    if not all(c in valid for c in seq.upper()):
        return None

    return seq.upper()


def load_restriction_db(csv_path: str | Path) -> list[dict]:
    enzymes = []
    skipped = []
    csv_path = Path(csv_path)

    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['Name'].strip()
            raw_seq = row['Seq'].strip()

            consensus = _parse_recognition_seq(raw_seq)
            if consensus is None:
                skipped.append(name)
                continue

            regex_str = _iupac_to_regex(consensus)
            try:
                pattern = re.compile(regex_str)
            except re.error:
                skipped.append(name)
                continue

            enzymes.append({
                'name': name,
                'consensus': consensus,
                'pattern': pattern,
            })

    return enzymes


def check_restriction_sites(spacer: str, enzymes: list[dict]) -> list[str]:
    spacer_upper = spacer.upper()
    found = []
    for enzyme in enzymes:
        if enzyme['pattern'].search(spacer_upper):
            found.append(enzyme['name'])
    return found



@lru_cache(maxsize=1)
def get_default_restriction_db() -> list[dict]:
    assets_dir = Path(__file__).parent.parent / 'assets'
    csv_path = assets_dir / 'restriction_parsing_result.csv'
    return load_restriction_db(csv_path)
