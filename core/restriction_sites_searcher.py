import csv
import re
from pathlib import Path
from functools import lru_cache


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
    'N': '[ATGC]', # aNy
}

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


def check_restriction_sites(
        sequence: str,
        enzymes: list[dict],
        target_seq : str,
        radius = [100,250,500]
        ) -> list[str]:
    sequence_upper = sequence.upper()
    found = []
    print('тут ищутся сайты рестрикции',sequence_upper)
    for enzyme in enzymes:
        
        if enzyme['pattern'].search(sequence_upper):
            current_enzyme = enzyme['name']

            start = target_seq.find(sequence_upper)
            center = start + len(sequence_upper) // 2
            target_seq_len = len(target_seq)

            current_enzyme_w_count = current_enzyme + '('
            for i, rad in enumerate(radius):
                window_start = max(0, center - rad)
                window_end = min(target_seq_len, center + rad)
                window = target_seq[window_start:window_end]
                current_count = len(enzyme['pattern'].findall(window))
                current_enzyme_w_count += str(current_count)
                if i + 1 != len(radius):
                    current_enzyme_w_count += '/'

            current_enzyme_w_count += ')'
            found.append(current_enzyme_w_count)

    return found

@lru_cache(maxsize=1)
def get_default_restriction_db() -> list[dict]:
    assets_dir = Path(__file__).parent.parent / 'assets'
    csv_path = assets_dir / 'restriction_parsing_result.csv'
    return load_restriction_db(csv_path)
