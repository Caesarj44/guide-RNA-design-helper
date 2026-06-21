import os
import re
import time
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
from ViennaRNA import RNA
import cairosvg
from PySide6.QtCore import QObject, Signal, QThread 


@dataclass
class GRNARecord:
    # Основные параметры gRNA
    sequence: str           # 20-нт спейсер (без PAM)
    pam: str                # PAM (4 нт: NGGN из CRISPOR)
    strand: str = ''        # 'sense' | 'antisense'
    doench16: float = 0.0   # Doench'16 score из CRISPOR (0–100)
    low_doench: bool = False # True если doench16 < порога (по умолчанию 50)

    # Скоринговые поля
    gc_freq: int = 0            # % GC всего спейсера
    gc_count_10to20: int = 0    # кол-во GC на позициях 10–20
    last_four_pur: int = 0      # кол-во пуринов на позициях 17–20 (0..4)
    access_18_to_20: int = 0    # взвешенная доступность позиций 18,19,20 (0..4; поз.20 весит 2)
    c_not_at_3: float = 0.0     # штраф/бонус: C нежелателен в позиции 3
    g_not_at_16: float = 0.0
    c_at_16: float = 0.0
    g_or_a_at_20: float = 0.0
    c_at_18: float = 0.0
    g_not_at_14: float = 0.0
    pam_n_score: float = 0.0    # бонус/штраф за N в PAM (первый нт PAM)

    # Фильтрующие флаги (666 = дисквалификация)
    oligo_t: int = 0            # 666 если есть TTTT
    gcc_at_16to20: int = 0      # 666 если GCC на позициях 16–20
    seven_links: int = 0        # 666 если ≥8 подряд связанных нт в структуре
    twelve_links: int = 0       # 666 если >12 связанных нт всего в структуре

    # Новые поля
    tgg_flag: bool = False      # True если спейсер начинается с TGG
    restriction_sites: list = field(default_factory=list)  # список рестриктаз

    # PNG вторичной структуры (путь к файлу)
    structure_png_path: str = ''
    # dot-bracket строка вторичной структуры
    structure_dot_bracket: str = ''
    # MFE (ккал/моль)
    structure_mfe: float = 0.0

    @property
    def is_disqualified (self) -> bool:
            return (self.oligo_t == 666 or
                    self.gcc_at_16to20 == 666 or
                    self.seven_links == 666 or
                    self.twelve_links == 666)


def gc_structure(sequence: str) -> tuple[int,int,int]:
    gc_count = sequence.count('G') + sequence.count('C')
    gc_freq = int(round(gc_count / len(sequence) * 100, 0))
    gc_10_20 = sequence[9:20].count('G') + sequence[9:20].count('C')
    last_four_pur = sequence[16:20].count('A') + sequence[16:20].count('G')
    return gc_10_20, gc_freq, last_four_pur


def compute_secondary_structure(
    spacer: str,
    grna_tail: str,
    output_dir: str | Path,
    index: int,
) -> tuple[str, str, float]:
    full_seq = spacer + grna_tail
    fc = RNA.fold_compound(full_seq)
    structure, mfe = fc.mfe()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = str(output_dir / f'{index}.png')
    
    with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as tmp:
        svg_path = tmp.name
    RNA.svg_rna_plot(full_seq,structure, svg_path)
    time.sleep(0.05)
    cairosvg.svg2png(url=svg_path, write_to=png_path)
    
    if os.path.exists(svg_path):
        os.remove(svg_path)

    return(png_path, structure, float(mfe))

compute_secondary_structure(grna_tail="GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT",
                            spacer="ACACTGATAGTTTAAACTGAAGG",
                            output_dir='output/pngs',
                            index=0)



DEFAULT_CONFIGS = {
    'min_gc': 40,
    'max_gc': 70,
    'doench_threshold': 50,
    'forv_strain': 0,
    'rev_strain': 1,
    'last_four_pur_1': 0.25,
    'last_four_pur_2': 0.5,
    'last_four_pur_3': 0.75,
    'last_four_pur_4': 1.0,
    'access_to_20': 2,
    'access_to_19': 1,
    'access_to_18': 1,
    'C_not_at_3': 0.5,
    'G_not_at_16': 0.5,
    'C_at_16': 0.5,
    'G_or_A_at_20': 0.5,
    'C_at_18': 0.5,
    'G_not_at_14': 0.5,
    'PAMs_A': 0,
    'PAMs_T': -1,
    'PAMs_G': 1,
    'PAMs_C': 1,
}

def analyze_single_grna(
    spacer: str,
    pam: str,
    guide_id: str,
    doench16: float,
    structure_dot: str,
    png_path: str,
    mfe: float,
    configs: dict,
    enzymes: list,
    doench_threshold: float = 50.0,
) -> GRNARecord:
    
    rec = GRNARecord(
        sequence=spacer,
        pam=pam,
        doench16=doench16,
        low_doench=(doench16 < doench_threshold),
        structure_dot_bracket=structure_dot,
        structure_png_path=png_path,
        structure_mfe=mfe,
    )

    # Цепь прямая или обратная
    if str(guide_id).endswith('forw'):
        rec.strand = 'sense'
    elif str(guide_id).endswith('rev'):
        rec.strand = 'antisense'

    # GC состав
    rec.gc_count_10to20, rec.gc_freq, rec.last_four_pur = gc_structure(spacer)

    # TGG флаг
    rec.tgg_flag = spacer.startswith('TGG')

    # Олиго-T
    if 'TTTT' in spacer:
        rec.oligo_t = 666

    # GCC на позициях 16–20 (индексы 15:20)
    if 'GCC' in spacer[15:20]:
        rec.gcc_at_16to20 = 666

    # Позиционные критерии
    if spacer[2] != 'C':
        rec.c_not_at_3 = configs['C_not_at_3']
    if spacer[15] != 'G':
        rec.g_not_at_16 = configs['G_not_at_16']
    if spacer[15] == 'C':
        rec.c_at_16 = configs['C_at_16']
    if spacer[19] in ('G', 'A'):
        rec.g_or_a_at_20 = configs['G_or_A_at_20']
    if spacer[17] == 'C':
        rec.c_at_18 = configs['C_at_18']
    if spacer[13] != 'G':
        rec.g_not_at_14 = configs['G_not_at_14']

    # PAM N-нуклеотид (первый нт PAM)
    pam_n = pam[0] if pam else ''
    if pam_n == 'G':
        rec.pam_n_score = configs['PAMs_G']
    elif pam_n == 'C':
        rec.pam_n_score = configs['PAMs_C']
    elif pam_n == 'T':
        rec.pam_n_score = configs['PAMs_T']
    elif pam_n == 'A':
        rec.pam_n_score = configs['PAMs_A']

    # Доступность позиций 18, 19, 20 (индексы 17, 18, 19 в dot-bracket)
    if structure_dot and len(structure_dot) >= 20:
        if structure_dot[19] not in ('(', ')'):
            rec.access_18_to_20 += configs['access_to_20']
        if structure_dot[18] not in ('(', ')'):
            rec.access_18_to_20 += configs['access_to_19']
        if structure_dot[17] not in ('(', ')'):
            rec.access_18_to_20 += configs['access_to_18']

    # Структурные фильтры (по первым 20 нт dot-bracket)
    if structure_dot and len(structure_dot) >= 20:
        spacer_dot = structure_dot[:20]
        if '((((((((' in spacer_dot:
            rec.seven_links = 666
        paired_count = spacer_dot.count('(') + spacer_dot.count(')')
        if paired_count > 12:
            rec.twelve_links = 666

    # Сайты рестрикции
    rec.restriction_sites = check_restriction_sites(spacer, enzymes)

    return rec

