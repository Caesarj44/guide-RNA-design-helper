import re
from dataclasses import dataclass

import pandas as pd

from core.doench16_scorer import predict_doench16


# test_seq = 'TGGCAGGATATATTGTGGTGTAAACAAATTGACGCTTAGACAACTTAATAACACATTGCGGACGTTTTTAATGTACTGAATTAACGCCGAATTAATTCGGGGGATCTGGATTTTAGTACTGGATTTTGGTTTTAGGAATTAGAAATTTTATTGATAGAAGTATTTTACAAATACAAATACATACTAAGGGTTTCTTATATGCTCAACACATGAGCGAAACCCTATAGGAACCCTAATTCCCTTATCTGGGAACTACTCACACATTATTATGGAGAAACTCGAAGATCCGTCGAGCTTGTCGATCGACAGATCCGGTCGGCATCATAACTTCGTATAGCATACATTATACGAAGTTAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGCTAGAGCAGCTTGAGCTTGGATCAGATTGTCGTTTCCCGCCTTCAGTTTAAACTATCAGTGTTTGACAGGATATATTGGCGGGTAAACCTAAGAGAAAAGAGCGTTTA'
# target_seq = test_seq
# #pam_seq = 'NNGRRT'
# pam_seq = 'NGG'
# spacer_len = 20

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

COMPLINETARY_MAP = {
    'A' : 'T',
    'T' : 'A',
    'G' : 'C',
    'C' : 'G',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'K' : 'M',
    'M' : 'K',
    'B' : 'V',
    'V' : 'B',
    'D' : 'H',
    'H' : 'D',
    'N' : 'N' 
}

@dataclass
class LocalGRNARecord:
    index : int
    spacer : str
    pam : str
    start : int
    end : int
    lenght : int
    strand : str
    chrom : str
    doench16 : float
    doench14 : float
    sequence_with_flanking : str
    distance_to_cut_from_end : int



def find_grna_in_sequence (
        pam_seq: str,
        spacer_len: int,
        target_seq: str,
        offset : int = 0,
        chrom : str = None
):
    pam_len = len(pam_seq)
    pam_forw = ''
    pam_rev = ''
    target_seq = target_seq.upper()

    for letter in pam_seq:
        pam_forw += IUPAC_MAP[letter]
    #print(pam_forw)

    for letter in pam_seq: 
        pam_rev = IUPAC_MAP[COMPLINETARY_MAP[letter]] + pam_rev
    #print(pam_rev)

    regular_pam_forw = re.compile(rf'(?=([ACGT]{{{spacer_len}}}{pam_forw}))')
    regular_pam_rev = re.compile(rf'(?=({pam_rev}[ATGC]{{{spacer_len}}}))')
    records = []

    counter = 0
    for i in regular_pam_forw.finditer(target_seq):
        counter += 1
        founded_seq = i.group(1)
        start_position = i.start()
        current_pam = founded_seq[spacer_len:spacer_len+pam_len]
        current_spacer = founded_seq[0:spacer_len]
        current_sequence_with_flanking = target_seq[offset+start_position-4:offset+start_position+26]
        if len(current_sequence_with_flanking) != 30:
            current_sequence_with_flanking = None
        #print(current_sequence_with_flanking, len(current_sequence_with_flanking))
        #print(f'----{founded_seq}---')

        if current_sequence_with_flanking is not None:
            current_doench16 = round(100*predict_doench16(current_sequence_with_flanking))
        else:
            current_doench16 = 'NotEnoughFlankSeq'

        rec = _create_new_grna_record(
            index = counter,
            spacer=current_spacer,
            pam=current_pam,
            chrom=chrom,
            start=offset + start_position,
            end=offset + start_position  + spacer_len,
            strand='+',
            lenght=pam_len+spacer_len,
            distance_to_cut_from_end= 3,
            doench16=current_doench16,
            sequence_with_flanking = current_sequence_with_flanking
        )
        records.append(rec)
    #print('reverse str')
    for i in regular_pam_rev.finditer(target_seq):
        counter +=1
        founded_seq = i.group(1)
        start_position = i.start()
        current_pam = _reverse_compliment(founded_seq[:pam_len])
        current_spacer = _reverse_compliment(founded_seq[pam_len:])
        current_sequence_with_flanking = _reverse_compliment(target_seq[offset+start_position-3:offset+start_position+27])
        if len(current_sequence_with_flanking) != 30:
            current_sequence_with_flanking = None
        #print(current_sequence_with_flanking, len(current_sequence_with_flanking))
        #print(f'----{founded_seq}---')
        if current_sequence_with_flanking is not None:
            current_doench16 = round(100*predict_doench16(current_sequence_with_flanking))
        else:
            current_doench16 = 'NotEnoughFlankSeq'

        rec = _create_new_grna_record(
            index=counter,
            spacer=current_spacer,
            pam=current_pam,
            chrom=chrom,
            start=offset + start_position + pam_len,
            end=offset + start_position  + spacer_len + pam_len,
            strand='-',
            lenght=pam_len+spacer_len,
            distance_to_cut_from_end=3,
            doench16=current_doench16,
            sequence_with_flanking = current_sequence_with_flanking
        )
        records.append(rec)

    records.sort(key=lambda r: r.start)
    #print(records)
    records_df = pd.DataFrame(records)
    return records_df


def _reverse_compliment(sequence):
    reversed_seq = ''
    for i in sequence:
        reversed_seq = COMPLINETARY_MAP[i] + reversed_seq
    return reversed_seq

def _create_new_grna_record(
    index : int,
    spacer : str,
    pam : str,
    start : int,
    end : int,
    strand : int,
    chrom : str,
    lenght : int,
    distance_to_cut_from_end,
    sequence_with_flanking : str,
    doench16: float = 0.0,
    doench14:float = 0.0,
    

):
    grna_rec = LocalGRNARecord(
        index=index,
        spacer=spacer,
        pam=pam,
        start=start,
        end=end,
        strand=strand,
        chrom=chrom,
        lenght=lenght,
        distance_to_cut_from_end=distance_to_cut_from_end,
        doench16=doench16,
        doench14=doench14,
        sequence_with_flanking=sequence_with_flanking
    )
    return grna_rec

# find_grna_in_sequence(
#     pam_seq=pam_seq,
#     spacer_len=spacer_len,
#     target_seq=target_seq
#     )