from dataclasses import dataclass,field
from pathlib import Path

import re

from ViennaRNA import RNA
import pandas as pd

from PySide6.QtCore import(
    QObject,
    Signal,
)
from core.restriction_sites_searcher import check_restriction_sites, get_default_restriction_db
from core.RNA_structure_plot import plot_rna

@dataclass
class GRNA_data(): 
    sequence : str 
    PAM : str
    strand : int = 0
    gc_freq : int = 0 
    doench_16 : float = 0.0

    png_2d_path : str = ''
    dot_bracked_structure : str = ''
    mfe : float = 0.0
    mfe_bad : bool = False
    oligoT : bool = False
    oligoG : bool = False
    oligoC : bool = False
    gc_count_10to20 : float = 0.0
    last_four_pur : float = 0.0
    access_18_to_20 : float = 0.0
    access_51_to_53: float = 0.0
    eight_links : float = 0.0
    twelve_links : float = 0.0
    C_not_at_3 : float = 0.0
    G_not_at_16 : float = 0.0
    C_at_16 : float = 0.0
    G_or_A_at_20 : float = 0.0
    C_at_18 : float = 0.0 
    G_not_at_14 : float = 0.0
    GCC_not_at_16to20 : float = 0.0
    PAMs_N : int = 0
    restriction_sites : list = field(default_factory=list)
    total_score : float = 0.0


GRNA_DEFAULT_CONFIGS = {
    'mfe_good' : 1,
    'mfe_bad' : True,
    'min_gc': 40,
    'max_gc': 70,
    'gc_in_range': 2,
    'oligoT': True,
    'oligoG' : True,
    'oligoC' : True,
    'forv_strain': 0,
    'rev_strain': 1,
    'GC_freq': 0,
    'GC_count_10to20': 0,
    'last_four_pur_1': 0.25,
    'last_four_pur_2': 0.5,
    'last_four_pur_3': 0.75,
    'last_four_pur_4': 1,
    'access_to_20': 4,
    'access_to_19': 2,
    'access_to_18': 2,
    'access_to_51' :2,
    'access_to_52' :2,
    'access_to_53' :2,
    'eight_links': True,
    'twelve_links': True,
    'C_not_at_3': 0.5,
    'G_not_at_16': 0.5,
    'C_at_16': 0.5,
    'G_or_A_at_20': 0.5,
    'C_at_18': 0.5,
    'G_not_at_14': 0.5,
    'GCC_not_at_16to20': True,
    'PAMs_A': 0,
    'PAMs_T': -1,
    'PAMs_G': 1,
    'PAMs_C': 1
}


def grna_analyzer(
    spacer: str,
    pam: str,
    guide_id: str,
    strand: str,
    doench16: float,
    structure_dot: str,
    png_path : str,
    mfe: float,
    configs: dict,
    enzymes: list,
    all_target_seq : str,
    doench_threshold: float = 50.0,
):
    grna_record = GRNA_data(
        sequence= spacer,
        PAM=pam,
        doench_16= doench16,
        png_2d_path = png_path,
        dot_bracked_structure=structure_dot,
        mfe=mfe,

        
    )
    print(grna_record.doench_16, 'grna_record')
    if strand == 'forv':
        grna_record.strand = configs['forv_strain']
        grna_record.total_score += grna_record.strand
    elif strand =='rev':
        grna_record.strand = configs['rev_strain']
        grna_record.total_score += grna_record.strand
        
    grna_record.gc_count_10to20, grna_record.gc_freq, grna_record.last_four_pur = _gc_structure(spacer)

    if configs['min_gc'] <= grna_record.gc_freq * 100 <= configs['max_gc']:
        grna_record.total_score += configs['gc_in_range']

    ## need to compute optimal mfe borders
    # print(grna_record.mfe, 'minimal free energy')
    # if grna_record.mfe > -34:
    #     grna_record.total_score += configs['mfe_good']
    # elif grna_record.mfe < -36:
    #     grna_record.mfe_bad = configs['mfe_bad']

    if 'TTTT' in grna_record.sequence[0:20]:
        grna_record.oligoT = configs['oligoT']
    if 'GGGG' in grna_record.sequence[0:20]:
        grna_record.oligoG = configs['oligoG']
    # if 'CCCC' in grna_record.sequence[0:20]:
    #     grna_record.oligoC = configs['oligoC']
    
    if grna_record.last_four_pur == 1:
        grna_record.last_four_pur = configs['last_four_pur_1']
        grna_record.total_score += grna_record.last_four_pur
    elif grna_record.last_four_pur == 2:
        grna_record.last_four_pur = configs['last_four_pur_2']
        grna_record.total_score += grna_record.last_four_pur
    elif grna_record.last_four_pur == 3:
        grna_record.last_four_pur = configs['last_four_pur_3']
        grna_record.total_score += grna_record.last_four_pur
    elif grna_record.last_four_pur == 4:
        grna_record.last_four_pur = configs['last_four_pur_4']
        grna_record.total_score += grna_record.last_four_pur
        
        
    if grna_record.sequence[2:3] != 'C': #цитозин не желателен в 3 положении
        grna_record.C_not_at_3 = configs['C_not_at_3']   
        grna_record.total_score += configs['C_not_at_3'] 

    if grna_record.sequence[15:16] != 'G': #гуанин нежелателен в 16 положении
        grna_record.G_not_at_16 = configs['G_not_at_16'] 
        grna_record.total_score += configs['G_not_at_16']  

    if grna_record.sequence[15:16] == 'C': #в 16 положении желателен цитозин
        grna_record.C_at_16 = configs['C_at_16'] 
        grna_record.total_score += configs['C_at_16']  

    if grna_record.sequence[19:20] == 'G' or grna_record.sequence[19:20] == 'A': #аденин или гуанин желателен в 20 положении
        grna_record.G_or_A_at_20 = configs['G_or_A_at_20']
        grna_record.total_score += configs['G_or_A_at_20'] 

    if grna_record.sequence[17:18] == 'C': #цитозин желателен в 18 положении
        grna_record.C_at_18 = configs['C_at_18'] 
        grna_record.total_score += configs['C_at_18']  

    if grna_record.sequence[13:14] != 'G': #гуанин не желателен в 14 положении
        grna_record.G_not_at_14 = configs['G_not_at_14']
        grna_record.total_score += configs['G_not_at_14'] 

    if 'GCC' in grna_record.sequence[15:20] : #GCC не должно быть на участке 16-20
        grna_record.GCC_not_at_16to20 = configs['GCC_not_at_16to20']
        
    if grna_record.dot_bracked_structure[19] != '(' and grna_record.dot_bracked_structure[19] != ')': #20 нуклеотид очень желательно не должен быть связан
        grna_record.access_18_to_20 += configs['access_to_20'] 
        grna_record.total_score += configs['access_to_20']            

    if grna_record.dot_bracked_structure[18] != '(' and grna_record.dot_bracked_structure[18] != ')':#19 нуклеотид желательно не должен быть связан
        grna_record.access_18_to_20 += configs['access_to_19']
        grna_record.total_score += configs['access_to_19']

    if grna_record.dot_bracked_structure[17] != '(' and grna_record.dot_bracked_structure[17] != ')':#18 нуклеотид желательно не должен быть связан
        grna_record.access_18_to_20 += configs['access_to_18']
        grna_record.total_score += configs['access_to_18']


    if grna_record.dot_bracked_structure[50] != '(' and grna_record.dot_bracked_structure[50] != ')': #51 нуклеотид очень желательно не должен быть связан
        grna_record.access_51_to_53 += configs['access_to_51'] 
        grna_record.total_score += configs['access_to_51']           

    if grna_record.dot_bracked_structure[51] != '(' and grna_record.dot_bracked_structure[51] != ')':#52 нуклеотид желательно не должен быть связан
        grna_record.access_51_to_53 += configs['access_to_52']
        grna_record.total_score += configs['access_to_52']

    if grna_record.dot_bracked_structure[52] != '(' and grna_record.dot_bracked_structure[52] != ')':#53 нуклеотид желательно не должен быть связан
        grna_record.access_51_to_53 += configs['access_to_53']
        grna_record.total_score += configs['access_to_53']


    if '((((((((' in grna_record.dot_bracked_structure[0:20]: #проверка нет ли 8 связанных подряд
        grna_record.eight_links = configs['eight_links']

    if grna_record.dot_bracked_structure[0:20].count('(') + grna_record.dot_bracked_structure[0:20].count(')') >12: #проверка нет ли более 12 связанных вообще
        grna_record.twelve_links = configs['twelve_links']  

    if grna_record.PAM[0] == 'G':
        grna_record.PAMs_N = configs['PAMs_G']
        grna_record.total_score += grna_record.PAMs_N
    elif grna_record.PAM[0] == 'C':
        grna_record.PAMs_N = configs['PAMs_C']
        grna_record.total_score += grna_record.PAMs_N
    elif grna_record.PAM[0] == 'T':
        grna_record.PAMs_N = configs['PAMs_T']
        grna_record.total_score += grna_record.PAMs_N
    elif grna_record.PAM[0] == 'A':
        grna_record.PAMs_N = configs['PAMs_A']
        grna_record.total_score += grna_record.PAMs_N

    near_spacer_end = found_near_spacer_end_sequence(spacer=spacer,target_seq=all_target_seq)
    grna_record.restriction_sites = check_restriction_sites(spacer=near_spacer_end,enzymes=enzymes)

    return grna_record

def found_near_spacer_end_sequence(spacer, target_seq):
    file_path = 'guides_2xLoxPafterCre_pz9Athaliana-unknownLoc.xls'
    crispor_df = pd.read_excel(file_path)

    pair = {
        'A' : 'T',
        'T' : 'A',
        'U' : 'A',
        'G' : 'C',
        'C' : 'G'
    }
    spacer_rev = ''
    target_seq = crispor_df.iloc[0,1]
    target_seq_len = len(target_seq)
    spacer_len = len(spacer)

    for i in range(spacer_len ):
        spacer_rev += pair[spacer[-i-1]]

    is_found = False
    while not is_found:
        for shift in range(0,target_seq_len - (spacer_len + 1),1):
            #print(target_seq[shift:shift+spacer_len])
            current_shift = target_seq[shift:shift+spacer_len]
            if spacer == current_shift:
                print('___FOUND___')
                is_found = True
                coords = (shift,shift+spacer_len)
            elif spacer_rev == current_shift:
                print('REV___FOUND___REV')
                is_found = True
                coords = (shift,shift+spacer_len)

    near_spacer_end = target_seq[coords[1]-13:coords[1]+7]
    print('Поиск сайтов рестрикции в:', near_spacer_end)
    return near_spacer_end

def _gc_structure(sequence):
    print('start computing gc....')
    gc_c = sequence.count('G') + sequence.count('C')
    gc_f = int(round(gc_c / len(sequence)*100,0))
    gc_c_10_20 = sequence[9:20].count('G') + sequence[9:20].count('C')
    last_four_nuc = sequence[16:20].count('A') + sequence[16:20].count('G')
    return(gc_c_10_20, gc_f,last_four_nuc)


#создание svg картинки вторичной структуры gRNA
def compute_rna_2d_structure(
        sequence: str,
        gRNA_tail_sequence:str,
        output_dir:str,
        index:str,
        ):
    print('start computing 2d structure....')
    full_seq = sequence + gRNA_tail_sequence
    fc = RNA.fold_compound(full_seq)
    dot_bracked_structure, mfe = fc.mfe()

    output_dir=Path(output_dir)
    output_dir.mkdir(parents = True,exist_ok = True)
    png_path = str(output_dir / f'{index+1}.png')
    
    plot_rna(full_seq,dot_bracked_structure,png_path)

    return (png_path, dot_bracked_structure, float(mfe))


class AnalyzerMachine(QObject):
    error_occurred = Signal(str)
    status_message = Signal(str)
    result_ready = Signal(list)
    progress = Signal(int)
    def __init__(
        self,
        crispor_path: str,
        grna_tail: str,
        output_dir: str,
        configs: dict | None = None,
        doench_border: int = 50,
    ):
        super().__init__()
        self.crispor_path = crispor_path
        self.grna_tail = grna_tail
        self.output_dir = output_dir
        configs = configs or GRNA_DEFAULT_CONFIGS.copy()
        self.doench_threshold = doench_border
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def try_run(self):
        print('trying run...')
        try:
            self._run_analysis()
        except Exception as e:
            self.error_occurred.emit(str(e))

    def _run_analysis(self):
        print('start analysis....')
        if not self.crispor_path:
            self.error_occurred.emit("Путь к файлу не найден")
            return
        if not Path(self.crispor_path).exists():
            self.error_occurred.emit('Ошибка, Выбранный файл не существует')
            return

        crispor_data_frame = pd.read_excel(self.crispor_path,skiprows=8)
        all_target_seq = pd.read_excel(self.crispor_path).iloc[0,1]
        enzymes = get_default_restriction_db()

        grna_count = len(crispor_data_frame)
        grna_result_list = []

        for i, row in crispor_data_frame.iterrows():
            if self._cancelled:
                    return
            print(f'analysis row {i}...')
            target_seq = str(row['targetSeq'])
            print('finding spacer')
            spacer = target_seq[:20]
            print('finding pam')
            pam = target_seq[20:24]
            guide_id_plus_strand = str(row.get('#guideId',''))
            guide_id = int(re.search(r"\d+", guide_id_plus_strand).group())
            strand = re.search(r"[A-Za-z]+$", guide_id_plus_strand).group()
            print(row["Doench '16-Score"])

            if row["Doench '16-Score"] == 'NotEnoughFlankSeq':
                doench16 = 0
            else:
                doench16 = int(row["Doench '16-Score"])
            print(doench16)


            print('trying compute 2d')
            self.status_message.emit(f'Структура {i + 1}/{grna_count}: {spacer}')

            png_path, dot_bracket, mfe = compute_rna_2d_structure(
                sequence=spacer,
                gRNA_tail_sequence=self.grna_tail,
                output_dir=self.output_dir,
                index=i,
            )
            # Анализ
            print('trying analysis grna')
            grna_result = grna_analyzer(
                spacer=spacer,
                pam=pam,
                guide_id=guide_id,
                structure_dot=dot_bracket,
                configs= GRNA_DEFAULT_CONFIGS,
                mfe=mfe,
                png_path = png_path,
                strand=strand,
                doench16=doench16,
                enzymes=enzymes,
                all_target_seq=all_target_seq
            )
            print('trying appending...')
            grna_result_list.append(grna_result)
            progress_val = int((i + 1) / grna_count * 100)
            self.progress.emit(progress_val)
        print('done computing')
        self.result_ready.emit(grna_result_list)

