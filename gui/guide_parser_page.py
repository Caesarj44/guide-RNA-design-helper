from pathlib import Path
import os

from PySide6.QtWidgets import (
    QWidget,
    QMessageBox,
    QVBoxLayout,
    QLineEdit,
    QLabel,
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QFileDialog,
    QProgressBar,
    QRadioButton,
)

from PySide6.QtCore import(
    Qt,
    QThread,
    QRegularExpression
)
from PySide6.QtGui import (
    QRegularExpressionValidator
)

import pandas as pd
import requests
from core.guide_in_sequence_parser import find_grna_in_sequence
from core.config_manager import load_configs
from core.crispor_analyser import AnalyzerMachine
from gui.grna_table import GrnaTableWidget
from gui.guide_map import GuideMapView

class GuideParserPage(QWidget):
    
    def __init__(self):
        super().__init__()
        self.results : list = []
        self._machinery: AnalyzerMachine | None = None
        self._thread: QThread | None = None
        self._all_guides: list = []
        self._pre_pair_guides: list | None = None
        self._pair_mode_active: bool = False
        self._syncing: bool = False
        self._target_seq: str = ""
        self._cut_distance: int = 3
        self._build_gui()

    def _build_gui(self):
        self.what_is_chosen = ''


        self.root_layout = QVBoxLayout()
        self.root_layout.setAlignment(Qt.AlignTop)
        self.parsing_settings = QGroupBox('Настройки поиска направляющих РНК')
        self.parsing_settings_layout = QVBoxLayout(self.parsing_settings)

        self.row1 = QHBoxLayout(self.parsing_settings)
        self.session_name_editline = QLineEdit()
        self.session_name_editline.setPlaceholderText('... --- ...')

        self.fasta_style_button = QRadioButton('.fasta файл')
        self.string_style_button = QRadioButton('Своя последовательность')
        self.ncbi_locus_style_button = QRadioButton('NCBI локус')
        self.crispor_style_button = QRadioButton('CRISPOR')

        self.fasta_style_button.toggled.connect(self._choose_parsing_type)
        self.string_style_button.toggled.connect(self._choose_parsing_type)
        self.ncbi_locus_style_button.toggled.connect(self._choose_parsing_type)
        self.crispor_style_button.toggled.connect(self._choose_parsing_type)

        self.start_parsing_button = QPushButton()
        self.start_parsing_button.setText('▶  Анализ')
        self.start_parsing_button.clicked.connect(self.start_parsing)


        self.row1.addWidget(QLabel('Название сессии:'))
        self.row1.addWidget(self.session_name_editline)
        self.row1.addWidget(self.fasta_style_button)
        self.row1.addWidget(self.string_style_button)
        self.row1.addWidget(self.ncbi_locus_style_button)
        self.row1.addWidget(self.crispor_style_button)
        self.row1.addStretch()
        self.row1.addWidget(self.start_parsing_button)
        
        self.row2 =  QHBoxLayout()
        self.fasta_file_path = QLineEdit()
        self.fasta_file_path.setPlaceholderText('путь до .fasta файла')
        self.fasta_file_path.setVisible(True)

        self.fasta_finder_button = QPushButton()
        self.fasta_finder_button.setText('Найти файл')
        self.fasta_finder_button.clicked.connect(self._browse_fasta_file)
        self.fasta_finder_button.setVisible(True)

        self.custom_sequence_input = QLineEdit()
        self.custom_sequence_input.setVisible(False)
        self.custom_sequence_input.setPlaceholderText('Введи свою последовательность. Только ATGC. Не более 5kbp.')
        self.custom_sequence_input.setValidator(
            QRegularExpressionValidator(QRegularExpression('^[ATGCatgc]*$'))
        )

        self.ncbi_locus_id_input = QLineEdit()
        self.ncbi_locus_id_input.setVisible(False)
        self.ncbi_locus_id_input.setPlaceholderText('id локуса...')
        self.ncbi_locus_id_input.setMaximumWidth(110)

        self.crispor_file_path = QLineEdit()
        self.crispor_file_path.setPlaceholderText('путь до выдачи CRISPOR')
        self.crispor_file_path.setVisible(False)

        self.crispor_finder_button = QPushButton()
        self.crispor_finder_button.setText('Найти файл')
        self.crispor_finder_button.setVisible(False)
        self.crispor_finder_button.clicked.connect(self._browse_crispor_file)
        
        self.coords_text = QLabel('Координаты (необязательно):')
        self.left_coord = QLineEdit()
        self.left_coord.setPlaceholderText('левый')
        self.left_coord.setVisible(True)
        self.left_coord.setMaximumWidth(80)
        self.right_coord = QLineEdit()
        self.right_coord.setPlaceholderText('правый')
        self.right_coord.setVisible(True)
        self.right_coord.setMaximumWidth(80)

        self.row2.addWidget(self.fasta_file_path)
        self.row2.addWidget(self.fasta_finder_button)
        self.row2.addWidget(self.custom_sequence_input)
        self.row2.addWidget(self.ncbi_locus_id_input)
        self.row2.addWidget(self.coords_text)
        self.row2.addWidget(self.left_coord)
        self.row2.addWidget(self.right_coord)
        self.row2.addWidget(self.crispor_file_path)
        self.row2.addWidget(self.crispor_finder_button)
        self.row2.addStretch()

        self.row3 = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(False)

        self.status_label = QLabel('')
        self.status_label.setObjectName('status_label')
        self.status_label.setVisible(False)

        self.text_progress = QLabel('Прогресс:')
        self.text_progress.setVisible(False)

        self.row3.addStretch()
        self.row3.addWidget(self.text_progress)
        self.row3.addWidget(self.progress_bar, 3)
        self.row3.addWidget(self.status_label, 2)
        self.row3.addStretch()

        self.genome_browser = GuideMapView()
        self.genome_browser.setVisible(False)
        self.genome_browser.guide_clicked.connect(self._on_guide_clicked)
        self.genome_browser.setMaximumHeight(100)

        self.grna_table = GrnaTableWidget()
        self.grna_table.setVisible(False)
        self.grna_table.guide_selected.connect(self._on_table_guide_selected)
        self.grna_table.guides_deleted.connect(self._on_guides_deleted)
        self.grna_table.pair_mode_toggled.connect(self._on_pair_mode_toggled)
        self.grna_table.pair_filter_reset.connect(self._on_pair_reset)

        self.parsing_settings_layout.addLayout(self.row1)
        self.parsing_settings_layout.addLayout(self.row2)
        self.parsing_settings_layout.addLayout(self.row3)
        self.root_layout.addWidget(self.parsing_settings)
        self.root_layout.addWidget(self.genome_browser)
        self.root_layout.addWidget(self.grna_table)
        self.setLayout(self.root_layout)

        self.fasta_style_button.setChecked(True)


    def _choose_parsing_type(self, checked):
        if self.fasta_style_button.isChecked():
            self.fasta_file_path.setVisible(True)
            self.fasta_finder_button.setVisible(True)

            self.custom_sequence_input.setVisible(False)

            self.ncbi_locus_id_input.setVisible(False)

            self.crispor_finder_button.setVisible(False)
            self.crispor_file_path.setVisible(False)

            self.coords_text.setVisible(True)
            self.left_coord.setVisible(True)
            self.right_coord.setVisible(True)
            self.what_is_chosen = 'fasta_file'
            self.chosen_type = 'Custom_searcher'

        elif self.string_style_button.isChecked():
            self.fasta_file_path.setVisible(False)
            self.fasta_finder_button.setVisible(False)

            self.custom_sequence_input.setVisible(True)
            
            self.ncbi_locus_id_input.setVisible(False)

            self.crispor_finder_button.setVisible(False)
            self.crispor_file_path.setVisible(False)

            self.coords_text.setVisible(False)
            self.left_coord.setVisible(False)
            self.right_coord.setVisible(False)
            self.what_is_chosen = 'custom_string'
            self.chosen_type = 'Custom_searcher'
        
        elif self.ncbi_locus_style_button.isChecked():
            self.fasta_file_path.setVisible(False)
            self.fasta_finder_button.setVisible(False)

            self.custom_sequence_input.setVisible(False)
            
            self.ncbi_locus_id_input.setVisible(True)

            self.crispor_finder_button.setVisible(False)
            self.crispor_file_path.setVisible(False)

            self.coords_text.setVisible(True)
            self.left_coord.setVisible(True)
            self.right_coord.setVisible(True)
            self.what_is_chosen = 'NCBI'
            self.chosen_type = 'Custom_searcher'
        
        elif self.crispor_style_button.isChecked():
            self.fasta_file_path.setVisible(False)
            self.fasta_finder_button.setVisible(False)

            self.custom_sequence_input.setVisible(False)
            
            self.ncbi_locus_id_input.setVisible(False)

            self.crispor_finder_button.setVisible(True)
            self.crispor_file_path.setVisible(True)

            self.coords_text.setVisible(False)
            self.left_coord.setVisible(False)
            self.right_coord.setVisible(False)
            self.what_is_chosen = 'CRISPOR'
            self.chosen_type = 'CRISPOR_analyzer'            

    def start_parsing(self):
        print('parsing...')
        print(self.what_is_chosen)
        self.grna_table.setVisible(False)
        self.genome_browser.setVisible(False)
        #to remove
        #self.test_target_seq = 'TGGCAGGATATATTGTGGTGTAAACAAATTGACGCTTAGACAACTTAATAACACATTGCGGACGTTTTTAATGTACTGAATTAACGCCGAATTAATTCGGGGGATCTGGATTTTAGTACTGGATTTTGGTTTTAGGAATTAGAAATTTTATTGATAGAAGTATTTTACAAATACAAATACATACTAAGGGTTTCTTATATGCTCAACACATGAGCGAAACCCTATAGGAACCCTAATTCCCTTATCTGGGAACTACTCACACATTATTATGGAGAAACTCGAAGATCCGTCGAGCTTGTCGATCGACAGATCCGGTCGGCATCATAACTTCGTATAGCATACATTATACGAAGTTAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGCTAGAGCAGCTTGAGCTTGGATCAGATTGTCGTTTCCCGCCTTCAGTTTAAACTATCAGTGTTTGACAGGATATATTGGCGGGTAAACCTAAGAGAAAAGAGCGTTTA'

        if self.what_is_chosen == 'custom_string':

            #if self.custom_sequence_input.text():
            self.target_seq = self.custom_sequence_input.text()
            #else:
                #self.target_seq = self.test_target_seq

        elif self.what_is_chosen == 'fasta_file':
            if not self.fasta_file_path.text():
                QMessageBox.warning(self, 'Ошибка', 'Выберите сначала файл')
                return
            if not Path(self.fasta_file_path.text()).exists():
                QMessageBox.warning(self, 'Ошибка', 'Выбранный файл не существует')
                return
            
            self.target_seq, fasta_header = self._read_fasta_file(path=self.fasta_file_path.text())


        elif self.what_is_chosen == 'NCBI':
            if not self.ncbi_locus_id_input.text():
                QMessageBox.warning(self, 'Ошибка', 'Введите ID')
                return
            
            self.target_seq = self._get_sequence_from_ncbi(self.ncbi_locus_id_input.text())

            if self.target_seq is None:
                QMessageBox.warning(self, 'Ошибка', 'Ошибка при поиске в NCBI')
                return            

        elif self.what_is_chosen == 'CRISPOR':
            self.path_crispor = self.crispor_file_path.text()
            if not self.path_crispor:
                QMessageBox.warning(self, 'Ошибка', 'Выберите сначала файл')
                return
            if not Path(self.path_crispor).exists():
                QMessageBox.warning(self, 'Ошибка', 'Выбранный файл не существует')
                return
            self.target_seq = pd.read_excel(self.crispor_file_path.text()).iloc[0,1]

        if len(self.target_seq) < 50:
            QMessageBox.warning(self,'Ошибка','Последовательность слишком короткая')
            return

        left_coord_value = 0
        right_coord_value = len(self.target_seq)
        if self.left_coord.text() or self.right_coord.text():

            if self.left_coord.text():
                left_coord_value = int(self.left_coord.text())
            if self.right_coord.text():
                if int(self.right_coord.text()) > len(self.target_seq):
                    QMessageBox.warning(self,'Ошибка','Правая координата больше длины последовательности')
                    return
                right_coord_value = int(self.right_coord.text())

            if (right_coord_value < left_coord_value + 50):
                QMessageBox.warning(self,'Ошибка','Установите правильные координаты')
                return
            
            self.target_seq = self.target_seq[left_coord_value:right_coord_value]

        if len(self.target_seq) > 5000:
            QMessageBox.warning(self, 'Ошибка', 'Длина последовательности больше 5000 п.н.')
            return

        self.session_name = ''
        session_name_text = self.session_name_editline.text().strip()

        if session_name_text:
            self.session_name = session_name_text

        elif self.what_is_chosen == 'CRISPOR':
            self.session_name = os.path.splitext(os.path.basename(self.path_crispor))[0]
        
        elif self.what_is_chosen == 'fasta_file':
            self.session_name = fasta_header
        
        elif self.what_is_chosen == 'NCBI':
            self.session_name = str(self.ncbi_locus_id_input.text()) + '_guides'

        else:
            for i in range(100):
                candidate = Path(__file__).parent.parent / f'output/gRNA result {i+1}'
                if not candidate.exists():
                    self.session_name = f'gRNA result {i+1}'
                    break
            else:
                self.session_name = 'gRNA result out of range'

        print(f'test{self.session_name}test')
        print(f'session name {self.session_name}')

        self.output_dir = Path(__file__).parent.parent / f'output/{self.session_name}'
        print(f'output dir {self.output_dir}')

        if load_configs('preset_gRNA_tails', load_configs('main_configs', 'selected_tail')) is not None:
            self.selected_cas_system = load_configs('preset_gRNA_tails', load_configs('main_configs', 'selected_tail'))

        elif load_configs('custom_gRNA_tails', load_configs('main_configs', 'selected_tail')) is not None:
            self.selected_cas_system = load_configs('custom_gRNA_tails', load_configs('main_configs', 'selected_tail'))['tail']

        self.tail = self.selected_cas_system['tail']
        self.spacer_len = self.selected_cas_system['spacer_len']
        self.pam_seq = self.selected_cas_system['pam']


        if self.chosen_type == 'CRISPOR_analyzer':
            self.input_grna_data_frame =pd.read_excel(self.crispor_file_path.text(),skiprows=8)
        elif self.chosen_type == 'Custom_searcher':
            self.input_grna_data_frame = find_grna_in_sequence(pam_seq=self.pam_seq,spacer_len=self.spacer_len, target_seq=self.target_seq)



        self.grna_table.clear()
        self.results.clear()
        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self.text_progress.setVisible(True)
        self.status_label.setVisible(True)
        self.start_parsing_button.setEnabled(False)
        self.grna_table.set_controls_enabled(False)


        if self._thread is not None and self._thread.isRunning():
            if self._machinery is not None:
                self._machinery.cancel()
            self._thread.quit()
            self._thread.wait(30)
        if self._thread is not None:
            self._thread.deleteLater()
            self._thread = None
        if self._machinery is not None:
            self._machinery.deleteLater()
            self._machinery = None

        self._thread = QThread()
        self._machinery = AnalyzerMachine(
            grna_tail=self.tail,
            output_dir=self.output_dir,
            configs=None,
            doench_border=50,
            type_of_machine=self.chosen_type,
            input_grna_data_frame=self.input_grna_data_frame,
            target_seq=self.target_seq
        )

        self._machinery.moveToThread(self._thread)
        self._thread.started.connect(self._machinery.try_run)
        self._machinery.progress.connect(self.progress_bar.setValue)
        self._machinery.status_message.connect(self.status_label.setText)
        self._machinery.result_ready.connect(self._analysis_done)
        self._thread.start()

    def _analysis_done(self, results):
        print('done parsing')
        self._all_guides = results
        self._pre_pair_guides = None
        self._pair_mode_active = False
        self._target_seq = self.target_seq
        self.grna_table.pair_mode_button.setChecked(False)
        self.grna_table.pair_reset_button.setVisible(False)
        self._refresh_views()
        self.grna_table.setVisible(True)
        self.progress_bar.setVisible(False)
        self.status_label.setVisible(False)
        self.text_progress.setVisible(False)
        self.grna_table.set_controls_enabled(True)
        self.start_parsing_button.setEnabled(True)
        self.genome_browser.setVisible(True)

    def _refresh_views(self):
        self.grna_table.populate(self._all_guides)
        self.genome_browser.set_guides(
            self._all_guides,
            region_start=0,
            region_end=len(self._target_seq),
            target_seq=self._target_seq,
            cut_distance=self._cut_distance,
        )

    def _cut_site(self, guide):
        cut_dist = self._cut_distance
        if guide.strand == '+':
            return guide.end - cut_dist
        else:
            return guide.start + cut_dist

    def _on_guide_clicked(self, idx):
        if self._pair_mode_active:
            self._apply_pair_filter(idx)
        else:
            self._syncing = True
            self.genome_browser.highlight_guide(idx)
            self.grna_table.select_guide_by_index(idx)
            self._syncing = False
        self.grna_table.pair_mode_button.setChecked(False)

    def _on_table_guide_selected(self, idx):
        if self._syncing:
            return
        self.genome_browser.highlight_guide(idx)

    def _on_guides_deleted(self, remaining_indices):
        remaining_set = set(remaining_indices)
        self._all_guides = [guide for guide in self._all_guides if guide.indx in remaining_set]
        self._refresh_views()

    def _on_pair_mode_toggled(self, active):
        self._pair_mode_active = active
        # if not active:
        #     self.grna_table.pair_reset_button.setVisible(False)

    def _on_pair_reset(self):
        if self._pre_pair_guides is not None:
            self._all_guides = list(self._pre_pair_guides)
            self._pre_pair_guides = None
            self._refresh_views()
        self.grna_table.pair_reset_button.setVisible(False)

    def _apply_pair_filter(self, selected_idx: int):
        selected = next((guide for guide in self._all_guides if guide.indx == selected_idx), None)
        if selected is None:
            return
        
        target_len = len(self._target_seq) if self._target_seq else 0
        selected_cut = self._cut_site(selected)
        if self._pre_pair_guides is None:
            self._pre_pair_guides = list(self._all_guides)

        filtered = []
        for guide in self._all_guides:
            if guide.indx == selected_idx:
                filtered.append(guide)
                continue

            if guide.start < selected.start:
                continue

            dist = abs(self._cut_site(guide) - selected_cut)
            if dist < 50:
                continue
            if dist > 2000:
                continue

            filtered.append(guide)

        self._all_guides = filtered
        self._refresh_views()
        self.grna_table.pair_reset_button.setVisible(True)

    def _browse_crispor_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, 'Выбери файл CRISPOR', '',
            'Excel файлы (*.xlsx *.xls);;Все файлы (*)',
        )
        if path:
            self.crispor_file_path.setText(path)

    def _browse_fasta_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, 'Выбери файл CRISPOR', '',
            'Excel файлы (*.fasta *.fa);;Все файлы (*)',
        )
        if path:
            self.fasta_file_path.setText(path)
    
    def _read_fasta_file(self,path):
        header = ''
        sequence = ''
        with open (path,'r') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    header = line[1:] 
                elif not line.startswith('>'):
                    sequence = line
                            
        return sequence, header
    
    def _get_sequence_from_ncbi(self, id, database = 'nucleotide',rettype ='fasta'):
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'id' : id,
            'db' : database,
            'rettype' : rettype,
            'retmode' : 'text'
        }
        
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            if rettype == "fasta":
                lines = response.text.strip().split('\n')
                #header = lines[0]
                sequence = ''.join(lines[1:])
                print(sequence)
                return sequence
            else:
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Ошибка запроса: {e}")
            return None