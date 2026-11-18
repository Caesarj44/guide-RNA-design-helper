from pathlib import Path
import os
from datetime import datetime

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
    QTextEdit,
    QButtonGroup,
    QGridLayout
)

from PySide6.QtCore import(
    Qt,
    QThread,
    QRegularExpression
)
from PySide6.QtGui import (
    QRegularExpressionValidator,
    QPixmap
)

from ViennaRNA import RNA

from core.RNA_structure_plot import plot_rna

OUTPUT_FOLDER = (
    Path(__file__).resolve().parent.parent
    / "output"
    / "single_RNA_structure"
)

class TDStructureImager (QWidget):

    def __init__(self):
        super().__init__()
        self._structure_pixmap = None
        self._build_gui()

    def _build_gui(self):
        self.root_layout = QVBoxLayout()
        self.root_layout.setAlignment(Qt.AlignTop)
        self.RNA_info_settings_groupbox = QGroupBox('Фолдинг РНК')
        self.RNA_info_settings_layout = QVBoxLayout(self.RNA_info_settings_groupbox)

# first row with two column
        self.input_and_settings_row_layout = QHBoxLayout()
# left part
        self.left_settings_groupbox = QGroupBox('Настройки моделирования')
        self.left_settings_layout = QVBoxLayout(self.left_settings_groupbox)
    # model type row
        self.model_type_row_layout = QHBoxLayout()
        self.model_type_row_layout.addWidget(QLabel("Тип модели построения"))
        self.centroid_model_radiobutton = QRadioButton('Centroid')
        self.MFE_model_radiobutton = QRadioButton('MFE')
        self.MEA_model_radiobutton = QRadioButton('MEA')
        self.model_type_radiobutton_group = QButtonGroup(self)
        self.model_type_radiobutton_group.addButton(self.centroid_model_radiobutton)
        self.model_type_radiobutton_group.addButton(self.MFE_model_radiobutton)
        self.model_type_radiobutton_group.addButton(self.MEA_model_radiobutton)
        self.model_type_row_layout.addWidget(self.centroid_model_radiobutton)
        self.model_type_row_layout.addWidget(self.MFE_model_radiobutton)
        self.model_type_row_layout.addWidget(self.MEA_model_radiobutton)
        # params types rows
        self.params_type_row_layout = QHBoxLayout()
        self.params_type_row_layout.addWidget(QLabel("Тип модели построения"))
        self.parameter_set_dna_mathews1999 = QRadioButton('Mathews1999')
        self.parameter_set_dna_mathews2004 = QRadioButton('Mathews2004')
        self.parameter_set_rna_andronescu2007 = QRadioButton('Andronescu2007')
        self.parameter_set_rna_langdon2018 = QRadioButton('Langdon2018')
        self.parameter_set_rna_turner1999 = QRadioButton('Turner1999')
        self.parameter_set_rna_turner2004 = QRadioButton('Turner2004')
        self.params_group = QButtonGroup(self)
        self.params_group.addButton(self.parameter_set_dna_mathews1999)
        self.params_group.addButton(self.parameter_set_dna_mathews2004)
        self.params_group.addButton(self.parameter_set_rna_andronescu2007)
        self.params_group.addButton(self.parameter_set_rna_langdon2018)
        self.params_group.addButton(self.parameter_set_rna_turner1999)
        self.params_group.addButton(self.parameter_set_rna_turner2004)

        # 2 buttons per row
        self.params_grid_layout = QGridLayout()
        self.params_grid_layout.addWidget(self.parameter_set_rna_turner2004, 0, 0)
        self.params_grid_layout.addWidget(self.parameter_set_rna_turner1999, 0, 1)
        self.params_grid_layout.addWidget(self.parameter_set_rna_andronescu2007, 1, 0)
        self.params_grid_layout.addWidget(self.parameter_set_rna_langdon2018, 1, 1)
        self.params_grid_layout.addWidget(self.parameter_set_dna_mathews1999, 2, 0)
        self.params_grid_layout.addWidget(self.parameter_set_dna_mathews2004, 2, 1)

        # all layouts adding
        self.left_settings_layout.addLayout(self.model_type_row_layout)
        self.left_settings_layout.addLayout(self.params_grid_layout)

# right part
        self.right_info_gRNA_layout = QVBoxLayout()

        self.first_right_row = QHBoxLayout()
        self.start_button = QPushButton()
        self.start_button.setText('▶  Построение')
        self.start_button.setProperty('class','start')
        self.start_button.clicked.connect(self._start_grna_modeling)
        self.first_right_row.addWidget(QLabel('Введите последовательность gRNA в окно ниже. До 2000 н.'))
        self.first_right_row.addStretch()
        self.first_right_row.addWidget(self.start_button)
        self.sequence_input = QLineEdit()
        self.sequence_input.setValidator(QRegularExpressionValidator(QRegularExpression('^[ATGCatgc]*$')))
        self.sequence_input.setMaxLength(2000)

        self.right_info_gRNA_layout.addLayout(self.first_right_row)
        self.right_info_gRNA_layout.addWidget(self.sequence_input)

        self.input_and_settings_row_layout.addWidget(self.left_settings_groupbox)
        self.input_and_settings_row_layout.addLayout(self.right_info_gRNA_layout)

# 2nd row with results
        self.second_row_with_results_groupbox = QGroupBox()
        self.second_row_with_results = QHBoxLayout(self.second_row_with_results_groupbox)
# left part
        self.left_result_groupbox = QGroupBox("Вторичная структура")
        self.left_result_groupbox.setObjectName('title')
        self.left_result_layout = QVBoxLayout(self.left_result_groupbox)
    # image
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setMinimumSize(400, 400)
        self.left_result_layout.addWidget(self.image_label)
        self.sequence_length_label = QLabel("Длина: 0 н.")
        self.left_result_layout.addWidget(self.sequence_length_label)
        

        self.second_row_with_results.addWidget(self.left_result_groupbox)
        
        self.RNA_info_settings_layout.addLayout(self.input_and_settings_row_layout)
        self.RNA_info_settings_layout.addWidget(self.second_row_with_results_groupbox)
        self.root_layout.addWidget(self.RNA_info_settings_groupbox)
        self.setLayout(self.root_layout)
        self.MFE_model_radiobutton.setChecked(True)
        self.parameter_set_rna_turner2004.setChecked(True)
        self.second_row_with_results_groupbox.setVisible(False)
#right part
        self.right_result_groupbox = QGroupBox('Характеристика')
        self.right_result_layout = QVBoxLayout(self.right_result_groupbox)
    #thermodynamics rows
        self.result_label = QLabel()
        self.right_result_layout.addWidget(self.result_label)

        self.second_row_with_results.addWidget(self.right_result_groupbox)

    def _start_grna_modeling(self):
        self.second_row_with_results_groupbox.setVisible(False)
        #rna_sequence = self.sequence_input.text()
        rna_sequence = 'GATCAATCACTTGCCTTCAGGTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT'
        if len(rna_sequence) == 0:
            QMessageBox.warning(
                self,
                'Ошибка',
                f'Введите последовательность'
            )
            return
        self.sequence_length_label.setText(f"Длина: {len(rna_sequence)} н.")

        if self.parameter_set_dna_mathews1999.isChecked():
            RNA.params_load_DNA_Mathews1999()
        elif self.parameter_set_dna_mathews2004.isChecked():
            RNA.params_load_DNA_Mathews2004()
        elif self.parameter_set_rna_andronescu2007.isChecked():
            RNA.params_load_RNA_Andronescu2007()
        elif self.parameter_set_rna_langdon2018.isChecked():
            RNA.params_load_RNA_Langdon2018()
        elif self.parameter_set_rna_turner1999.isChecked():
            RNA.params_load_RNA_Turner1999()
        elif self.parameter_set_rna_turner2004.isChecked():
            RNA.params_load_RNA_Turner2004()
        
        rna_sequence = rna_sequence.replace('T','U')


        fc = RNA.fold_compound(rna_sequence)
        mfe_structure, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        _, ensemble_energy = fc.pf()

        if self.centroid_model_radiobutton.isChecked():
            model_type = 'Centroid'
            dot_bracket_structure, centroid_distance = fc.centroid()
            model_specific = f"Среднее расстояние до ансамбля: {centroid_distance:.2f} нт"

        elif self.MFE_model_radiobutton.isChecked():
            model_type = 'MFE'
            dot_bracket_structure = mfe_structure
            model_specific = ""

        elif self.MEA_model_radiobutton.isChecked():  
            model_type = 'MEA'
            gamma = 1.0
            dot_bracket_structure, mea_score = fc.MEA(gamma)
            model_specific = f"MEA score (γ = {gamma}): {mea_score:.2f}<br>"

        structure_energy = fc.eval_structure(dot_bracket_structure)
        ensemble_defect = fc.ensemble_defect(dot_bracket_structure)
        freq_mfe = fc.pr_structure(mfe_structure)

        self.result_label.setText(
            f"<b>GS: <b> {(rna_sequence.count('G')+rna_sequence.count('C'))/len(rna_sequence)*100:.0f}%<br>"
            f"<b>Модель:</b> {model_type}<br>"
            f"<b>MFE:</b> {mfe:.2f} ккал/моль<br>"
            f"<b>DBS:<b> {dot_bracket_structure}<br>"
            f"<b>Энергия предсказанной структуры:</b> {structure_energy:.2f} ккал/моль<br>"
            f"<b>Свободная энергия ансамбля:</b> {ensemble_energy:.2f} ккал/моль<br>"
            f"<b>Частота MFE в ансамбле:</b> {freq_mfe:.3f}<br>"
            f"<b>Ensemble defect:</b> {ensemble_defect:.3f}<br>"
            + model_specific
        )


        
        time_now = datetime.now().strftime("%H:%M:%S %d,%m,%y")

        OUTPUT_FOLDER.mkdir(parents=True,exist_ok=True)
        output_dir = OUTPUT_FOLDER / f'{time_now}'
        plot_rna(
            sequence=rna_sequence, 
            structure= dot_bracket_structure,
            output_path=output_dir
        )

        self._display_structure_image(str(output_dir)+'.png')
        
        print(rna_sequence)
        print(model_type)
        self.second_row_with_results_groupbox.setVisible(True)


    def _update_image_label(self):
        if self._structure_pixmap is None:
            return
        scaled = self._structure_pixmap.scaled(
            self.image_label.size(),
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation,
        )
        self.image_label.setPixmap(scaled)

    def _display_structure_image(self, image_path):
        pixmap = QPixmap(str(image_path))
        if pixmap.isNull():
            QMessageBox.warning(
                self,
                'Ошибка',
                f'Не удалось загрузить изображение:\n{image_path}'
            )
            return
        self._structure_pixmap = pixmap
        self._update_image_label()
