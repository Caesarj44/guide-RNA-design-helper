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
    QTextEdit
)

from PySide6.QtCore import(
    Qt,
    QThread,
    QRegularExpression
)
from PySide6.QtGui import (
    QRegularExpressionValidator
)

class TDStructureImager (QWidget):

    def __init__(self):
        super().__init__()
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
        self.model_type_row_layout = QHBoxLayout()
    # model type row
        self.model_type_row_layout.addWidget(QLabel("Тип модели построения"))
        self.centroid_model_radiobutton = QRadioButton('Centroid')
        self.MFE_model_radiobutton = QRadioButton('MFE')
        self.model_type_row_layout.addWidget(self.centroid_model_radiobutton)
        self.model_type_row_layout.addWidget(self.MFE_model_radiobutton)

    # all layouts adding
        self.left_settings_layout.addLayout(self.model_type_row_layout)

# right part
        self.right_info_gRNA_layout = QVBoxLayout()

        self.first_right_row = QHBoxLayout()
        self.start_button = QPushButton()
        self.start_button.setText('▶  Построение')
        self.start_button.clicked.connect(self._start_grna_modeling)
        self.first_right_row.addWidget(QLabel('Введите последовательность gRNA в окно ниже. До 2000 п.н.'))
        self.first_right_row.addStretch()
        self.first_right_row.addWidget(self.start_button)
        self.sequence_input = QLineEdit()
        self.sequence_input.setValidator(QRegularExpressionValidator(QRegularExpression('^[ATGCatgc]*$')))
        self.sequence_input.setMaxLength(2000)

        self.right_info_gRNA_layout.addLayout(self.first_right_row)
        self.right_info_gRNA_layout.addWidget(self.sequence_input)


        self.input_and_settings_row_layout.addWidget(self.left_settings_groupbox)
        self.input_and_settings_row_layout.addLayout(self.right_info_gRNA_layout)
        self.RNA_info_settings_layout.addLayout(self.input_and_settings_row_layout)
        self.root_layout.addWidget(self.RNA_info_settings_groupbox)
        self.setLayout(self.root_layout)
        self.MFE_model_radiobutton.setChecked(True)

    def _start_grna_modeling(self):
        # grna_sequence = self.sequence_input.text()
        grna_sequence = 'TGCTCCATACAAGCCAACCAGTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT'
        if self.centroid_model_radiobutton.isChecked:
            model_type = 'Centroid'
        elif self.MFE_model_radiobutton.isChecked:
            model_type = 'MFE'
        print(grna_sequence)
        print(model_type)