import sys
import json
import os
import csv
import pandas as pd
import time
import cairosvg
from dataclasses import dataclass, field
from pathlib import Path
from PySide6.QtCore import Qt, QSize, QSettings
from PySide6.QtWidgets import (QApplication, QLabel, QWidget, QListWidget, QListWidgetItem,
                               QPushButton, QVBoxLayout, QHBoxLayout,QStackedWidget,QFormLayout,
                               QCheckBox,QLineEdit,QComboBox,QSpinBox, QFileDialog)
from PySide6.QtGui import QFontDatabase, QFont,QIcon    

@dataclass
class gRNA_sequence():
    sequence : str 
    PAM : str
    oligoT : int = 0
    strain : int = 0
    GC_freq : int = 0 
    GC_count_10to20 : int = 0 
    last_four_pur : int = 0 
    access_18_to_20 : int = 0
    seven_links : int = 0
    twelve_links : int = 0
    C_not_at_3 : int = 0 
    G_not_at_16 : int = 0 
    C_at_16 : int = 0 
    G_or_A_at_20 : int = 0 
    C_at_18 : int = 0 
    G_not_at_14 : int = 0
    GCC_not_at_16to20 : int = 0
    PAMs_N : int = 0
    restriction_sites : list = ()
    total_score : int = 0

DEFAULT_CONFIGS = {'main_configs': {
    "dark_theme": False,
    "default_gRNA_tail" :'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTT',

}
}

CONFIGS_FILE = Path(__file__).parent / "configs.json"

def load_configs(section=None, key=None):
    if not CONFIGS_FILE.exists():
        save_configs(DEFAULT_CONFIGS)
        data = DEFAULT_CONFIGS.copy()
    else:
        try:
            with open(CONFIGS_FILE, "r", encoding="utf-8") as f:
                data = json.load(f)
        except json.JSONDecodeError:
            data = DEFAULT_CONFIGS.copy()

    for sec, values in DEFAULT_CONFIGS.items():
        if sec not in data:
            data[sec] = values
        elif isinstance(values, dict):
            for k, v in values.items():
                if k not in data[sec]:
                    data[sec][k] = v

    if section is None:
        return data
    if key is None:
        return data.get(section, {})
    return data.get(section, {}).get(key)

def save_configs(data: dict):
    with open(CONFIGS_FILE, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def update_config(section: str, key: str, value):
    current_CONFIGS = load_configs()
    if section not in current_CONFIGS:
        current_CONFIGS[section] = {}
    current_CONFIGS[section][key] = value
    save_configs(current_CONFIGS)

def load_style(filename):
    with open("assets/" + filename, "r", encoding="utf-8") as f:
        return f.read()

#==== страница основная
class MainPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Главная страница"))
        self.setLayout(layout)

#==== страница с настройками
class SettingsPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()

        # Чекбокс тёмной темы
        self.dark_mode = QCheckBox()
        self.dark_mode.stateChanged.connect(self._toggle_theme)
        layout.addRow("Тёмная тема:", self.dark_mode)

        # Поле хвоста gRNA
        self.gRNA_tail_input = QLineEdit()
        self.gRNA_tail_input.editingFinished.connect(self._save_tail)
        layout.addRow("Хвост gRNA по умолчанию:", self.gRNA_tail_input)

        self.setLayout(layout)
        self._load_settings()

    def _load_settings(self):
        self.dark_mode.blockSignals(True)
        self.dark_mode.setChecked(load_configs("main_configs", "dark_theme"))
        self.dark_mode.blockSignals(False)

        self.gRNA_tail_input.setText(
            load_configs("main_configs", "default_gRNA_tail")
        )

        self._apply_theme(load_configs("main_configs", "dark_theme"))

    def _toggle_theme(self, state):
        is_dark = (state == Qt.Checked.value)
        self._apply_theme(is_dark)
        update_config("main_configs", "dark_theme", is_dark)

    def _apply_theme(self, is_dark: bool):
        app = QApplication.instance()
        app.setStyleSheet(load_style("dark.qss") if is_dark else load_style("light.qss"))

    def _save_tail(self):
        tail = self.gRNA_tail_input.text()
        update_config("main_configs", "default_gRNA_tail", tail) 

#==== страница скоринга CRISPOR
class CRISPORPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()
        self.start_analysis = QPushButton('Начать анализ')
        self.find_file = QPushButton('Найти файл')
        self.find_file.clicked.connect(self._open_dialog)
        self.CRISPOR_input = QLineEdit()
        self.CRISPOR_input.setPlaceholderText('Путь до файла')
        self.gRNA_tail_input = QLineEdit()
        self.gRNA_tail_input.setPlaceholderText('Используется хвост по умолчанию')

        gRNA_tail_row = QHBoxLayout()
        gRNA_tail_row.addWidget(QLabel('Хвост gRNA:'),1)
        gRNA_tail_row.addWidget(self.gRNA_tail_input,10)

        find_file_row = QHBoxLayout()
        find_file_row.addWidget(self.CRISPOR_input,4)
        find_file_row.addWidget(self.find_file,1)

        CRISPOR_config_row = QHBoxLayout()
        CRISPOR_config_row.addWidget(QLabel('Настройки CRISPOR-анализа'))

        layout.addRow(find_file_row)
        layout.addRow(gRNA_tail_row)
        layout.addRow(CRISPOR_config_row)
        layout.addRow(self.start_analysis)
        self.setLayout(layout)

    def _open_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Выбери файл",
            "",
            "EXCEL файлы (*.xlsx *.xlx);;Все файлы (*)"
        )
        if path:
            self.CRISPOR_input.setText(path)

    def _run_analysis(self):
        path = self.CRISPOR_input.text()
        if not path:
            return



class Widget(QWidget):
    def __init__(self, parent=None):
        super(Widget, self).__init__(parent)

        self.menu_widget = QListWidget()
        self.menu_widget.setUniformItemSizes(True)
        self.menu_widget.setFixedWidth(50)
        self.menu_widget.setIconSize(QSize(40,40))
        self.menu_widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff) # отключение бесящего скролбара

#==== добавление пункта меню с домашней страницей
        home_item = QListWidgetItem()
        home_item.setSizeHint(QSize(48, 48))
        home_item.setIcon(QIcon('assets/home48x48.ico'))
        self.menu_widget.addItem(home_item) 

#==== добавление пункта меню с настройками
        settings_item = QListWidgetItem()
        settings_item.setTextAlignment(Qt.AlignCenter)
        settings_item.setSizeHint(QSize(48,48))
        settings_item.setIcon(QIcon('assets/settings48x48.ico'))
        self.menu_widget.addItem(settings_item)
#==== добавление пункта меню с данными из CRISPOR
        CRISPOR_item = QListWidgetItem()
        CRISPOR_item.setTextAlignment(Qt.AlignCenter)
        CRISPOR_item.setSizeHint(QSize(48,48))
        CRISPOR_item.setText('CRISPOR')
        self.menu_widget.addItem(CRISPOR_item)
#==== добавление пункта меню с собственным анализом gRNA

#==== создание стека страниц для переключения боковой панелью
        self.stack = QStackedWidget()
        self.stack.addWidget(MainPage())      # индекс 0
        self.stack.addWidget(SettingsPage())  # индекс 1
        self.stack.addWidget(CRISPORPage())  # индекс 2
        self.menu_widget.currentRowChanged.connect(self.stack.setCurrentIndex)

        text_widget = QLabel("lorem ipsum")

        content_layout = QVBoxLayout()
        content_layout.addWidget(text_widget)

        self.main_widget = QWidget()
        self.main_widget.setLayout(content_layout)

        self.main_layout = QHBoxLayout()
        self.main_layout.addWidget(self.menu_widget)
        self.main_layout.addWidget(self.stack)
        self.setLayout(self.main_layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)

#==== загрузка сохраненных настроек
    cfg = load_configs()
    if cfg['main_configs']["dark_theme"]:
        app.setStyleSheet(load_style("dark.qss"))
    else:
        app.setStyleSheet(load_style("light.qss"))

    w = Widget()
    w.resize(800, 600)
    w.show()

    sys.exit(app.exec())