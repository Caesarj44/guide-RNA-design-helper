from pathlib import Path

from PySide6.QtCore import (
    Qt,
    Signal,
)
from PySide6.QtWidgets import(
    QWidget,
    QVBoxLayout,
    QApplication,
    QFormLayout,
    QCheckBox,
    QGroupBox
)


from core.config_manager import load_configs,load_style,update_config
from gui.tail_selector_combo_box import TailComboBox

class SettingsPage(QWidget):
    configs_changed = Signal(dict)

    def __init__(self, parent = None):
        super().__init__(parent)
        self._build_gui()
    
    def _build_gui(self):
        root_layout = QVBoxLayout(self)
        root_layout.setContentsMargins(16, 16, 16, 16)
        root_layout.setSpacing(12)
        
        # Создание подокошка с настройками интерфейса
        group_gui = QGroupBox('Интерфейс')
        form_gui = QFormLayout(group_gui)
        self.dark_mode = QCheckBox()
        self.dark_mode.stateChanged.connect(self._toggle_theme)
        form_gui.addRow('Тёмная тема:', self.dark_mode)

        # Создание подокошка с настройками CRISPOR анализатора
        # В том числе и выпадающее окошко с кастомными gRNA
        group_grna  = QGroupBox('Настройки gRNA')
        form_grna = QFormLayout(group_grna)
        custom_tails = load_configs('custom_gRNA_tails')
        self.grna_tail_list = TailComboBox()
        self.grna_tail_list.main_list(custom_tails)
        form_grna.addRow('Хвост gRNA:', self.grna_tail_list)

        root_layout.addWidget(group_gui)
        root_layout.addWidget(group_grna)
        root_layout.addStretch()
        self.setLayout(root_layout)
        self._load_settings()
    
    def _load_settings(self):
        cfg = load_configs()
        main = cfg.get('main_configs', {})

        # Тема
        self.dark_mode.blockSignals(True)
        self.dark_mode.setChecked(main.get('dark_theme', False))
        self.dark_mode.blockSignals(False)
        self._apply_theme(main.get('dark_theme', False))

    def _toggle_theme(self, state):
        is_dark = (state == Qt.Checked.value)
        self._apply_theme(is_dark)
        update_config('main_configs', 'dark_theme', is_dark)
        self._emit_configs()

    def _apply_theme(self, is_dark: bool):
        assets = Path(__file__).parent.parent / 'assets'
        qss_file = 'dark.qss' if is_dark else 'light.qss'
        style = load_style(str(assets / qss_file))
        QApplication.instance().setStyleSheet(style)

    def _emit_configs(self):
        cfg = load_configs()
        self.configs_changed.emit(cfg.get('main_configs', {}))



