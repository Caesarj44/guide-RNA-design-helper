
from PySide6.QtWidgets import (
    QComboBox,
    QMessageBox,
    QVBoxLayout,
    QDialog,
    QFormLayout,
    QLineEdit,
    QDialogButtonBox,
    QLabel,
)

from PySide6.QtCore import(
    QRegularExpression,
    Signal,
)
from PySide6.QtGui import (
    QRegularExpressionValidator,
)
from core.config_manager import load_configs, update_config

PRESET_TAILS = load_configs('preset_gRNA_tails')
_ADD_CUSTOM_MARKER = '__add_custom__'

class TailComboBox(QComboBox):
    
    def __init__(self):
        super().__init__()
        self.setEditable(False)
        self.setInsertPolicy(QComboBox.NoInsert)
        self._custom_names = set()

        self.currentIndexChanged.connect(self._index_changed)

    def main_list(self, custom_tails: dict):
        self.blockSignals(True)
        self.clear()
        self._custom_names.clear()

        # Пресеты
        for name,keys  in PRESET_TAILS.items():
            self.addItem(f'{keys['spacer_len']}bp, {keys['pam']}, {name}')
        # Разделитель
        self.insertSeparator(self.count())
        # Кастомные 
        for name, seq in custom_tails.items():
            self.addItem(name, seq)
            self._custom_names.add(name)
        # Разделитель
        self.insertSeparator(self.count())
        self.addItem('✚  Добавить свой хвост...', _ADD_CUSTOM_MARKER)
        last_selected_tail = load_configs('main_configs', 'selected_tail')
        if last_selected_tail:
            idx = self.findText(last_selected_tail)
            if idx >= 0:
                self.setCurrentIndex(idx)
        self.blockSignals(False)
        
    
    def _index_changed(self, idx):
        if idx < 0:
            return
        data = self.itemData(idx)
        if data == _ADD_CUSTOM_MARKER:
            self.setCurrentIndex(0)
            self._add_custom_tail()
            return
        update_config('main_configs','selected_tail',self.itemText(idx))
        print(load_configs('main_configs','selected_tail'))

    def _add_custom_tail(self):
        dialog = AddTailDialog(self)
        if dialog.exec() == AddTailDialog.Accepted:
            name = dialog.name_input.text().strip()
            seq = dialog.seq_input.text().strip().upper()

            if not name or not seq:
                return
            if not all(i in 'ATGC' for i in seq):
                QMessageBox.warning(self, 'Ошибка', 'Последовательность должна содержать только ATGC')
                return
            if name in PRESET_TAILS:
                QMessageBox.warning(self, 'Ошибка', f'Название «{name}» уже занято пресетом')
                return

            insert_pos = self.count() - 2
            self.insertItem(insert_pos, name, seq)
            self._custom_names.add(name)
            self.setCurrentIndex(insert_pos)

            update_config('custom_gRNA_tails', name, seq)
            update_config('main_configs','selected_tail',name)


    def get_tail_sequence(self) -> str:
        idx = self.currentIndex()
        if idx >= 0 and self.itemData(idx) is not None and self.itemData(idx) != _ADD_CUSTOM_MARKER:
            return self.itemData(idx)
        return self.currentText().strip().upper()

    def is_custom(self, idx: int) -> bool:
        name = self.itemText(idx)
        return name in self._custom_names

    custom_tail_added = Signal(str, str)  
    custom_tail_removed = Signal(str)

class AddTailDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Добавить хвост gRNA')
        self.setMinimumWidth(400)
        self.setModal(True)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(16)

        # Пояснение
        hint = QLabel('Введите название и последовательность хвоста gRNA.\nДопустимы только символы ATGC')
        hint.setProperty('class', 'secondary')
        layout.addWidget(hint)

        # Поля ввода
        form = QFormLayout()
        form.setSpacing(12)

        self.name_input = QLineEdit()
        self.name_input.setPlaceholderText('Например: My Cas9')
        self.name_input.setMaxLength(50)
        form.addRow('Название:', self.name_input)

        self.seq_input = QLineEdit()
        self.seq_input.setPlaceholderText('ATGCGTACGT...')
        self.seq_input.setValidator(
            QRegularExpressionValidator(QRegularExpression('^[ATGCatgc]*$'))
        )
        form.addRow('Последовательность:', self.seq_input)

        self.pam_input = QLineEdit()
        self.pam_input.setValidator(
            QRegularExpressionValidator(QRegularExpression('^[ATGCRYSWKMBDHVNatgcryswkmbdhvn]*$'))
        )
        self.pam_input.setPlaceholderText('NGG')
        form.addRow('PAM-сайт',self.pam_input)

        self.spacer_len_input = QLineEdit()
        self.spacer_len_input.setPlaceholderText('20')
        form.addRow('Длина спейсера',self.spacer_len_input)

        self.cut_site_input = QLineEdit()
        self.cut_site_input.setPlaceholderText('3')
        self.cut_site_input.setValidator(
            QRegularExpressionValidator(QRegularExpression('^[0-9]$'))
        )
        form.addRow("Сайт разреза от 3' конца", self.cut_site_input)

        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._validate_and_accept)
        buttons.rejected.connect(self.reject)

        ok_btn = buttons.button(QDialogButtonBox.Ok)
        ok_btn.setProperty('class', 'primary')
        ok_btn.setText('Добавить')
        cancel_btn = buttons.button(QDialogButtonBox.Cancel)
        cancel_btn.setText('Отмена')
        layout.addWidget(buttons)

    def _validate_and_accept(self):
        name = self.name_input.text().strip()
        seq = self.seq_input.text().strip().upper()
        pam = self.pam_input.text().strip().upper()
        spacer_len = self.spacer_len_input.text().strip()
        cut_site = self.cut_site_input.text().strip()

        if not name:
            self.name_input.setFocus()
            return
        if not seq:
            self.seq_input.setFocus()
            return
        if not pam:
            self.pam_input.setFocus()
            return
        if not spacer_len:
            self.spacer_len_input.setFocus()
            return        
        if not cut_site:
            self.cut_site_input.setFocus()
            return       

        self.accept()

