from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QFormLayout,
    QLabel,
    QGroupBox

)


class HomePage(QWidget):
    def __init__(self, parent = None):
        super().__init__(parent)
        self._build_gui()

    def _build_gui(self):
        root_layout = QVBoxLayout(self)
        root_layout.setContentsMargins(16, 16, 16, 16)
        root_layout.setSpacing(12)
        
        # Создание подокошка с основной инфой по программе
        group_info = QGroupBox('Основная информация о программе')
        form_info = QFormLayout(group_info)
        self.program_name_and_version = QLabel()
        self.program_name_and_version.setText('GRDH v1.1.1\nguide RNA design helper')
        form_info.addRow(self.program_name_and_version)
        
        group_authors = QGroupBox('Авторы')
        form_authors = QFormLayout(group_authors)
        self.authors = QLabel()
        self.authors.setText('Лукин А.Д.\nФТА???\nПНВ???\nЛаборатория биоинженерии растений ИЦИГ СО РАН')
        form_authors.addRow(self.authors)

        root_layout.addWidget(group_info)
        root_layout.addWidget(group_authors)
        root_layout.addStretch()
        self.setLayout(root_layout)