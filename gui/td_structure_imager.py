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

class TDStructureImager (QWidget):

    def __init__(self):
        super().__init__()
        self._build_gui()

    def _build_gui(self):
        pass