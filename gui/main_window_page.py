from pathlib import Path

from PySide6.QtCore import (
    Qt,
    QSize
)
from PySide6.QtWidgets import(
    QWidget,
    QHBoxLayout,
    QListWidgetItem,
    QListWidget,
    QStackedWidget,
)
from PySide6.QtGui import (
    QIcon
)

from gui.settings_page import SettingsPage
#from gui.CRISPOR_scoring_page import CrisporAnalyzerPage
from gui.home_page import HomePage
from gui.guide_parser_page import GuideParserPage
from gui.td_structure_imager import TDStructureImager

from core.config_manager import load_configs

ASSETS = Path(__file__).parent.parent / "assets"

MENU_ITEMS = [
    ('cambodian_resort', 'Home'),
    ('gear_symmetry', 'Settings'),
    ('RNA_2_var', 'Guide parser'),
    ('RNA_cal', '2d imager'),
]


class MainPage(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(1200, 800)
        self.dark_theme = bool(load_configs('main_configs', 'dark_theme'))
        self._build_ui()

    def _icon(self, base_name: str) -> QIcon:
        variant = 'dark' if self.dark_theme else 'light'
        return QIcon(str(ASSETS / f'{base_name}_{variant}.ico'))

    def apply_theme_icons(self):
        for item, base_name in self._icon_items:
            item.setIcon(self._icon(base_name))

    def _on_theme_changed(self, dark: bool):
        self.dark_theme = dark
        self.apply_theme_icons()

    def _build_ui(self):
        self.menu_widget = QListWidget()
        self.menu_widget.setObjectName('Menu')
        self.menu_widget.setUniformItemSizes(True)
        self.menu_widget.setFixedWidth(50)
        self.menu_widget.setIconSize(QSize(40, 40))
        self.menu_widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self._icon_items = []
        for base_name, tooltip in MENU_ITEMS:
            item = QListWidgetItem()
            item.setSizeHint(QSize(52, 52))
            item.setToolTip(tooltip)
            item.setTextAlignment(Qt.AlignCenter)
            self.menu_widget.addItem(item)
            self._icon_items.append((item, base_name))

        self.home_item = self._icon_items[0][0]
        self.settings_item = self._icon_items[1][0]
        self.guide_parser_item = self._icon_items[2][0]
        self.td_structure_imager_item = self._icon_items[3][0]

        self.apply_theme_icons()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setAlignment(Qt.AlignLeft)

        self.home_widget = HomePage()
        self.settings_page = SettingsPage()
        #self.crispor_page = CrisporAnalyzerPage()
        self.guide_parser_page = GuideParserPage()
        self.td_structure_page = TDStructureImager()

        self.settings_page.theme_changed.connect(self._on_theme_changed)

        self.stack = QStackedWidget()
        self.stack.addWidget(self.home_widget)
        self.stack.addWidget(self.settings_page)
        #self.stack.addWidget(self.crispor_page)
        self.stack.addWidget(self.guide_parser_page)
        self.stack.addWidget(self.td_structure_page)

        self.menu_widget.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.menu_widget.setCurrentRow(0)

        layout.addWidget(self.menu_widget)
        layout.addWidget(self.stack)
        self.setLayout(layout)
