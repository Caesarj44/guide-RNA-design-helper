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
from gui.CRISPOR_scoring_page import CrisporAnalyzerPage
from gui.home_page import HomePage

ASSETS = Path(__file__).parent.parent / "assets"
print(ASSETS)

class MainPage(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(1200,800)
        self._build_ui()

    def _build_ui(self):
        #=== боковое меню
        self.menu_widget = QListWidget()
        self.menu_widget.setObjectName('Menu')
        self.menu_widget.setUniformItemSizes(True)
        self.menu_widget.setFixedWidth(50)
        self.menu_widget.setIconSize(QSize(40,40))
        self.menu_widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff) # отключение бесящего скролбара
        

        def _add_gui_menu_item(icon_name: str, tooltip: str):
            item = QListWidgetItem()
            item.setSizeHint(QSize(52,52))
            icon_path = str(ASSETS / icon_name)
            item.setIcon(QIcon(icon_path))
            item.setToolTip(tooltip)
            item.setTextAlignment(Qt.AlignCenter)
            self.menu_widget.addItem(item)
        
        _add_gui_menu_item('home48x48.ico', 'home_page')
        _add_gui_menu_item('settings48x48.ico', 'settings_page')
        _add_gui_menu_item('settings16x16.ico', 'crispor_page')
    
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setAlignment(Qt.AlignLeft)
        

        self.home_widget = HomePage()
        self.settings_page = SettingsPage()
        self.crispor_page = CrisporAnalyzerPage()
        
        self.stack = QStackedWidget()
        self.stack.addWidget(self.home_widget)
        self.stack.addWidget(self.settings_page)
        self.stack.addWidget(self.crispor_page)
        

        self.menu_widget.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.menu_widget.setCurrentRow(0)
        
        layout.addWidget(self.menu_widget)
        layout.addWidget(self.stack)
        self.setLayout(layout)
