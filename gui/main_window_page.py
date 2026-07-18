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

ASSETS = Path(__file__).parent.parent / "assets"
print(ASSETS)

class MainPage(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(1200,800)
        self._build_ui()
        # self.home_item : QListWidgetItem() 
        # self.settings_item : QListWidgetItem()
        # self.guide_parser_item : QListWidgetItem() 
        # self.test_home_item : QListWidgetItem()
        

    def _build_ui(self):
        self.menu_widget = QListWidget()
        self.menu_widget.setObjectName('Menu')
        self.menu_widget.setUniformItemSizes(True)
        self.menu_widget.setFixedWidth(50)
        self.menu_widget.setIconSize(QSize(40,40))
        self.menu_widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        

        def _add_gui_menu_item(icon_name: str, tooltip: str):
            item = QListWidgetItem()
            item.setSizeHint(QSize(52,52))
            icon_path = str(ASSETS / icon_name)
            item.setIcon(QIcon(icon_path))
            item.setToolTip(tooltip)
            item.setTextAlignment(Qt.AlignCenter)
            self.menu_widget.addItem(item)
    
        self.home_item = QListWidgetItem()
        self.home_item.setSizeHint(QSize(52,52))
        home_icon_path = str(ASSETS / 'cambodian_resort_var2.ico')
        self.home_item.setIcon(QIcon(home_icon_path))
        self.home_item.setToolTip('Home')
        self.home_item.setTextAlignment(Qt.AlignCenter)
        self.menu_widget.addItem(self.home_item)


        self.settings_item = QListWidgetItem()
        self.settings_item.setSizeHint(QSize(52,52))
        settings_icon_path = str(ASSETS / 'gear_symmetry_var2.ico')
        self.settings_item.setIcon(QIcon(settings_icon_path))
        self.settings_item.setToolTip('Settings')
        self.settings_item.setTextAlignment(Qt.AlignCenter)
        self.menu_widget.addItem(self.settings_item)

        self.guide_parser_item = QListWidgetItem()
        self.guide_parser_item.setSizeHint(QSize(52,52))
        guide_parser_icon_path = str(ASSETS / 'RNA_2_var.ico')
        self.guide_parser_item.setIcon(QIcon(guide_parser_icon_path))
        self.guide_parser_item.setToolTip('Guide parser')
        self.guide_parser_item.setTextAlignment(Qt.AlignCenter)
        self.menu_widget.addItem(self.guide_parser_item)



        #_add_gui_menu_item('home48x48.ico', 'Home')
        #_add_gui_menu_item('gear_symmetry_var2.ico', 'Settings')
        #_add_gui_menu_item('settings16x16.ico', 'crispor_page')
        #_add_gui_menu_item('gear_symmetry_var1.ico', 'Guide parser')
    
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setAlignment(Qt.AlignLeft)
        

        self.home_widget = HomePage()
        self.settings_page = SettingsPage()
        #self.crispor_page = CrisporAnalyzerPage()
        self.guide_parser_page = GuideParserPage()
        
        self.stack = QStackedWidget()
        self.stack.addWidget(self.home_widget)
        self.stack.addWidget(self.settings_page)
        #self.stack.addWidget(self.crispor_page)
        self.stack.addWidget(self.guide_parser_page)
        

        self.menu_widget.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.menu_widget.setCurrentRow(0)
        
        layout.addWidget(self.menu_widget)
        layout.addWidget(self.stack)
        self.setLayout(layout)
