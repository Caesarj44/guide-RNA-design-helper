import sys
from pathlib import Path

from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt

from gui.main_window_page import MainPage
from core.config_manager import load_configs, load_style

# Добавляем корень проекта в sys.path (для запуска из любой директории)
ROOT = Path(__file__).parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))



ASSETS = ROOT / 'assets'    

def main():
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )

    app = QApplication(sys.argv)
    app.setApplicationName('guide RNA design helper v1.1.1')
    app.setApplicationVersion('1.1.1')
    app.setOrganizationName('Nerv')

    cfg = load_configs()
    is_dark = cfg.get('main_configs', {}).get('dark_theme', False)
    print(cfg)
    print(is_dark)
    qss_file = 'dark.qss' if is_dark else 'light.qss'
    style = load_style(str(ASSETS / qss_file))
    if style:
        app.setStyleSheet(style)
    
    window = MainPage()
    window.showMaximized()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
