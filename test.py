import sys
from PySide6.QtCore import Qt, QSize, QSettings
from PySide6.QtWidgets import (QApplication, QLabel, QWidget, QListWidget,
                                QListWidgetItem, QPushButton, QVBoxLayout,
                                QHBoxLayout, QStackedWidget, QFormLayout,
                                QCheckBox)
from PySide6.QtGui import QIcon


# ================================================================
# Загрузка тем из файлов
# ================================================================

def load_style(filename: str) -> str:
    """Читает QSS файл и возвращает строку со стилями."""
    try:
        with open(filename, "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        print(f"Файл стилей {filename} не найден, используется пустой стиль.")
        return ""


LIGHT_STYLE = load_style("light.qss")
DARK_STYLE  = load_style("dark.qss")


# ================================================================
# Страницы
# ================================================================

class MainPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Главная страница"))
        self.setLayout(layout)


class SettingsPage(QWidget):
    def __init__(self):
        super().__init__()
        layout = QFormLayout()

        # --- Чекбокс тёмной темы ---
        self.dark_mode_cb = QCheckBox()
        self.dark_mode_cb.stateChanged.connect(self._toggle_theme)
        layout.addRow("Тёмная тема:", self.dark_mode_cb)

        self.setLayout(layout)

        # Загружаем сохранённое значение и применяем тему при старте
        self._load_settings()

    def _load_settings(self):
        """Читает сохранённые настройки и применяет тему."""
        s = QSettings("МояЛаба", "GRDH")
        is_dark = s.value("dark_mode", False, type=bool)

        # Блокируем сигнал чтобы setChecked не вызвал _toggle_theme дважды
        self.dark_mode_cb.blockSignals(True)
        self.dark_mode_cb.setChecked(is_dark)
        self.dark_mode_cb.blockSignals(False)

        # Применяем тему вручную
        self._apply_theme(is_dark)

    def _toggle_theme(self, state):
        """Вызывается при изменении чекбокса."""
        is_dark = (state == Qt.Checked.value)
        self._apply_theme(is_dark)

        # Сохраняем выбор
        s = QSettings("МояЛаба", "GRDH")
        s.setValue("dark_mode", is_dark)

    def _apply_theme(self, is_dark: bool):
        """Применяет нужный QSS ко всему приложению."""
        app = QApplication.instance()
        app.setStyleSheet(DARK_STYLE if is_dark else LIGHT_STYLE)


# ================================================================
# Главное окно
# ================================================================

class Widget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        # --- Боковое меню ---
        self.menu_widget = QListWidget()
        self.menu_widget.setUniformItemSizes(True)
        self.menu_widget.setFixedWidth(50)
        self.menu_widget.setIconSize(QSize(40, 40))
        self.menu_widget.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # Пункт: Главная (индекс 0)
        home_item = QListWidgetItem()
        home_item.setSizeHint(QSize(48, 48))
        home_item.setIcon(QIcon("home48x48.ico"))
        self.menu_widget.addItem(home_item)

        # Пункт: Настройки (индекс 1)
        settings_item = QListWidgetItem()
        settings_item.setSizeHint(QSize(48, 48))
        settings_item.setIcon(QIcon("settings48x48.ico"))
        self.menu_widget.addItem(settings_item)

        # --- Стек страниц ---
        self.stack = QStackedWidget()
        self.stack.addWidget(MainPage())     # индекс 0 — совпадает с home_item
        self.stack.addWidget(SettingsPage()) # индекс 1 — совпадает с settings_item

        # Переключение страниц по клику на меню
        self.menu_widget.currentRowChanged.connect(self.stack.setCurrentIndex)

        # --- Главный layout ---
        main_layout = QHBoxLayout()
        main_layout.addWidget(self.menu_widget)
        main_layout.addWidget(self.stack)
        self.setLayout(main_layout)


# ================================================================
# Точка входа
# ================================================================

if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = Widget()
    w.resize(800, 600)
    w.show()

    sys.exit(app.exec())