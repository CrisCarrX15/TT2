from PySide2.QtWidgets import QApplication, QDialog, QVBoxLayout, QLabel
from PySide2.QtGui import QPalette, QColor, QFont
from PySide2.QtWidgets import QTextEdit, QPushButton
from PySide2.QtCore import Qt

class DialogMessage(QDialog):
    def __init__(self, text, title, color):
        super().__init__()
        self.setWindowTitle(title)
        self.setWindowModality(Qt.ApplicationModal)
        self.setMinimumWidth(400)  # Establecer el ancho mínimo del diálogo

        layout = QVBoxLayout()

        text_edit = QTextEdit()
        text_edit.setPlainText(text)
        text_edit.setReadOnly(True)  # Hacer el QTextEdit de solo lectura
        text_edit.setFont(QFont("Arial", 12))  # Establecer un tamaño de fuente más grande
        # Set the background color using QSS
        self.setStyleSheet(f"background-color: {color};")
        layout.addWidget(text_edit)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

        self.setLayout(layout)

