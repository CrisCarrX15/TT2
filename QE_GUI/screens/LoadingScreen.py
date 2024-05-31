from PySide2.QtWidgets import QDialog, QVBoxLayout, QLabel
from PySide2.QtCore import Qt

class LoadingDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Loading")
        self.setWindowModality(Qt.ApplicationModal)  # Lock the entire app
        # Adjust window flags to remove close button
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowCloseButtonHint)

        layout = QVBoxLayout()
        layout.addWidget(QLabel("Running Quantum ESPRESSO, please wait..."))
        self.setLayout(layout)
        self.setFixedSize(400, 100)