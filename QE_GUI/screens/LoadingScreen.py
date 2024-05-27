from PySide2.QtWidgets import QDialog, QVBoxLayout, QLabel

class LoadingDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Cargando")
        self.setModal(True)
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Ejecutando Quantum ESPRESSO, por favor espere..."))
        self.setLayout(layout)
        self.setFixedSize(300, 100)
