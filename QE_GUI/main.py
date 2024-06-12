from PySide2.QtWidgets import QApplication
import sys
import os

# Agregar el directorio QE_GUI a la ruta de búsqueda de módulos
sys.path.append(os.path.dirname(__file__))

from file_operations.verify_qe import check_quantum_espresso
from screens.DialogMessage import DialogMessage
from screens.main_screen import run_main


def main():
    qe_installed = check_quantum_espresso('pw.x')

    if qe_installed:
        run_main()
    else:
        app = QApplication(sys.argv)
        dialog = DialogMessage('Please install Quantum Espresso to continue', 'Error', '#DC3545')
        dialog.exec_()
        sys.exit(0)   



if __name__ == "__main__":
    main()
    """
    qe_installed = check_quantum_espresso('pw.x')

    if qe_installed:
        run_main()
    else:
        app = QApplication(sys.argv)
        dialog = DialogMessage('Please install Quantum Espresso to continue', 'Error', '#DC3545')
        dialog.exec_()
        sys.exit(0)"""