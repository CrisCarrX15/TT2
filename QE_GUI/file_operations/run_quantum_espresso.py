import sys
import os
import subprocess

class RunQuantumEspresso():

    def run_qe_process(self, file_path_in, file_path_out):
            print(f"Execute process")
            subprocess.run(f'pw.x < {file_path_in} > {file_path_out}', shell=True)