import subprocess

def check_quantum_espresso(command='pw.x'):
    installed = False
    try:
        # Run the command and redirect the output to avoid waiting for input
        result = subprocess.run(command, input="quit\n", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=5)
        if result.returncode == 0 or "Program PWSCF" in result.stdout:
            installed = True
        else:
            installed = False
    except subprocess.TimeoutExpired:
        installed = True  # We assume it is installed if there is a timeout
    except FileNotFoundError:
        installed = False

    return installed