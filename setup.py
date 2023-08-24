import sys
import subprocess

#main
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'argparse'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tomli==2.0.1'])
# optimization script
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy==1.23.5'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas==2.0.0'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'datetime'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'healpy'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'ligo.skymap==1.0.6'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'astropy==5.1'])
# skymaps
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib==3.6.0'])
# GUI
#subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tk'])
#subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'PIL'])

