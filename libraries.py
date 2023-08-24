import sys, os
import subprocess

#main
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'argparse'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tomli==2.0.1'])
# optimization script
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numpy==1.23.5'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'pandas==2.0.0'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'datetime'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib==3.6.0'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'healpy'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'ligo.skymap==1.0.6'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'utm'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'h5py'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'numexpr'])
subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tables'])

# folders
os.mkdir("datasets")
os.mkdir("results")
