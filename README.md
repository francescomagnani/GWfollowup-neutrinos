# GWfollowup-neutrinos
Hi, welcome to this repository created from my master's thesis in Nuclear and Subnuclear Physics at the University of Bologna.

On the 24th of May 2023, the existing web of operating gravitational wave (GW) interferometers LIGO/KAGRA/Virgo started the new data-taking run (the run O4). More info in: https://observing.docs.ligo.org/plan/
The aim of the new run is to further understand GW from neutron star binaries as well as to detect for the first time GWs emitted by black hole binaries and burst events, such as supernovae.
A powerful tool is Multi-Messenger Astronomy, which has already proven to the world the possibility of studying cosmic environments using different messengers: cosmic rays, photons, neutrinos, and GWs.

For my master's thesis, I've worked in a neutrino telescope doing what is called GW follow-up, i.e. the analysis of neutrinos detected in correspondence with a GW event. 
The procedure adopted within the collaboration is to look for a single neutrino event at 3sigma of significance, in a time window of 1000 s around the GW arrival time. To do so, the background of the observation has to be reduced to less than 2.7e-3.

With the present application, it is possible to select a GW event and perform cuts trying to reduce the background to the 3 sigma threshold and detect a single ON event. There are 4 possible variables that can be used for the selection, either parallelly or singularly. Just remember that the more complex the analysis the higher the execution time.

# Beta testing
The software has been tested on Linux machines, with python3.x versions. For other configurations, I cannot guarantee it works. Therefore, I strongly suggest to create an Ubuntu 20.4 virtual machine (Desktop ISO image: https://releases.ubuntu.com/20.04/ubuntu-20.04.6-desktop-amd64.iso) or to use the Windows Linux Subsystem (WSL) on Windows machines.

# Prerequisites
1. Python3.x installed, not older. 
2. Your Python version should already have the following libraries installed: sys, os, subprocess, operator, itertools, random, time, math, PIL, and tkinter (sudo apt-get install python3.10-tk, [10MB]). If not present, please install them.
3. pip installed: sudo apt install python3-pip [199MB]
4. XXXX of free memory
   
# How to install and start the application
1. Download the GitHub repository. For memory reasons, the simulation dataset of the events that are used by the application has been uploaded on the drive. Download and save it in the "datasets" folder: https://drive.google.com/file/d/1Vh-vh0Ph1sfnGNVSEv5OiO5FMQJxKSUz/view?usp=drive_link
2. Since the software uses different libraries, let's install them by starting the setup.py file with python: python setup.py
4. Now you're ready to start the application. Inside the main folder digit: python3.x GWfollowup.py --config configuration/optimizationConfig.toml --iter 0

# How to use the application


