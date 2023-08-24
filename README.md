# GWfollowup-neutrinos
Hi, welcome to this repository created from my master's thesis in Nuclear and Subnuclear Physics at the University of Bologna.

On the 24th of May 2023, the existing web of operating gravitational wave (GW) interferometers LIGO/KAGRA/Virgo started the new data-taking run (the run O4). More info in: https://observing.docs.ligo.org/plan/.
The aim of the new run is to further understand GW from neutron star mergers as well as to detect GWs emitted by black hole binaries and burst events such as supernovae.

Multi-Messenger Astronomy is a powerful tool which has already proven the possibility of studying cosmic environments using different messengers: cosmic rays, photons, neutrinos, and GWs.

For my master's thesis, I've worked in a neutrino telescope doing what is called GW follow-up, i.e. the analysis of neutrinos detected in correspondence with a GW event. 
The procedure adopted within the collaboration is to look for a single neutrino event at 3sigma of significance, in a time window of 1000 s around the GW arrival time. To do so, the background of the observation has to be reduced to less than 2.7e-3.

With the present application, it is possible to select a GW event and perform cuts trying to reduce the background to the 3 sigma threshold and detect ON events, even if it is hard to get any events post-cuts and if you do, you might detect the first neutrino ever in correspondence with a GW. Read below to know how to properly run the software.

# Beta testing
The software has been tested on Ubuntu 20.04 with python3.8.10. For other configurations, I cannot guarantee it works. 
Therefore, I strongly suggest to create an Ubuntu 20.4 virtual machine (Desktop ISO image: https://releases.ubuntu.com/20.04/ubuntu-20.04.6-desktop-amd64.iso) or to use the Windows Linux Subsystem (WSL) on Windows machines installing Ubuntu from the play store.

# Prerequisites
0. XXXX of free memory.
1. Before to install any library do: "sudo apt update", "sudo apt upgrade", and if needed also "sudo apt autoremove"
2. Python3.8.10 installed. Your Python version should already have the following libraries installed: sys, os, subprocess, operator, itertools, random, time, and math.
3. Install two more libraries: PIL "sudo apt-get install python3-pil.imagetk" [?], and tkinter "sudo apt-get install python3.10-tk" [10MB].
5. Install pip "sudo apt install python3-pip" [199MB], and update it "pip install --upgrade pip" (I work with python3-pip23.1.2)
   
# How to install and start the software
1. Download the GitHub repository. For memory reasons, the dataset with the events has been uploaded on the drive. Download and save it in the "datasets" folder: https://drive.google.com/file/d/1Vh-vh0Ph1sfnGNVSEv5OiO5FMQJxKSUz/view?usp=drive_link
2. Install python libraries needed; in the software folder run: python3.x libraries.py
4. Now you're ready to start the application. Inside the main folder digit: python3.x GWfollowup.py --config configuration/optimizationConfig.toml --iter 0

# How to use the application


