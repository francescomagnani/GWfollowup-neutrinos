# GWfollowup-neutrinos
Hi, welcome to this repository created from my master's thesis in Nuclear and Subnuclear Physics at the University of Bologna.

On the 24th of May 2023, the existing web of operating gravitational wave (GW) interferometers LIGO/KAGRA/Virgo started the new data-taking run (the run O4). More info in: https://observing.docs.ligo.org/plan/.
The aim of O4 is to further understand GW from neutron star mergers as well as to detect GWs emitted by black hole binaries and burst events such as supernovae.

**This is not a software to detect GWs**, instead it exploits the dynamic and powerful tool of Multi-Messenger Astronomy (MMA). MMA is a field of astronomy, which takes advantage of all the 4 types of cosmic carriers to investigate the cosmos: cosmic rays, photons, neutrinos, and GWs.
**This software** simulates the behaviour of what will soon become the world's largest neutrino telescope in performing neutrino searches after the emission of a GW alert from LIGO (Laser Interferometer Gravitational-wave Observatory). The way of working is similar to the procedure adopted within the collaboration: a single 3sigma-significance neutrino event is searched for in a time window of 1000 s around the GW alert arrival time. 
The analysis is performed entirely in Python, also with the use of a machine learning technique called Binary Decision Tree, as explained in the "Concept" section.

# The Concept
Performing MMA studies allows one to possibly detect more information regarding the astrophysical system observed, e.g. stars, galaxies, black holes. In particular, the use of neutrino in astronomy is highly profitable since they can escape very dense environment, hence, proving the interiors of exotic objects in the cosmos. Still, no-one has ever detected neutrinos in correspondence of a GW, which would result in a great new discovery.

Therefore, we exploit the power of an underwater neutrino telescope to perform what is called the GW follow-up, i.e. the search for neutrinos after a GW event. 
Usually, neutrino telescopes perform analyses with neutrino candidates coming from _under_ the experiment (the _up-going sky_), i.e. particles which have traversed Earth. This choice is justificated by the extremely high atmospheric muon contamination, which is instead suppressed by the Earth's thickness. The problem of atmospheric muons is that they enter the detector and leave a similar track to that of neutrino events.
Here, instead, **I** wanted to **perform GW follow-up with events coming from _above_ the experiment** (the _down-going sky_), i.e. events which have crossed the Earth's atmosphere and approximately 3 km of water to reach the detector on the seabed. This would increase the visibility of the detector and help to investigate sources which pass from the up-going to the down-going sky in the moment of the emission (_transient sources_), hence are difficult to study.

Contrarily from the _up-going_ sky, the muon contamination is so high that strong cut on datasets are needed to move closer to the 3 sigma significance threshold. For this reason, a study on the strongest variables to use in the suppression of the background has been performed and are: likelihood of the track, track length, track number of hits in the detector, and bdt score. For each event, the last variable describes its probability of being a muon or a neutrino crossing the detector, computed by a BDT model trained on simulated events.

# Beta testing
This software has been tested on Ubuntu 20.04 with python3.8.10. For other configurations, its functionality is not guaranteed.
Therefore, I strongly suggest to create an Ubuntu 20.4 virtual machine (Desktop ISO image: https://releases.ubuntu.com/20.04/ubuntu-20.04.6-desktop-amd64.iso) or to use the Windows Linux Subsystem (WSL) on Windows machines installing Ubuntu from the play store.

# Prerequisites
0. Bewteen 500-1000 MB of free memory.
1. It is always a good time to run: <sup>sudo apt update</sup>, <sup>sudo apt upgrade</sup>, and if needed also <sup>sudo apt autoremove</sup>.
2. _Python3.8.10_ installed. Your Python version should already have the following libraries installed: _sys, os, subprocess, operator, itertools, random, time_, and _math_.
3. Install three more libraries:
	+ _astropy_: <sup>sudo apt-get install python3-astropy</sup>,
	+ _PIL_: <sup>sudo apt-get install python3-pil.imagetk</sup>,
	+ _tkinter_: <sup>sudo apt-get install python3-tk<sup>.
5. Install _pip_ <sup>sudo apt install python3-pip</sup>, and update it with <sup>pip install --upgrade pip</sup> (I work with _python3-pip23.1.2_)
6. In the case OpenSSL and cryptography packages are already installed in your Python version, you might have a problem with the installation of the software. Please, remove _OpenSSL_ packages manually:
	+ <sup>sudo rm -rf /usr/local/lib/python3.8/dist-packages/OpenSSL</sup>
	+ <sup>sudo rm -rf /usr/local/lib/python3.8/dist-packages/pyOpenSSL...VERSION...</sup>
	+ <sup>sudo rm -rf /home/YOU/.local/lib/python3.8/site-packages/OpenSSL</sup>
	+ <sup>sudo rm -rf /home/YOU/.local/lib/python3.8/site-packages/pyOpenSSL...VERSION...</sup>  
Then re-install _OpenSSL_ (version 22.0.0 works for me): <sup>pip install pyOpenSSL==22.0.0</sup>  
Finally, force re-installation of cryptography to the version 38.0.4 (it works fine): <sup>pip install --force-reinstall "cryptography==38.0.4"</sup>
   
# How to install and start the software
1. Download the GitHub repository. For memory reasons, the dataset with the telescope events has been uploaded on the drive. Download and save it in the "datasets" folder: **https://drive.google.com/file/d/1Vh-vh0Ph1sfnGNVSEv5OiO5FMQJxKSUz/view?usp=drive_link**
2. Install python libraries needed by the software by running: <sup>python3.x libraries.py</sup>
4. Now you're ready to start the application. Inside the main folder digit: <sup>python3.x GWfollowup.py --config configuration/optimizationConfig.toml --iter 0</sup>

# Play with the software
![app1](https://github.com/francescomagnani/GWfollowup-neutrinos/assets/75760916/dd525b37-1379-4f59-8653-b5f3628a7480)

According to the variables selected, the software performs the optimization of a certain number of initial cuts (specified in the _phase space_
