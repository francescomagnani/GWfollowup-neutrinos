# GWfollowup-neutrinos
Hi, welcome to this repository created from my master's thesis in Nuclear and Subnuclear Physics at the University of Bologna.

On the 24th of May 2023, the existing network of operating gravitational wave (GW) interferometers LIGO/KAGRA/Virgo started the new data-taking run (the run O4). More info in: https://observing.docs.ligo.org/plan/.
The aim of O4 is to further understand GW from neutron star mergers as well as to detect GWs emitted by black hole binaries and burst events such as supernovae.

**This is not a software to detect GWs**, instead it exploits the dynamic and powerful tool of Multi-Messenger Astronomy (MMA). MMA is a field of astronomy, which takes advantage of all of the 4 types of cosmic carriers to investigate the cosmos: cosmic rays, photons, neutrinos, and GWs.  
**This software simulates the behaviour of what will soon become the world's largest neutrino telescope in performing neutrino searches after the emission of a GW alerts by LIGO** (Laser Interferometer Gravitational-wave Observatory). The way of working is to detect a single 3sigma-significance neutrino event in a time window of 1000 s around the GW alert arrival time. 
The analysis is entirely performed in Python, both with the use of classical and _machine learning_ techniques.

# The Concept
Performing MMA studies allows one to possibly detect more information regarding the astrophysical system observed, e.g. stars, galaxies, black holes. In particular, the use of neutrino in astronomy is highly profitable since they can escape very dense environments due to their low cross-section, hence, proving the interiors of exotic objects in the cosmos. Still, no-one has ever detected neutrinos in correspondence of a GW.

Therefore, here I exploit the power of an underwater neutrino telescope to perform what is called the **GW follow-up**, i.e. the search for neutrinos after a GW event. 
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
When you start the software the following window will appear.  
<img src="https://github.com/francescomagnani/GWfollowup-neutrinos/assets/75760916/dd525b37-1379-4f59-8653-b5f3628a7480" width="400">  
Firstly, click on the **Browse file** button and select one of the three files already present in the _alerts_ folder, from the window that will open. This will allow the software to read and interpret the selected GW alert (one can even upload a ".fits" file of its own).
<img src="https://github.com/francescomagnani/GWfollowup-neutrinos/assets/75760916/12405692-cf79-497a-ac86-dee850e317dc" width="400">
Wait for the upload of the picture on the main window. Here are **shown all the down-going events present between December 5 and 18, 2022**.  
<img src="https://github.com/francescomagnani/GWfollowup-neutrinos/assets/75760916/5145dfcd-1e29-46b9-83d4-1cf9812a119f" width="400">  
  
Subsequently, you may define the parameters for the analysis.  
Select at least one of the optimization variables. Remember that the more the selected variables the higher the execution time. The software will perform cuts on the selected variables iteratively, moving within the ranges specified in the configuration file, i.e. _optimizationConfig.toml_, here also shown:
- likelihood: [100,750), step = 25
- tracklength [m]: [200,950), step = 20
- n. of hits: [50,1000), step = 50
- bdt score: [log<sub>10</sub>(-5),log<sub>10</sub>(-0.1)), step = 0.1  
Since some values in these ranges are strong enough to cut all the data present in the dataset (read **The Concept** section to understand why strong cuts are needed), the algorithm will start the selection by cutting on safer ranges, which will then be optimized on the broader ranges above. The safe ranges are shown below:
- likelihood: [100,300)
- tracklength: [300,550)
- n. of hits: [50,200)
- bdt score: [log<sub>10</sub>(-0.3),log<sub>10</sub>(-0.1))
  
Finally, to make the algorithm faster, one can set the _phase space_ and _order_ parameters accordingly. The first one tells the software which is the fraction of values to be considered among the safe ranges. If _phase space_ = 1, all the values in the safe ranges will be considered; if _phase space_ = 0.5 only half of them will be used. The _order_ parameter tells the software in which order the variables will be optimized. If _order_ = 1, all the possible orders will be considered, e.g. with 2 variables the orders are 2: first variable A, then B; or vice-versa; if _order_ = 0.5 only the order _variable A then variable B_ will be considered.  
  
At this point you are ready to start the software by pressing the **Start** button. At the end of the optimization, the image will be updated **showing only the events present in a 1000 s window around the GW event** after the cuts.
<img src="https://github.com/francescomagnani/GWfollowup-neutrinos/assets/75760916/e598754e-1126-44d8-9c65-0016641dfeb7" width="400">
