import sys, os
import operator
from itertools import permutations, product
import random
import numpy as np
import pandas as pd
from datetime import datetime
import time
import healpy as hp
from ligo.skymap.io import fits
from astropy.units import Quantity
from astropy.time import Time

from selection.skymaps import get_region_90, extend_region_90, convert_region_to_local, define_zenith_bands, plot_skymap_healpix


def enablePrint():
    sys.stdout = sys.__stdout__

def blockPrint():
    sys.stdout = open(os.devnull, 'w')



class Opt:
    def __init__(self,configu, it):
        self.config = configu
        self.iteration = it
        self.outdir = configu['OUTPUT_DIR']
        self.indir = configu['INPUT_DIR']

    def read_GW(self,filenam):
        # load the events from some existing h5 file
        self.dataframe = pd.read_hdf(self.config['DATA_DIR'])
        # compute livetimes
        self.livetime_sig = self.config['SIGNAL_lvt'] # [s] ITERATION 0
        self.livetime_bkg = self.config['BKG_lvt'] # [s]
        fits_file = filenam
        self.ID = "GW"+fits_file.split('_')[0].split('GW')[3]
        # make result folder
        self.res_fold = './results/'+self.ID
        if not os.path.exists(self.res_fold):
            os.mkdir(self.res_fold)
        # read alert map
        fits_header = fits.read_sky_map(fits_file, nest=False)[1]
        self.data_str = "2022-12-18T23:50:00" # fake date
        self.evt_date = Time(self.data_str, format="isot")
        self.skymap = fits.read_sky_map(fits_file, nest=False)[0]
        self.skymap = hp.ud_grade(self.skymap, 64)  # downgrade a bit the resolution to make it quicker
        nside = hp.npix2nside(len(self.skymap))
        # equatorial coordinates into pixels
        self.dataframe["pixel_eq"] = hp.ang2pix(
            nside,
            self.dataframe["trackfit_ra"].to_numpy(),
            self.dataframe["trackfit_dec"].to_numpy(),
            lonlat=True,
        )
        # definition of the start and end times of the search
        self.timesig_start = self.evt_date + Quantity(self.config["TIMEWINDOWS"][self.iteration][0])
        self.timesig_stop = self.evt_date + Quantity(self.config["TIMEWINDOWS"][self.iteration][1])
        self.timebkg_start = self.evt_date + Quantity(self.config["TIMEWINDOW_BKG"][0])
        self.timebkg_stop = self.evt_date + Quantity(self.config["TIMEWINDOW_BKG"][1])

        # get 90% region (list of pixels to consider)
        self.region90 = get_region_90(self.skymap)
        # get extended 90% region = RoI (list of pixels to consider)
        self.region90ext = extend_region_90(self.region90, nside, self.config["ANGULAR_EXTENSION"])
        # convert it to local coordinates (skymap with p_i = how often pixel i is inside the RoI over the signal time window)
        self.skymap_region90ext_local = convert_region_to_local(
            self.region90ext, nside, self.timesig_start, self.timesig_stop, self.config["DETECTOR"])
        # define the corresponding theta bands (OFF bands)
        self.pixelbands, self.omegaoff, self.omegaon = define_zenith_bands(self.skymap_region90ext_local, size_bands_deg=5)
        self.nbands = len(self.pixelbands)

        # get pixel in local coordinates for each neutrino event
        phi = np.arctan2(self.dataframe['dir_y'].to_numpy(), self.dataframe['dir_x'].to_numpy())
        theta = np.arccos(self.dataframe['dir_z'].to_numpy())
        self.dataframe["pixel_local"] = hp.ang2pix(nside, theta, phi)

        # mask on and off on data
        self.mask_OFF = (
            (self.dataframe['mjd'] >= Time(self.timebkg_start,format='isot').to_value('mjd')) 
            & (self.dataframe['mjd'] <= Time(self.timebkg_stop,format='isot').to_value('mjd'))
            &  (~self.dataframe["pixel_eq"].isin(self.region90ext)) 
        )
        self.mask_ON = (
            (self.dataframe["mjd"] >= Time(self.timesig_start,format='isot').to_value('mjd'))
            & (self.dataframe["mjd"] <= Time(self.timesig_stop,format='isot').to_value('mjd'))
            & (self.dataframe["pixel_eq"].isin(self.region90ext))
        )
        self.mask_allON = (
            (self.dataframe["mjd"] >= Time(self.timesig_start,format='isot').to_value('mjd'))
            & (self.dataframe["mjd"] <= Time(self.timesig_stop,format='isot').to_value('mjd'))
        )

        self.outfig_skymap = os.path.join(self.res_fold, f"Skymap_{self.config['NAME']}_{self.ID}_iter_{self.iteration}.png")
        plot_skymap_healpix(self.skymap, self.region90ext, self.dataframe[self.mask_allON], self.evt_date, self.outfig_skymap, self.ID) # print and save img



    # DEFINE SCORE FUNCTION
    def score_function(self,omegaonn_i,omegaofff_i,dataframee): 
        nevts_i=len(dataframee)
        bkg_i = nevts_i * (omegaonn_i/omegaofff_i) * (self.livetime_sig/self.livetime_bkg)
        return 1 / (bkg_i+1), nevts_i, bkg_i

    def condition_function(self,bkg_thr_i):
        return 1 / (bkg_thr_i+1)


    def epsilon_op(self,algo_mod,a,b):
        if algo_mod == "boost_cut":
            return operator.ge(a,b)
        else:
            return operator.gt(a,b)
        # feel free to program other options

    # OPTIMIZATION CYCLEs -------------------------------------------------------------------------------------------------------------------------------------------
    def opt_selection(self,bkg_thr_i,omegaonnn_i,omegaoffff_i,dataframeee,mode_print):
        if mode_print == "debug":
            enablePrint()
        else:
            blockPrint()
        # open the configuration file: take variables, ranges, saferanges.
        print("Starting optimization")
        # controls?
        #variables = list(self.config["VARIABLES"]) # optimizing variables
        print("variables: ",self.variables)
        null_res = list(np.zeros(len(self.variables))) # in case the algo doesn't find a good cut, it returns this
        print("null list: ",null_res)
        if int(self.config["GENERAL"]["linearity"]) != 1:
            print("Apologies human =( \n but this kind of linearity has not been developed yet. \n If you want to give the model another linearity (different than 1) you have to develop the corresponding code.")
            enablePrint()
            return null_res
        perm = list(permutations(np.arange(0,len(self.variables),1))) # order in which the variables will be optimized
        perms = []
        for p in perm:
            perms.append(list(p))
        S = len(perms) # all the possible orders of the selected variables
        print("permutations, possible stories = ",perms,S)
        directions = {">" : operator.gt, "<" : operator.lt}
        saferanges = []
        ranges = []
        for rang in self.variables:
            print("rang = ",rang,self.config["VARIABLES"][rang])
            saferanges.append(list(np.round(np.arange(self.config["VARIABLES"][rang][1],self.config["VARIABLES"][rang][2],self.config["VARIABLES"][rang][4]),1))) # saferanges from where to select the random cuts
            ranges.append(list(np.round(np.arange(self.config["VARIABLES"][rang][0],self.config["VARIABLES"][rang][3],self.config["VARIABLES"][rang][4]),1))) # whole ranges where to optimize variables
        print("saferanges = ",saferanges,len(saferanges))
        print("ranges = ",ranges,len(ranges))
        R = 1 # it will be the number of possible random initial cuts to try
        for i in range(len(saferanges)):
            R *= len(saferanges[i])
        print("number of possible initial safe cuts", R)
        # set all the starting cuts + variables' optimization orders
        perm_safe = list(product(*saferanges))
        perms_safe = []
        for i in perm_safe:
            perms_safe.append(list(i))
        cuts_rnd = []
        #cuts_rnd_order = []
        print("starting point = ",self.stpoint)
        print("orders = ",self.ord)
        for j in range(int(R*float(self.stpoint))):
            newcut = random.choice(perms_safe)
            print("new random: ",newcut)
            perms_safe.remove(newcut)
            print("remove from available safecuts: ",perms_safe)
            stories = random.sample(perms,int(S*float(self.ord)))
            print("orders: ",stories)
            for st in stories:
                cuts_rnd.append([newcut,st])
        print("safe cuts chosen: ",cuts_rnd)
        # condition if condition mod
        epsiloncondition = ''
        if self.config["GENERAL"]["modality"] ==  "condition":
            epsiloncondition = self.condition_function(bkg_thr_i)
        print("condition modality on: ",epsiloncondition)
        
        # initial settings
        cuts = [] # will contain all the already tried cuts + variables orders
        epsilon = 0 # the initial efficiency function
        epsilon_best = 0 # the best efficiency function reached (the minimum)
        c_best = []
        nbkg_best = 0
        bkg_best = 0

        print("settings applied, let's start the algorithm...")

        # algo --------------------------------------------------------------------------------------------------------------- M I L L E N N I A ---------------
        for evt in cuts_rnd:
            # settings for the epoch cycle
            u = 8 # will tell the direction onto which to cut (up or down along a certain variable's range)
            pins = []
            cut_labs = [] # will contain the best cut found every millennium
            
            # start the new millennium with random cuts (then improve it through the epoch cycle)
            for i in range(0,len(saferanges)):
                pins.append([1,1]) # both the directions are available (and possible to optimize)

            # start with a cut + story
            cut_labs = evt
            print("cut + order to apply: ",cut_labs)
            #apply the first cut
            cuts.append(cut_labs) # update overall cuts
            print("apply initial cut + order ",cut_labs[0],cut_labs[1])
            print("dataframe size before cut: ",dataframeee.shape[0]) # (pixelbands_i,livetime_sigg,livetime_bkgg,omegaonn_i,omegaofff_i,dataframee,bkg_sta,bkg_sto,region90extt): 
            epsilon_pre, numberbkgeventspre, bkgeventspre = self.score_function(omegaonnn_i,omegaoffff_i,dataframeee)
            print("dataframe nbkg, bkg, nsgn, sgn PRE CUT: ",numberbkgeventspre, bkgeventspre)
            dataframeee_cut = dataframeee
            for i in cut_labs[1]:
                print(self.config["VARIABLES"][self.variables[i]][5],self.variables[i],cut_labs[0][i])
                if self.config["VARIABLES"][self.variables[i]][6] == "log10":
                    print("log10")
                    dataframeee_mill = dataframeee_cut.loc[directions[self.config["VARIABLES"][self.variables[i]][5]](np.log10(dataframeee_cut[self.variables[i]]),cut_labs[0][i])]
                    dataframeee_cut = dataframeee_mill
                else:
                    print("no log ")
                    dataframeee_mill = dataframeee_cut.loc[directions[self.config["VARIABLES"][self.variables[i]][5]](dataframeee_cut[self.variables[i]],cut_labs[0][i])]
                    dataframeee_cut = dataframeee_mill
                # fare <> e fare log
            print("dataframe size post cut: ",dataframeee_mill.shape[0])
            # compute epsilon 
            #1 / (bkg_i+1), nevts_i, bkg_i, n_ON_i, sgn_i
            epsilon, numberbkgevents, bkgevents = self.score_function(omegaonnn_i,omegaoffff_i,dataframeee_mill)
            print("score function = ",epsilon, numberbkgevents, bkgevents)
            if epsilon == 1: # if cut too strong, reset it
                epsilon = 0
                print("cut too strong!")

            epsilon_M = epsilon # this will remain fixed for the whole duration (we can perform checks on this)

            # select the optimizing variable and the direction
            for index in cut_labs[1]: # starts epoch cycle -------------------------------------------------------- S T O R Y -------------------------------------
                print("next variable: ",self.variables[index])

                if self.config["GENERAL"]["mod"] == "boost_cut": # if with boost mod, then try to boost all the cuts, also the bad starting cuts (those with the first epsilon already 1 (put to 0), which cuts away all the bkg). We cannot give a null cut
                    print("boost cut modality.")
                    if ( (epsilon_M == 0) & (epsilon == 0) ):
                        pins[index]=[1,0]
                        print("boost mod, set: ",pins[index])
                if self.config["GENERAL"]["mod2"] == "online": # see if there's a mod called 'online'. In this case we cannot perform ALL the possible cuts, nor even boost each cut too much, bc it has to be fast. So ends it if counter > x
                    print("mod2 = online")
                    online_counter = 0
                    cut_boosted = []
                print("pin variable: ",pins[index])

                while pins[index] != [0,0]: # optimize the same variable both up and down if improvements
                    print("entered while cycle")

                    # check if the current values is a range extreme
                    if ranges[index].index(cut_labs[0][index]) == 0: # lower extreme
                        pins[index][0] = 0
                        print("lower extreme: ",cut_labs[0][index],ranges[index],pins[index])
                    elif ranges[index].index(cut_labs[0][index]) == len(ranges[index])-1: # upper extreme
                        pins[index][1] = 0
                        print("upper extreme: ",cut_labs[0][index],ranges[index],pins[index])
                    print("extreme check done.")
                    # check if "direction", "new variable", or if it's time to change variable
                    print("pins per variabile: ",index,pins[index])

                    # check if counter is high, in the case of boost condition
                    if self.config["GENERAL"]["mod2"] == "online":
                        if online_counter > 5:
                            pins[index][u] = 0
                            u = 8 # stop it, it's enough, we can't lose more time
                            print("I've tried to boost this cut, but it didn't improve in 5 steps: ",cut_boosted)
                    print("counter checked/not checked")

                    if pins[index] == [1,1]:
                        u = random.randint(0,1)
                        print("normal variable, choose randomly ",u)
                    elif pins[index] == [1,0]:
                        u = 0
                        print("pins [1,0] so improve down, u = ",u)
                    elif pins[index] == [0,1]:
                        u = 1
                        print("pins [0,1] so improve up, u = ",u)
                    else: # [0,0] so it's not a good variable anymore, change var
                        u = 8 # reset u
                        print("stopped to improve. Reset: ",u)
                    print("pins checked.")
                    
                    if u != 8:
                        momentarycut = cut_labs[0][:] # copy the old cut
                        print("copy cut: ",momentarycut,cut_labs[0])
                        if u == 0:
                            print("cut direction: ",u)
                            momentarycut[index] = ranges[index][ranges[index].index(cut_labs[0][index])-1]
                            print("substitute lower value: ",momentarycut,ranges[index],ranges[index][ranges[index].index(cut_labs[0][index])-1])
                        elif u == 1:
                            print("cut direction: ",u)
                            momentarycut[index] = ranges[index][ranges[index].index(cut_labs[0][index])+1]
                            print("substitute upper value: ",momentarycut,ranges[index],ranges[index][ranges[index].index(cut_labs[0][index])+1])
                        # apply cut
                        print("dataframe size before cut: ",dataframeee.shape[0])
                        dataframeee_cut_epoch = dataframeee
                        for i in cut_labs[1]:
                            print(self.config["VARIABLES"][self.variables[i]][5],self.variables[i],momentarycut[i])
                            if self.config["VARIABLES"][self.variables[i]][6] == "log10":
                                print("log10")
                                dataframeee_epoch = dataframeee_cut_epoch.loc[directions[self.config["VARIABLES"][self.variables[i]][5]](np.log10(dataframeee_cut_epoch[self.variables[i]]),momentarycut[i])]
                                dataframeee_cut_epoch = dataframeee_epoch
                            else:
                                dataframeee_epoch = dataframeee_cut_epoch.loc[directions[self.config["VARIABLES"][self.variables[i]][5]](dataframeee_cut_epoch[self.variables[i]],momentarycut[i])]
                                dataframeee_cut_epoch = dataframeee_epoch
                            # fare <> e fare log
                        print("dataframe size post cut: ",dataframeee_epoch.shape[0])
                        # compute epsilon
                        epsilon_epoch, numberbkgevents_e, bkgevents_e = self.score_function(omegaonnn_i,omegaoffff_i,dataframeee_epoch)
                        print("score function = ",epsilon_epoch, numberbkgevents_e)
                        if epsilon_epoch == 1:
                            epsilon_epoch = 0
                            print("too strong optimization.",epsilon_epoch)
                
                        # check if optimized
                        if self.epsilon_op(self.config["GENERAL"]["mod"],epsilon_epoch,epsilon):
                            if operator.eq(epsilon_epoch,epsilon): # if it didn't change however we want to try it for x steps more
                                print("epsilon didn't change: ",epsilon_epoch,epsilon)
                                if self.config["GENERAL"]["mod2"] == "online": # already 1 trial, note it
                                    print("mod2 = online (increase counter)")
                                    cut_boosted.append(cut_labs[0])
                                    print("added one boosted cut ",cut_boosted)
                                    online_counter += 1
                            else:
                                print("improved: ",epsilon_epoch,epsilon)
                            epsilon = epsilon_epoch # update the score function
                            cut_labs[0] = momentarycut[:] # update the best cut (for this starting point, millennium)
                            numberbkgevents = numberbkgevents_e
                            bkgevents = bkgevents_e
                            #numbersgnevents = numbersgnevents_e
                            #sgnevents = sgnevents_e
                            # check improving direction
                            if int(self.config["GENERAL"]["linearity"]) == 1:
                                if u == 0:
                                    print("ok, keep lowering this variable")
                                    pins[index] = [1,0]
                                elif u == 1:
                                    print("ok, keep raising the variable")
                                    pins[index] = [0,1]
                                else:
                                    print("u cannot be here!")
                                    enablePrint()
                                    return null_res, "Segmentation_fault: linearity == 1 but u != 0 or 1: "+str(u), numberbkgeventspre, bkgeventspre, numberbkgevents, bkgevents, epsilon, 0
                            else: # non-linear
                                print("Apologies human =( \n but this kind of linearity has not been developed yet. \n If you want to give the model another linearity (different than 1) you have to develop the corresponding code.")
                                enablePrint()
                                return null_res, "Segmentation_fault: linearity != 1 hasn't been developed yet: "+str(int(self.config["GENERAL"]["linearity"])), numberbkgeventspre, bkgeventspre, numberbkgevents, bkgevents, epsilon, 0

                            # if condition modality's on: check if it's been satisfied
                            if self.config["GENERAL"]["modality"] ==  "condition":
                                print("condition modality on. Let's check if it is satisfied")
                                if epsilon_epoch > epsiloncondition: 
                                    print("satisfied.")
                                    enablePrint()
                                    return cut_labs[0], "YES", numberbkgeventspre, bkgeventspre, numberbkgevents, bkgevents, epsilon, epsiloncondition
                                #else:
                                #    print("no condition.")
                            elif ( (self.config["GENERAL"]["modality"] !=  "no_condition") & (self.config["GENERAL"]["modality"] !=  "condition") ):
                                print("Apologies human =( \n but this modality has not been developed yet. \n If you want to teach the model another modality (different than 'condition') you have to develop the corresponding code.")
                                enablePrint()
                                return null_res, "Segmentation_fault: modality different than 'no_condition' or 'condition' hasn't been developed yet: "+self.config["GENERAL"]["modality"], numberbkgeventspre, bkgeventspre, numberbkgevents, bkgevents, epsilon, 0
                        else:
                            pins[index][u] = 0
                            print("This technique didn't work: pins = ",pins[index])

                    else:
                        print("It doesn't improve anymore, so let's change strategy.")
            print("\n\n")
            if epsilon > epsilon_best: # it improved, optimized
                c_best = cut_labs[:]
                epsilon_best = epsilon
                nbkg_best = numberbkgevents
                bkg_best = bkgevents
                #nsgn_best = numbersgnevents
                #sgn_best = sgnevents
                
        print("best cut ever: ",c_best)
        if c_best:
            if self.config["GENERAL"]["modality"] ==  "condition":
                print("Best found. condition on: save results even if it didn't reach it!")
                enablePrint()
                return c_best[0], "no", numberbkgeventspre,bkgeventspre,nbkg_best,bkg_best, epsilon_best, epsiloncondition
            else:
                print("Best found. condition off: save results.")
                enablePrint()
                return c_best[0], "no condition applied", numberbkgeventspre,bkgeventspre,nbkg_best,bkg_best, epsilon_best, 0
        else:
            if self.config["GENERAL"]["modality"] ==  "condition":
                print("No best found. condition on. save results even if it didn't reach it!")
                print(self.variables)
                print(null_res)
                print(numberbkgeventspre)
                print(bkgeventspre)
                print(epsiloncondition)
                enablePrint()
                return null_res, "no best cut has been found. Consider to change initial conditions.", numberbkgeventspre, bkgeventspre, 0, 0, 0, epsiloncondition
            else:
                print("No best found. condition off. save results.")
                enablePrint()
                return null_res, "no best cut has been found. Consider to change initial conditions.", numberbkgeventspre,bkgeventspre,0,0,0, 0


    def opt_results(self,results_selection,condition,numb_bkg_pre,val_bkg_pre,numb_bkg,val_bkg,scf,scf_condition,mod_print):
        #check that if opt_selection failed we can't do this!
        print("saving results")
        if mod_print == "debug":
            enablePrint()
        else:
            blockPrint()
        report_text = ["GENERAL INFO","\n","Alert name: "+str(self.ID),"\n","Alert date: "+str(self.data_str),"\n\n","SELECTION INFO","\n"]
        print("report_text titoli")
        #variables = list(self.config["VARIABLES"]) # optimizing variables
        var_ensemble = ''
        var_dir = ''
        safe_ranges = ''
        good_ranges = ''
        steps = ''
        prob_c = 0
        prob_s = 0
        print("reading variables info")
        for v in self.variables:
            var_ensemble += v+" "
            print("var_ensemble")
            safe_ranges += str(self.config["VARIABLES"][v][1])+":"+str(self.config["VARIABLES"][v][2])+" "
            print("safe")
            good_ranges += str(self.config["VARIABLES"][v][0])+":"+str(self.config["VARIABLES"][v][3])+" "
            print("ranges")
            steps += str(self.config["VARIABLES"][v][4])+" "
            print("steps")
            var_dir += self.config["VARIABLES"][v][5]+" "
            print("var_dir")
        print("got variables info")
        prob_c = str(self.stpoint)
        print("proba_c")
        prob_s = str(self.ord)
        print("proba_s")
        report_text.append("Optimization variables: "+var_ensemble)
        report_text.append("\n")
        report_text.append("Starting ranges: "+safe_ranges)
        report_text.append("\n")
        report_text.append("Optimization ranges: "+good_ranges)
        report_text.append("\n")
        report_text.append("Variables' granularity: "+steps)
        report_text.append("\n")
        report_text.append("Cut directions: "+var_dir)
        report_text.append("\n")
        report_text.append("livetime ON [s]: "+str(self.livetime_sig)+"\t livetime OFF [s]: "+str(self.livetime_bkg))
        report_text.append("\n")
        report_text.append("Starting point probability: "+prob_c)
        report_text.append("\n")
        report_text.append("Optimization order probability: "+prob_s)
        report_text.append("\n")
        # P(best_cut) = P(best order | best rnd) * P(find best cut) = P(best rnd | best order) * P(find best cut)
        # P(find best cut) = granularity(var1)*gran(var2)*...
        print("saved variables info")
        if self.config["GENERAL"]["modality"] == "condition":
            if self.config["GENERAL"]["stop"] == 0.9973072703700011:
                report_text.append("Algorithm modality: 3 sigma significance")
        else:
            report_text.append("Algorithm modality: ",self.config["GENERAL"]["modality"])
        print("modality saved")
        report_text.append("\n")
        report_text.append("\n")
        bkg_tot = 0.
        nbkg_tot = 0
        for b in range(self.nbands):
            report_text.append("BAND: "+str(b))
            report_text.append("\n")
            print("BAND: "+str(b))
            report_text.append("solid angle ON: "+str(self.omegaon[b])+"\t solid angle OFF: "+str(self.omegaoff[b]))
            report_text.append("\n")
            print("solid angle ON: "+str(self.omegaon[b])+"\t solid angle OFF: "+str(self.omegaoff[b]))
            report_text.append("number of OFF events pre-cut: "+str(numb_bkg_pre[b])+"\t expected background pre-cut: "+str(val_bkg_pre[b]))
            report_text.append("\n")
            print("number of OFF events pre-cut: "+str(numb_bkg_pre[b])+"\t expected background pre-cut: "+str(val_bkg_pre[b]))
            report_text.append("best cut: "+str(results_selection[b]))
            report_text.append("\n")
            print("best cut: "+str(results_selection[b]))
            print("bkg", numb_bkg[b], val_bkg[b])
            bkg_tot += val_bkg[b]
            nbkg_tot += numb_bkg[b]
            print("bkg summed ",bkg_tot)
            report_text.append("number of OFF events: "+str(numb_bkg[b])+"\t expected background: "+str(val_bkg[b]))
            report_text.append("\n")
            print("number of OFF events: "+str(numb_bkg[b])+"\t expected background: "+str(val_bkg[b]))
            # score f
            report_text.append("score function: "+str(scf[b]))
            report_text.append("\n")
            print("score function: "+str(scf[b]))
            # condition score f
            print("condition score function: "+str(scf_condition[b]))
            report_text.append("condition score function: "+str(scf_condition[b]))
            report_text.append("\n")
            print("Condition satisfied: ",condition[b])
            report_text.append("Condition satisfied: "+str(condition[b]))
            report_text.append("\n")
            report_text.append("--------------------------------------------------------------------------")
            report_text.append("\n")
            print("prepare new band")

        report_text.append("\n")
        report_text.append("Total expected background ("+str(nbkg_tot)+" events): "+str(bkg_tot))
        report_text.append("\n")
        print("Total expected background ("+str(nbkg_tot)+" events): "+str(bkg_tot))

        # write it into a file
        now = datetime.now()
        a = self.outdir+'/'+self.ID+"/"+'optimization_'+str(now.strftime("%d-%m-%Y_%H-%M"))+".txt"
        print("producing txt file in ",a)
        with open(self.outdir+'/'+self.ID+"/"+'optimization_'+str(now.strftime("%d-%m-%Y_%H-%M"))+".txt", 'w') as fw:
            for item in report_text:
                fw.writelines(item)
        enablePrint() # reactivate prints




    def optimization_procedure(self):
        self.total_omegaon = np.sum(self.omegaon)
        self.total_omegaoff = np.sum(self.omegaoff)
        # optimization
        time_now = time.time()
        results_selection = []
        condition = []
        numbbkgpre = []
        valbkgpre = []
        numbbkg = []
        valbkg = []
        scf = []
        scf_c = []
        self.mask_OFF2 = (self.dataframe['mjd'] >= Time(self.timebkg_start,format='isot').to_value('mjd')) & (self.dataframe['mjd'] <= Time(self.timebkg_stop,format='isot').to_value('mjd'))
        for iband in range(self.nbands):
            df_band = self.dataframe[self.mask_OFF2 & self.dataframe["pixel_local"].isin(self.pixelbands[iband])]
            #ratio_livetime = livetime_bkg / livetime_sig
            #ratio_omega = bands["omegaoff"][iband] / total_omegaon
            #bkg_threshold = ratio_omega * ratio_livetime * config["OPTIMIZATION"]["BKG_THRESHOLD"] # = MAX acceptable NUMB. OF EVENTS to have 3sigma overall
            # bkg = #events * omegaON_i/omegaOFF_i * livetimeON/livetimeOFF
            bkg_threshold_i = (self.omegaon[iband]/self.total_omegaon)*((1/self.config["GENERAL"]["stop"])-1)
            cuts, cond_i, nbkgp_i, bkgp_i, nbkg_i, bkg_i, scf_i, scf_c_i = self.opt_selection(bkg_threshold_i,self.omegaon[iband],self.omegaoff[iband],df_band,'no')
            results_selection.append(cuts)
            condition.append(cond_i)
            numbbkgpre.append(nbkgp_i)
            valbkgpre.append(bkgp_i)
            numbbkg.append(nbkg_i)
            valbkg.append(bkg_i)
            scf.append(scf_i)
            scf_c.append(scf_c_i)

        time_later = time.time()
        print("Optimization time [s]: ",time_later-time_now)

        #save
        time_now = time.time()
        self.opt_results(results_selection,condition,numbbkgpre,valbkgpre,numbbkg,valbkg,scf,scf_c,"no")
        time_later = time.time()
        print("Saving time [s]: ",time_later-time_now)

        # apply cuts
        self.dataframe = self.dataframe[self.mask_OFF2]
        finalcut = np.zeros(len(self.dataframe), dtype=bool)
        directions = {">" : operator.gt, "<" : operator.lt}
        for iband in range(self.nbands):
            band = self.dataframe["pixel_local"].isin(self.pixelbands[iband])
            for j in range(len(results_selection[iband])):
                if self.config["VARIABLES"][self.variables[j]][6] == "log10":
                    finalcut[band] = np.log10(self.dataframe[band][self.variables[j]]) > results_selection[iband][j]
                else:
                    finalcut[band] = self.dataframe[band][self.variables[j]] > results_selection[iband][j]
        self.dataframe_postcut = self.dataframe[finalcut]

        # plot result
        self.mask_allON2 = (
            (self.dataframe_postcut["mjd"] >= Time(self.timesig_start,format='isot').to_value('mjd'))
            & (self.dataframe_postcut["mjd"] <= Time(self.timesig_stop,format='isot').to_value('mjd'))
        )
        usedvars = ""
        for i in self.variables:
            usedvars += i+"_"
        self.outfig_skymap_postcut = os.path.join(self.res_fold, f"SkymapOptimization_{self.ID}_variables_{usedvars}_PS_{self.stpoint}_orders_{self.ord}_iter_{self.iteration}.png")
        plot_skymap_healpix(self.skymap, self.region90ext, self.dataframe_postcut[self.mask_allON2], self.evt_date, self.outfig_skymap_postcut, self.ID) # print and save img

        # update labels
        self.mask_ON2 = (
            (self.dataframe_postcut["mjd"] >= Time(self.timesig_start,format='isot').to_value('mjd'))
            & (self.dataframe_postcut["mjd"] <= Time(self.timesig_stop,format='isot').to_value('mjd'))
            & (self.dataframe_postcut["pixel_eq"].isin(self.region90ext))
        )
        n_on = len(self.dataframe_postcut[self.mask_ON2])
        exp_bkg = np.sum(valbkg)
        condition = exp_bkg < 2.7e-3
        return n_on, exp_bkg, self.outfig_skymap_postcut


                
        
            
            
    