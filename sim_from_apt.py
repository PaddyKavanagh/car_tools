#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @author: Migo (SRON), Christophe (CEA), Paddy (DIAS)
#
"""

Tools to create a list of MIRISim simulation objects from an APT file

This is heavily based on Migo's xml to APT scripts and Christophe's
run_mirisim_from_apt.py script.

The idea is that, with a list of simulation objects, one can used a scene object
and simulation object defined elsewhere to run a full set of MIRISim simulations
as defined in an APT file.

Requires Migo's aptxml and aptxml_underthehood scripts.

"""
#from __future__ import absolute_import, division, print_function

from pdb import set_trace as stop
import os, sys
import shutil

import subprocess
import aptxml

from mirisim.config_parser import *

##########
# INPUT
inputFileName = '../car_tools_test/apt_file.aptx'
##########


########
# Default values for the sim_config
name = "Default Simulation"
rel_obsdate = 0.0
scene = "scene.ini"
#IMA section
filter = "F1130W"
readDetect = 'FULL'
ima_mode = 'FAST'
ima_exposures = 5
ima_integrations = 4
ima_frames = 10
#MRS section
disperser = 'SHORT'
detector = 'SW'
mrs_mode = 'SLOW'
mrs_exposures = 5
mrs_integrations = 4
mrs_frames = 10

# TODO implement dithering here once mirisim catch up with the APT
# PARAMETERS THAT SHOULD NOT BE DEFINED HERE
# BUT ARE ANYWAY BECAUSE MIRISIM IS OUTDATED REGARDING DITHERING
StartInd=1
NDither=2,
DitherPat="ima_recommended_dither.dat"

#~ # MRS defaults
#~ name="Default Simulation", rel_obsdate=0.0, scene="scene.ini",
            #~ POP='MRS', ConfigPath='MRS_1SHORT', Dither=True, StartInd=1,
            #~ NDither=2, DitherPat="mrs_recommended_dither.dat",
            #~ filter="F1130W", readDetect='FULL', ima_mode='FAST',
            #~ ima_exposures=5, ima_integrations=4, ima_frames=10,
            #~ disperser='SHORT', detector='SW', mrs_mode='SLOW',
            #~ mrs_exposures=5, mrs_integrations=4, mrs_frames=10)

#~ # IMA defaults
#~ name="Default Simulation", rel_obsdate=0.0,
            #~ scene="scene.ini", POP='IMA', ConfigPath='IMA_FULL',
            #~ Dither=True, StartInd=1, NDither=2,
            #~ DitherPat="ima_recommended_dither.dat",
            #~ filter="F1130W", readDetect='FULL', ima_mode='FAST',
            #~ ima_exposures=5, ima_integrations=4, ima_frames=10,
            #~ disperser='SHORT', detector='SW', mrs_mode='SLOW',
            #~ mrs_exposures=5, mrs_integrations=4, mrs_frames=10)

#~ # LRS defaults
#~ name="Default Simulation", rel_obsdate=0.0, scene="scene.ini",
            #~ POP='IMA', ConfigPath='LRS_SLIT', Dither=True, StartInd=1,
            #~ NDither=2, DitherPat="lrs_recommended_dither.dat",
            #~ filter="P750L", readDetect='FULL', ima_mode='SLOW',
            #~ ima_exposures=2, ima_integrations=3, ima_frames=5,
            #~ disperser='SHORT', detector='SW', mrs_mode='SLOW',
            #~ mrs_exposures=5, mrs_integrations=4, mrs_frames=10)



##################
# We read the apt file using Migo stuff
splitName = os.path.splitext(inputFileName)
xmlName = splitName[0]+'.xml'

if not os.path.isfile(xmlName):
    try:
        dummy=subprocess.check_output(['unzip',inputFileName])
        ## -f to force overwrite ("freshen")
    except subprocess.CalledProcessError as e:
        print e
        print e.output
        raise
    inputFileName = splitName[0]+'.xml'


xml=aptxml.aptxml(xmlName)
#################


targetDictionaries = []
for targ in xml.targets:
    #print 'New target:'
    newDict={}
    for key in targ.tags:
        #print key, targ.readValue(key)
        newDict[key]=targ.readValue(key)
    targetDictionaries.append(newDict)


obsModes = [aptxml.mrsObservationMosaic,
            aptxml.miriCoronObservationMosaic, aptxml.miriImagerObservationMosaic,
            aptxml.miriLrsObservationMosaic, aptxml.nirspecIfuObservationMosaic]
observations={}
for mode in obsModes:
    observations[mode]=[]
for obs in xml.observations:
    #print obs.readValue('label')
    obsDict = {}
    for obsType in obsModes:
        try:
            newObs=obsType(obs)
            #print '\tis of type',obsType
            for key in newObs.tags:
                #print '\t%s\t%s'%(tag,newObs.readValue(tag))
                obsDict[key]=newObs.readValue(key)
            for listClass,listElements in zip(newObs.listClasses, newObs.lists):
                assert len(listElements) == 1
                for elem in listElements:
                    for key in listClass.tags:
                        #print '\t%s\t%s'%(key,elem.readValue(key))
                        obsDict[key]=elem.readValue(key)
            observations[obsType].append(obsDict)
        except:
            #print '\tis NOT of type',obsType
            pass

#print observations
#sys.exit(0)



# For MRS observations
sim_cfgs = []
for params in observations[aptxml.mrsObservationMosaic]:

    # Common section
    observation_prefix = params['label']
    POP = 'MRS'
    if params['Dither'] == None:
        Dither = False
    else:
        Dither = True
        #TODO
        #~ StartInd = params['']
        #~ NDither = params['']
        #~ DitherPat = params['']
        # 'DitherType': '4-Point'
        # 'OptimizedFor': 'ALL'

    # We define only the section we're interested in. The other section, unused, still needs
    # some parameters that are taken from the default ones (we don't care as
    # long as they exists)
    # MRS section
    disperser = params['Wavelength']
    detector = 'BOTH'
    mrs_mode = params['ReadoutPattern']
    mrs_exposures = params['Exposures']
    mrs_integrations = params['Integrations']
    mrs_frames = params['Groups']

    if disperser == 'ALL':
        for tmp in ['SHORT', 'MEDIUM', 'LONG']:
            ConfigPath = 'MRS_1{}'.format(tmp) # By default we take for channel 1 because this information isn't available in the APT

            sim_config = SimConfig.makeSim(name=name, rel_obsdate=rel_obsdate, scene=scene, POP=POP,
                                            ConfigPath=ConfigPath, Dither=Dither, StartInd=StartInd,
                                            NDither=NDither, DitherPat=DitherPat, filter=filter,
                                            readDetect=readDetect, ima_mode=ima_mode, ima_exposures=ima_exposures,
                                            ima_integrations=ima_integrations, ima_frames=ima_frames,
                                            disperser=tmp, detector=detector, mrs_mode=mrs_mode,
                                            mrs_exposures=mrs_exposures, mrs_integrations=mrs_integrations,
                                            mrs_frames=mrs_frames)
            sim_cfgs.append(sim_config)

    else:
        ConfigPath = 'MRS_1{}'.format(disperser) # By default we take for channel 1 because this information isn't available in the APT

        sim_config = SimConfig.makeSim(name=name, rel_obsdate=rel_obsdate, scene=scene, POP=POP,
                                        ConfigPath=ConfigPath, Dither=Dither, StartInd=StartInd,
                                        NDither=NDither, DitherPat=DitherPat, filter=filter,
                                        readDetect=readDetect, ima_mode=ima_mode, ima_exposures=ima_exposures,
                                        ima_integrations=ima_integrations, ima_frames=ima_frames,
                                        disperser=disperser, detector=detector, mrs_mode=mrs_mode,
                                        mrs_exposures=mrs_exposures, mrs_integrations=mrs_integrations,
                                        mrs_frames=mrs_frames)

        sim_cfgs.append(sim_config)


#print sim_cfgs
#print len(sim_cfgs)
#sys.exit(0)
