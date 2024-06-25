# Introduction

Repository for millicharged particle search using BDTs (MicroBooNE). Includes training of blips, applying bdts, making figures, obtaining sensitivities and limits. This includes two different types of channels, blips and Wirecell.

We use `Uproot` and `pandas` to handle `ROOT` file inputs and outputs, `xgboost` for BDT training and obtaing BDT scores for new samples, and `pyhf` for obtaining sensitivities and limits.


# Workflow

The starting point is the ntuples produced in the MicroBooNE gpvms/grid, that's a whole different workflow which we will not detail here, but you start with one of two ntuples with a definite structure.

Each of the following subheadings describes a stage in processing and each of these stages will produce image outputs alongside some root file.

## Directory structure and some terminology

The repository has a number of directories in it. 

- `models`
  
  Contains `.txt` and `.json` files relating to a specific BDT **model**. Models are a specific trained BDT, they are determined by the sample used for training, the **features** used for training, and the (xgboost) configuration used for the training. Model names are usually acompanied by a **tag** which usually specifies the _date_ in which they were trained (and hopefully in my notes there would be an entry with some more details about what makes that training special). 
  
  The directories `root` and `img` which will be documented shortly, will have a subdirectory for each different model.
  
  The name of the files within should be self-explanatory. the .json file is the BDT model itself which can be loaded into a new  xgboost session and applied to a different sample. The _featmap_ text file is the list of features that were used for training and will be checked when applying this trained BDT to a new sample. They are typically used in the macro previous to obtainint the BDT scores as a compatibility check: if a sample you want to test against your trained BDT doesn't have the same variables, it will be impossible to assign a BDT score to them. Finally the _dump_ text file I'm not sure what it does but when training xgboost gives you the option to dump, so it's probably some information that could be used in the event of needing to debug.

- `root`

	Contains each stage of the processed ROOT files for each different model (with the exception of systematics which has its own subdirectory (although this should change)).

- `img`

	Contains all of the figures generated in all of the processing (except for sensitivities/limits, which have their own subdirectory), separated by model.

- `macros`

	Contains macros which will be called from one of the main programs in the root directory. The goal (for organisation's sake) is not to put anything here that you would want to execute yourself by hand, although they should be able to be run in a standalone manner as well.

- `sensitivity`

	Here is where sensitivities and limits are obtained. This part of the code works more or less independently of the others as it was its own repository at some point. It runs pyhf taking as inputs histogram files made with the main code, and makes the final parameter-space figures.
	


## Preselection

Wirecell follows a different workflow which isn't included in this repository. For blips we grab the ntuple output from the gpvms/grid and its raw variables and perform some basic preselection based on ADC, then take _pairs of blips_ and preselect based on the angle variables with angle coordinates coming from the NuMI target. With the preselected events we construct new kinematic variables for blip-pairs.

The preselection has for objective to reduce the background as for n blips total we would have n(n-1)/2 blip pairs. For default detector thresholds we have 40 blips per event on average, and we are working with 10,000 events so that gives us 7.8 million pairs per mass point. We reduce this considerably with preselection while not removing signal.

## BDT training

## Applying BDT to samples

## Obtaining sensitivity/limits


# Samples

Here we briefly describe the samples (simulation, systematics, data) from the two different channels (blips, Wirecell) we consider. We briefly mention the things we can plot or are relevant interesting from each sample. 

## Blips

### Simulation

### Detector systematics

### Data

## Wirecell

For the Wirecell samples we are using a different BDT training method at the moment, which is hosted in the manchester machines, as we use gpu training. It uses macros written by Pawel which are not part of this repository.

### Simulation

We use the millicharge generator, then reco1 and a "reco2wc" stage which outputs wirecell ntuples without some of the extra weights which are used for genie. We need to use a tarball release for the reco2 stage as wirecell is designed to associate truth in neutrino events only.

### Detector systematics

The same as the Simulation sample except following the detector variation workflow, which in practice only changes one or a couple of fhicls compared to the regular workflow. We have 11 different variations however and we require big statistics for all of them as well as event coincidence, so just generating these samples is relatively time consuming.

### Data

The data Wirecell ntuples have already been processed by the production team, we should be able to apply our BDT to them and obtain scores for each event, so no heavy processing is needed.

Although we do not have any formal data-blinding procedure, we do try not to base our analysis on the data and for this reason we have left data analysis for last. Whenever looking at actual data we use only a fraction of it, and we may include here the macros used to obtain a fraction of events from a given data ntuple.
