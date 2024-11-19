# SUS Trigger Analysis Framework

This framework is created to check the SUS HLT rates. It uses the OMSRatesNtuples from STEAM Group to do rate checks between various stages, and plots groups of triggers all in the same location. 

## Setup
The required packages for the `Jupyter Notebook` is given in the `requirements.txt` file. 

```
pip3 install -r requirements.txt
```


For the lxplus, you need to create a new `CMSSW_14_X_X`

```
mkdir SUSTrigger
cd SUSTrigger/ 
cmsrel CMSSW_14_0_4
cd CMSSW_14_0_4/src
cmsenv
```

Then, pull the framework from github:

```
scl enable rh-git29 bash
git config --global user.name <yourUserName>
git config --global user.email <yourEmailAddress>
git config --global http.emptyAuth true
unset SSH_ASKPASS
git clone --recursive https://github.com/asimsek/SUSTriggerAnalysis
cd SUSTriggerAnalysis/
```
 

The most recent set of OMSRatesNtuples can be found here: 

```
/eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```

The trigger groups are defined in a json file: `triggerNames.json`. 


Run number is given to the framework and it can be found by looking at [OMS](https://cmsoms.cern.ch/cms/run_3/pp_fills_2024) and [PdmVRun3Analysis](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Year_2024).
Please select the run number with the highest statistics, taking into account the different eras and interventions.
Also, Make sure select the runs which has a similar amount of deadtime.
For more, see [2024 Golden JSON](https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json) and [OMS HLT Trigger Rates](https://cmsoms.cern.ch/cms/triggers/hlt_trigger_rates).

## Rate Analysis

Once you have defined your eras, run the `rateAnalysis.py` script to generate the plots. 

```
python3 rateAnalysis.py triggerNames.json eras/eraRate.json plots/STEAM_Nov2024 /eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```

## Rate Monitoring Plots 

Similarly to produce the trigger rate plots run the `rateMonitoring.py` script to generate the plots. 

```
python3 rateMonitoring.py triggerNames.json eras/eraMonitoring.json plots/STEAM_Nov2024 /eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```


If you require to filter trigger paths, find full trigger names in [SUS Trigger Review 2023](https://docs.google.com/spreadsheets/d/1bZl4qtq0FK1YO6wF73X49rLlnca6vZIuz204E47TPnk/edit?gid=1247874029#gid=1247874029) and/or [SUS Trigger Group Twiki Page](https://sus-wiki.docs.cern.ch/Trigger/) - SUS triggers.



