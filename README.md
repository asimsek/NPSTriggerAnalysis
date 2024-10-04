# SUS Trigger Analysis Framework

This framework is created to check the SUS HLT rates. It uses the OMSRatesNtuples from Silvio to do rate checks between various stages, and plots groups of triggers all in the same location. 

## Setup 

The most recent set of OMSRatesNtuples can be found here: 

```
/eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```

The trigger groups are defined in a json file: `triggerNames.json`. 


Run number is given to the framework and it can be found by looking at [OMS](https://cmsoms.cern.ch/cms/run_3/pp_fills_2024).
Please select the run number with the highest statistics, taking into account the different eras and interventions.
Also, Make sure select the runs which has a similar amount of deadtime.

## Rate Analysis

Once you have defined your eras, run the `rateAnalysis.py` script to generate the plots. 

```
python3 rateAnalysis.py triggerNames.json eras/eraRate.json plots/STEAM_Oct2024 /eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```

## Rate Monitoring Plots 

Similarly to produce the aggregate rate plots run the `rateMonitoring.py` script to generate the plots. 

```
python3 rateMonitoring.py triggerNames.json eras/eraMonitoring.json plots/STEAM_Oct2024 /eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/physics.root
```


If you want to filter for trigger paths in your own PAG/POG, see [SUS Trigger Review 2023](https://docs.google.com/spreadsheets/d/1bZl4qtq0FK1YO6wF73X49rLlnca6vZIuz204E47TPnk/edit?gid=1247874029#gid=1247874029) and [SUS Trigger Group Twiki Page](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSTriggerGroup).



