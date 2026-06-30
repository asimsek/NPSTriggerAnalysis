# NPS Trigger Analysis Framework

This framework is created to check the NPS HLT rates. It uses the `OMSRatesNtuples` from STEAM Group to do rate checks between various stages, and plots groups of triggers all in the same location.

## Setup

> [!IMPORTANT]
> The required packages for the `Jupyter Notebook` is given in the `requirements.txt` file. 

```
pip3 install -r requirements.txt
```


**For the lxplus, you need to create a new `CMSSW` (14\_X\_X or later)**

```
mkdir NPSTrigger
cd NPSTrigger/ 
cmsrel CMSSW_15_1_0_patch1
cd CMSSW_15_1_0_patch1/src
cmsenv
```

**Then, pull the framework from github:**

```
git config --global user.name <yourUserName>
git config --global user.email <yourEmailAddress>
git config --global http.emptyAuth true
unset SSH_ASKPASS
git clone --recursive https://github.com/asimsek/NPSTriggerAnalysis
cd NPSTriggerAnalysis/
```
 

**The most recent set of OMSRatesNtuples can be found here:**

```
/eos/cms/store/group/tsg/STEAM/OMSRateNtuple/2024/2024_physics_merged.root
/eos/user/s/sdonato/www/OMSRatesNtuple/OMSRatesNtuple/OMS_ntuplizer/2025_physics_merged.root
/eos/user/s/sdonato/www/OMSRatesNtuple/OMSRatesNtuple/OMS_ntuplizer/2026_physics_merged.root
```

The trigger groups are defined in a json file: `triggerNames.json`.

> [!TIP]
> Run number is given to the framework and it can be found by looking at [OMS](https://cmsoms.cern.ch/cms/run_3/pp_fills_2024) and [PdmVRun3Analysis](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Year_2024).<br>
> Please select the run number with the highest statistics, taking into account the different eras and interventions.<br>
> Also, Make sure select the runs which has a similar amount of deadtime.<br>
> For more, see [2024 Golden JSON](https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json) and [OMS HLT Trigger Rates](https://cmsoms.cern.ch/cms/triggers/hlt_trigger_rates).<br>

## Rate Analysis

Once you have defined your eras, run the `rateAnalysis.py` script to generate the plots. 

**NPS-Only Triggers:**

```
python3 rateAnalysis.py jsonFiles/triggerLists/2026/triggerNames_NPSOnly.json jsonFiles/eras/2026/eraRate.json plots/2026/NPSOnly_May2026 rootFiles/inputFiles.txt
```

## Rate Monitoring Plots

Similarly to produce the trigger rate plots run the `rateMonitoring.py` script to generate the plots. 

**NPS-Only Triggers:**

```
## 2025
python3 rateMonitoring.py jsonFiles/triggerLists/2025/triggerNames_NPSOnly.json jsonFiles/eras/2025/eraMonitoring.txt plots/2025/NPSOnly_June2025 rootFiles/inputFiles.txt --specialFills --l1seed

## 2026
python3 rateMonitoring.py jsonFiles/triggerLists/2026/triggerNames_NPSOnly.json jsonFiles/eras/2026/eraMonitoring.txt plots/2026/NPSOnly_June2026 rootFiles/inputFiles.txt --specialFills --l1seed

## Run3
python3 rateMonitoring.py jsonFiles/triggerLists/2026/triggerNames_NPSOnly.json jsonFiles/eras/Run3/eraMonitoring.txt plots/Run3/NPSOnly_June2026 rootFiles/inputFiles.txt --specialFills --l1seed
```

> [!NOTE]
> **More options for `rateMonitoring.py`:**<br>
> * `--l1seed`: Generate L1-seed trigger dictionary and run additional monitoring<br>
> * `--gRun`: Regenerate GRun.csv (requires --l1seed) before extracting seeds<br>
> * `--noHLT`: No HLT monitoring<br>
> * `--splitEraVersions`: Keep versioned eras separate in the rate-change matrix instead of merging, e.g. `2024Iv1`/`2024Iv2` rather than `2024I`.<br>
> * `--noRateMatrix`: Skip the adjacent-era rate-change matrix plot and its CSV outputs.<br>
> * `--noTrendPlots`: Skip all run/date trend PDFs.<br>
> * `--noRunPlots`: Skip run-number trend PDFs.<br>
> * `--noDatePlots`: Skip date trend PDFs.<br>
> * `--noL1TrendPlots`: (With `--l1seed`) skip only the L1 seed trend PDFs.<br>
> * `--npsRateBoxPercent`: adds a second panel with NPS / total rate (%)<br>
> * `--totalRateBranch BRANCH`: override denominator branch; default is Stream_HLTRates if available.<br>
> * `--totalRateFromAllHLT`: Reads all branches starting with "HLT_"<br>
> * `--specialFills [CONFIG]`: Build fill-number monitoring and special-fill/reference ratio plots. Without a path, uses `jsonFiles/specialFills/specialFills.json`.<br>

> [!TIP]
> Special-fill cases can define `fills`, `fillRange`, `referenceFills`, or `referenceFillRange`. If explicit reference fills are omitted, the reference is built from other selected fills in the same era, excluding all configured special fills.


## MISC

> [!TIP]
> If you require to filter/change trigger paths, find full trigger names in [NPS Trigger Review 2023](https://docs.google.com/spreadsheets/d/1bZl4qtq0FK1YO6wF73X49rLlnca6vZIuz204E47TPnk/edit?gid=1247874029#gid=1247874029) and/or [NPS Trigger Group Twiki Page](https://nps-wiki.docs.cern.ch/Trigger/) - NPS triggers.


