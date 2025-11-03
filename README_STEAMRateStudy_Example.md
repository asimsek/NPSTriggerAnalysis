# Example: STEAM HLT Rate Study


Please refer to the [STEAM page](https://cms-hlt-steam.docs.cern.ch/Rates/RatesGeneral/) for most up-to-date instructions.


## Step 1: Setting up the environment


```sh
cmsrel CMSSW_15_1_0_patch1
cd CMSSW_15_1_0_patch1/src
cmsenv
scram build -j 4
```

```sh
git clone https://github.com/cms-steam/SteamRatesEdmWorkflow.git
cd SteamRatesEdmWorkflow/Prod/
```

 Step 2: Generating the EDM ROOT files

```sh
hltGetConfiguration /dev/CMSSW_15_1_0/GRun --full --offline \
 --no-output --data --process MYHLT --type GRun \
 --prescale 2p0E34+ZeroBias+HLTPhysics \
 --globaltag 150X_dataRun3_HLT_v1 --max-events -1  > hlt.py
```

```sh
edmConfigDump hlt.py > hlt_config.py
```



```sh
voms-proxy-init --voms cms --valid 168:00
user=asimsek; cert=$(find /tmp/x509up_u* -user "$user" -print -quit); cp "$cert" "/afs/cern.ch/user/${user:0:1}/$user/private/"
```


#### Test Locally:


> [!TIP]
> You can use `createFileList.py` script to create a list from the given datasets and run numbers.


```sh
cmsRun run_steamflow_cfg.py nEvents=1 &> log.txt &
```

#### Send to Condor:

```sh
python3 cmsCondorData.py run_steamflow_cfg.py /afs/cern.ch/work/a/asimsek/private/TriggerRates/CMSSW_15_1_0_patch1/src/ /eos/cms/store/user/asimsek/NPSTriggerRates/01Nov2025 -n 1 -q workday -p /afs/cern.ch/user/a/asimsek/private/x509up_u75207
```

```sh
./sub_total.jobb
```


> [!TIP]
> You can check your jobs with `condor_q` command.


## Step 3: Running the counting script


```sh
cd $CMSSW_BASE/src/SteamRatesEdmWorkflow/Rates/
```

> [!WARNING]
> Don't forget to change the `inputFilesDir`, `cmsswDir`, `json_file`, `maps`, and `flavour` variables inside the `config_makeCondorJobsData.py`.
> Set `maps = "nomaps"` to use the `Draw.py` script. 
> Although the jobs on condor is usually fast, they might take a bit longer than `maps = "nomaps"` on HTCondor (`flavour = "tomorrow"`)


```sh
python3 config_makeCondorJobsData.py
```


#### Test Locally:

```sh
source Jobs/Job_0/sub_0.sh
```

#### Send to Condor:

```sh
./sub_total.jobb
```



## Step 4: Misc

#### Calculate HLT_PS:

> [!NOTE]
> Find the dataset and HLT prescale values from [OMS](https://cmsoms.cern.ch/cms/triggers/prescale?cms_run=398183)
> All **8** `Dataset_EphemeralHLTPhysics*` dataset has the same PS value (`PS=8`) — one of them enough.
> Prescale of the `HLT_EphemeralPhysics_v` is `220`.


```
HLT_PS = (PS) x (HLT_EphemeralPhysics)
HLT_PS = 8 x 220
HLT_PS = 1760
```

#### Calculate Lumi, N_LS, Avg. Pileup:

```
cd $CMSSW_BASE/src/SteamRatesEdmWorkflow/Rates/
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
python3 -m pip install --user brilws
```


```sh
brilcalc lumi -u 1e33/cm2s -b "STABLE BEAMS" -i Json/json_398183.txt
```

> [!WARNING]
> The normtag (`normtag_PHYSICS.json`) for your JSON may not be available yet! 
> If available, add: `--normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json`
> use `-u /fb` instead of `-u 1e33/cm2s` when needed.



#### Avg Instantaneous lumi:

> [!NOTE]
> **Divide the `total delivered luminosity (1e33/cm2s)` by the `number of LS (nls)` in your JSON.** 
> **This is your `average input luminosity (Lumi_In)`, which you should write into the config file.**
> The final output needs to be presented in `1e34 /cm^2/s`.


```
Lumi_In = 9357 / 427  ## 1e33 /cm^2/s
Lumi_In = 2.19        ## 1e34 /cm^2/s
```


#### Deadtime:

**The average deadtime can be calculated using the ratio of `total recorded luminosity` and `total delivered luminosity`.**

```
Avg_DeadTime = (1 - (total recorded lumi/total delivered lumi))*100
Avg_DeadTime = (1 - (8929.191427781/9357.778527154))*100
Avg_DeadTime = 4.58%
```


#### Avg Pileup:

> [!IMPORTANT]
> The [recommended cross section for Run 3](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Recommended_cross_section) = `69.2 mb (69.2e-27)`.
> However, use `80mb (80e-27)` as recommended [here](https://indico.cern.ch/event/1344500/contributions/5792898/attachments/2791865/4869203/HLT%20Rate%20Studies-2024_v2.pdf#page=5) — matches with the OMS PU.
> Find `Fill Number` from [OMS Run Report](https://cmsoms.cern.ch/cms/runs/report?cms_run=398183) - `cms_run=11190`
> Find `number of bunches` [OMS Fill Report](https://cmsoms.cern.ch/cms/fills/report?cms_fill=11190) - `NBunches = 2448`
> [LHC Revolution Frequency](https://lhc-machine-outreach.web.cern.ch/collisions.htm) = `11245 Hz`


```
pileup = (Lumi_In x pp_inelastic_xsec x NBunches) / (LHC_Frequency)
pileup = ((2.19e34 cm^-2 s^-1) x (80e-27 cm^2)) / (2448 x 11245 s^-1)
pileup = 63.64
```




## Step 5: Merging and scaling

>[!WARNING]
> Don't forget to change the `lumi_in`, `lumi_target`, `hlt_ps`, and `maps` variables inside the `config_mergeOutputsData.py`.
> `Lumi_In`, `HLT_PS`, `Deadtime`, etc. calculation steps are given below.
> Set `maps = "nomaps"` to use the `Draw.py` script.


```sh
python3 config_mergeOutputsData.py
```


