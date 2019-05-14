# BayesQcd

This is the package for running rebalance and smear on Ra2/b-style ntuples while on an LPC machine. 
## Set up code in a nobackup area (modify appropriately if you forked the repo)

```
cmsrel CMSSW_10_1_0
cd CMSSW_10_1_0/src
cmsenv
git clone https://github.com/sbein/BayesQcd/
cd BayesQcd/
mkdir jobs
mkdir output
mkdir pdfs
mkdir pdfs/ClosureTests
```

### run the rebalance and smear code
I'm skipping the steps needed to create the prior and smearing templates

```
python tools/MaximizePosteriorTM.py --fnamekeyword Summer16v3.GJets_DR-0p4_HT-600 --quickrun True
```

Generate plots overlaying observed and R&S histograms

```
python tools/closurePlotter.py <output file from previous step>
```

This creates pdfs and a root file with canvases. You'll notice your histograms will likely suffer from low statistics, which is why it's good to use the batch system heavily when doing this (iteration time can be about 20 minutes to an hour). 


### submit jobs to the condor batch system:

This script defaults to submitting one job per input file. Assuming you have a valid proxy, the following script will initiate a large submission. Notice the first command below cleans up the jobs directory. It is important to do this before you submit. The script also suggests you to delete the output directory where your root files are returned so it can have a clean start. 

```
bash tools/CleanJobs.sh
python tools/submitjobs.py --analyzer tools/MaximizePosteriorTM.py --fnamekeyword Summer16v3.GJets_DR-0p4_HT --quickrun True
```
The quickrun option set to true tells the script to only run over 10,000 events per file. This argument can be removed when you're ready to max out your statistics. Output files will be put in the local output/<keyword> directory matching the specified keyword for the filename. The status of the jobs can be checked with

```
condor_q |grep <your user name>
```

Once the jobs are done, a wrapper for the hadd routine can be called which also fits a spline to each response function:
```
python tools/mergeHistosFinalizeWeights.py output/<folder name>
```

This applies the 1/n(simulated) on top of your the lumi*xsec weights, the latter of which was applied in the analyzer script. It creates a single file on which you can run the closurePlotter.py script. 

### smoothen responses and prior distributions with splines
```
python tools/articulateSplines.py ResponseTemplates2017.root "output/Fall17MiniAODv2.TTJets/*Fall17MiniAODv2.TTJets*.root"
python tools/articulateSplines.py ResponseTemplates2016.root "output/Summer16.QCD_HT/*Summer16.QCD_HT*.root"
```

The output file (its name is the middle argument) should be put into the usefulthings directory:

```
mv ResponseTemplates2017.root usefulthings/
```

## Rebalance and smear code

Check that the rebalance and smear code is pointing to the desired response template file, and then run the R&S code:

```
python tools/RebalanceAndSmear.py --fnamekeyword RunIIFall17MiniAODv2.QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8_248
python tools/closureTests.py <output of previous script>
```
