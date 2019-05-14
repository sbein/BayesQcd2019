# BayesQcd

This is the package for running rebalance and smear on Ra2/b-style ntuples while on an LPC machine. 
## Set up code in a nobackup area (modify appropriately if you forked the repo)

```
cmsrel CMSSW_10_1_0
cd CMSSW_10_1_0/src
cmsenv
git clone https://github.com/sbein/BayesQcd/
cd BayesQcd/
```

### run a simple analyzer script

```
python tools/SimpleAnalyzer.py --fnamekeyword RunIIFall17MiniAODv2.TTJets
```

## generate jet the response templates and prior distributions

```
python tools/ResponseMaker.py --filenamekey Fall17MiniAODv2.QCD
```

### submit large batch of response functions:

If the first time:
```
mkdir output
```

Assuming you have a valid proxy, the following script will initiate a large submission 

```
python tools/submitjobs.py --analyzer tools/ResponseMaker.py --fnamekeyword Fall17MiniAODv2.QCD #(for 2017)
python tools/submitjobs.py --analyzer tools/ResponseMaker.py --fnamekeyword Summer16.QCD_HT #(for 2016)
```

Output files will be put in the local output/<keyword> directory matching the specified keyword for the filename. The status of the jobs can be checked with

```
condor_q |grep <your user name>
```

Once the jobs are done, a wrapper for the hadd routine can be called which also fits a spline to each response function:

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

if all goes well, you will now have a root file containing "closure test" plots comparing truth and method. To increase the statistics of the closure test, run over all the QCD files:

```
python tools/submitJobs.py --analyzer tools/RebalanceAndSmear.py --fnamekeyword RunIIFall17MiniAODv2.QCD_HT2000
```

The output files will be cached in a new folder created in the output directory. Hadd the files with the appropriate weights by calling:

```
python tools/mergeHistoFiles.py mergedHistos.root <path-with-wildcards-to-output-of-rands>
python tools/closureTests.py mergedHistos.root
```
