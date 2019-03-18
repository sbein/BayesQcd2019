#Welcome to the industrial age of Sam's rebalance and smear code. You're going to have a lot of fun!
import os,sys
from ROOT import *
from array import array
from glob import glob
from utils import *
#from ra2blibs import *
import time

###stuff that would be nice in a config file
mhtjetetacut = 5.0 # also needs be be changed in UsefulJet.h
lhdHardMetJetPtCut = 15.0
AnHardMetJetPtCut = 30.0
cutoff = 15.0
nCuts = 8
isdata = False
rebalancedMetCut = 100

##load in delphes libraries to access input collections:
gSystem.Load("/nfs/dust/cms/user/beinsam/RebalanceAndSmear/CMSSW_10_1_0/src/SampleProduction/delphes/build/libDelphes")
##load in UsefulJet class, which the rebalance and smear code uses
gROOT.ProcessLine(open('src/UsefulJet.cc').read())
exec('from ROOT import *')
gROOT.ProcessLine(open('src/BayesRandS.cc').read())
exec('from ROOT import *')


##read in command line arguments
defaultInfile_ = "../SampleProduction/delphes/delphes_gjet4.root"
#T2qqGG.root
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", type=int, default=1,help="analyzer script to batch")
parser.add_argument("-printevery", "--printevery", type=int, default=1000,help="short run")
parser.add_argument("-fin", "--fnamekeyword", type=str,default=defaultInfile_,help="file")
parser.add_argument("-bootstrap", "--bootstrap", type=str, default='0',help="boot strapping (0,1of5,2of5,3of5,...)")
parser.add_argument("-quickrun", "--quickrun", type=bool, default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
bootstrap = args.bootstrap
verbosity = args.verbosity
printevery = args.printevery
quickrun = args.quickrun
if quickrun: n2process = 1000000
else: n2process = 9999999999999
llhdHardMetThresh = 15
mktree = False
BTag_Cut = 0.5 #Delphes b-tagging binary
if bootstrap=='0': 
	bootstrapmode = False
	bootupfactor = 1
else: 
	bootstrapmode = True
	from random import randint
	thisbootstrap, nbootstraps = bootstrap.split('of')
	thisbootstrap, nbootstraps = int(thisbootstrap), int(nbootstraps)
	print 'thisbootstrap, nbootstraps', thisbootstrap, nbootstraps
	bootupfactor = nbootstraps


#Dictionary list of signal regions
regionCuts = {}
pi = 3.14159
Inf = 9999
#varlist =                            ['St',    'HardMet','NJets','BTags','DPhiPhoPho','NPhotons']
regionCuts['NoCuts']                = [[0,Inf],  [0,Inf], [0,Inf],[0,Inf],[-Inf,Inf],    [1,Inf]]
regionCuts['Baseline1Pho']          = [[200,Inf],[70,Inf],[0,Inf],[0,Inf],[-Inf,Inf],    [1,1]]
regionCuts['Baseline2Pho']          = [[200,Inf],[70,Inf],[0,Inf],[0,Inf],[-Inf,Inf],    [2,Inf]]


##declare and load a tree
c = TChain('Delphes')
for fname in inputFiles: c.Add(fname)
nentries = c.GetEntries()
c.Show(0)
n2process = min(n2process, nentries)
print 'n(entries) =', n2process

##feed the tree to delphes, set up which branches need to be used
treeReader = ExRootTreeReader(c)
numberOfEntries = treeReader.GetEntries()
branchHT = treeReader.UseBranch("ScalarHT")
branchJets = treeReader.UseBranch("Jet")
branchGenJet = treeReader.UseBranch("GenJet")
branchMissingET = treeReader.UseBranch("MissingET")
branchGenMissingET = treeReader.UseBranch("GenMissingET")
branchPhoton = treeReader.UseBranch("Photon")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")


varlist = ['St','HardMet','NJets','BTags', 'DPhiPhoPho', 'NPhotons']
indexVar = {}
for ivar, var in enumerate(varlist): indexVar[var] = ivar
indexVar[''] = -1
nmain = len(varlist)

def selectionFeatureVector(fvector, regionkey='', omitcuts='', omitcuts_dphi=''):
	iomits, iomits_dphi = [], []  
	for cut in omitcuts.split('Vs'): iomits.append(indexVar[cut])
	for i, feature in enumerate(fvector):
		if i==nmain: break
		if i in iomits: continue
		if not (feature>=regionCuts[regionkey][i][0] and feature<=regionCuts[regionkey][i][1]): 
			return False
	return True


templateFileName = 'usefulthings/llhd-prior-coarse.root'
ftemplate = TFile(templateFileName)
print 'using templates from',templateFileName
hPtTemplate = ftemplate.Get('hPtTemplate')
templatePtAxis = hPtTemplate.GetXaxis()
hEtaTemplate = ftemplate.Get('hEtaTemplate')
templateEtaAxis = hEtaTemplate.GetXaxis()
hHtTemplate = ftemplate.Get('hHtTemplate')
templateStAxis = hHtTemplate.GetXaxis()

##Create output file
infileID = fnamekeyword.split('/')[-1].replace('.root','')
newname = 'posterior-'+infileID+'.root'
fnew = TFile(newname, 'recreate')
print 'creating', newname
if mktree:
	treefile = TFile('littletreeLowHardMet'+fnamekeyword+'.root','recreate')
	littletree = TTree('littletree','littletree')
	prepareLittleTree(littletree)

hSt = TH1F('hSt','hSt',120,0,2500)
hSt.Sumw2()
hStWeighted = TH1F('hStWeighted','hStWeighted',120,0,2500)
hStWeighted.Sumw2()
hPassFit = TH1F('hPassFit','hPassFit',5,0,5)
hPassFit.Sumw2()
hTotFit = TH1F('hTotFit','hTotFit',5,0,5)
hTotFit.Sumw2()

GleanTemplatesFromFile(ftemplate)

histoStructDict = {}
for region in regionCuts:
	for var in varlist:
		histname = region+'_'+var
		histoStructDict[histname] = mkHistoStruct(histname, binning)

xsec_times_lumi_over_nevents = 1.0
t0 = time.time()
for ientry in range(n2process):


	if ientry%printevery==0:
		print "processing event", ientry, '/', n2process, 'time', time.time()-t0
	c.GetEntry(ientry)

	weight = xsec_times_lumi_over_nevents
	
	acme_objects = vector('TLorentzVector')()
	recophotons = vector('TLorentzVector')()
	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	if not len(branchPhoton)>0: continue
	for ipho, pho in enumerate(branchPhoton):
		if not pho.PT>30: continue
		if not abs(pho.Eta)<2.4: continue
		tlvpho = TLorentzVector()
		tlvpho.SetPtEtaPhiE(pho.PT, pho.Eta, pho.Phi, pho.E)
		recophotons.push_back(tlvpho)	
		acme_objects.push_back(tlvpho)	
		
	##some universal selection (filters, etc)
	# <insert any event filters here>
	if not int(recophotons.size())>0: continue
	
	if len(recophotons)>1: dphi = abs(recophotons[0].DeltaPhi(recophotons[1]))
	else: dphi = -1
		
	recoelectrons = vector('TLorentzVector')()
	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	for iel, el in enumerate(branchElectron):
		if not el.PT>30: continue
		if not abs(el.Eta)<2.4: continue
		tlvel = TLorentzVector()
		tlvel.SetPtEtaPhiE(el.PT, el.Eta, el.Phi, el.PT*TMath.CosH(el.Eta))
		recoelectrons.push_back(tlvel)
		acme_objects.push_back(tlvel)
		
	recomuons = vector('TLorentzVector')()
	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	for imu, mu in enumerate(branchMuon):
		if not mu.PT>30: continue
		if not abs(mu.Eta)<2.4: continue
		tlvmu = TLorentzVector()
		tlvmu.SetPtEtaPhiE(mu.PT, mu.Eta, mu.Phi, mu.PT*TMath.CosH(mu.Eta))
		recomuons.push_back(tlvmu)	
		acme_objects.push_back(tlvmu)				
	
	AcmeVector = TLorentzVector()
	AcmeVector.SetPxPyPzE(0,0,0,0)
	for obj in acme_objects: AcmeVector+=obj		
	
	_Templates_.AcmeVector = AcmeVector	
	
	##declare empty vector of UsefulJets (in c++, std::vector<UsefulJet>):
	recojets = vector('UsefulJet')()

	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	for ijet, jet in enumerate(branchJets):
		if not jet.PT>15: continue
		if not abs(jet.Eta)<5: continue
		tlvjet = TLorentzVector()
		tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))
		ujet = UsefulJet(tlvjet, jet.BTag)
		if calcMinDr(acme_objects, ujet, 0.3)<0.3: continue
		recojets.push_back(ujet)

	genjets_ = vector('TLorentzVector')()
	#build up the vector of jets using TLorentzVectors; 
	#this is where you have to interface with the input format you're using
	for ijet, jet in enumerate(branchGenJet):
		if not jet.PT>15: continue
		if not abs(jet.Eta)<5: continue
		tlvjet = TLorentzVector()
		tlvjet.SetPtEtaPhiE(jet.PT, jet.Eta, jet.Phi, jet.PT*TMath.CosH(jet.Eta))
		if calcMinDr(acme_objects, tlvjet, 0.3)<0.3: continue		
		genjets_.push_back(tlvjet)
	gSt = getHt(genjets_,AnHardMetJetPtCut)
	for obj in acme_objects: gSt+=obj.Pt()
	fillth1(hSt, gSt,1)
	
	
	matchedCsvVec = createMatchedCsvVector(genjets_, recojets)
	genjets = CreateUsefulJetVector(genjets_, matchedCsvVec)

	##a few global objects
	MetVec = TLorentzVector()
	MetVec.SetPtEtaPhiE(branchMissingET[0].MET,0,branchMissingET[0].Phi,branchMissingET[0].MET)

	##observed histogram
	tSt = getHt(recojets,AnHardMetJetPtCut)
	for obj in acme_objects: tSt+=obj.Pt()
	tHardMetVec = getHardMet(recojets,AnHardMetJetPtCut, mhtjetetacut)
	tHardMetVec-=AcmeVector
	tHardMetPt, tHardMetPhi = tHardMetVec.Pt(), tHardMetVec.Phi()
	tNJets = countJets(recojets,AnHardMetJetPtCut)
	tBTags = countBJets(recojets,AnHardMetJetPtCut)
	redoneMET = redoMET(MetVec, recojets, recojets)
	tMetPt,tMetPhi = redoneMET.Pt(), redoneMET.Phi()

	#mindphi = 4
	#for jet in recojets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(tHardMetVec)))	
	fv = [tSt,tHardMetPt,tNJets,tBTags,dphi, int(recophotons.size())]

	for regionkey in regionCuts:
		for ivar, varname in enumerate(varlist):
			hname = regionkey+'_'+varname
			if selectionFeatureVector(fv,regionkey,varname,''): 
				fillth1(histoStructDict[hname].Observed, fv[ivar], weight)

	if mktree:
		if fv[1]>=150 and fv[1]<=200 and fv[0]>500 and fv[2]>3:
			growTree(littletree, fv, jetPhis, weight)            
		

	fitsucceed = RebalanceJets(recojets)
	rebalancedJets = _Templates_.dynamicJets

	mSt = getHt(rebalancedJets,AnHardMetJetPtCut)
	for obj in acme_objects: mSt+=obj.Pt()
	mHardMetVec = getHardMet(rebalancedJets,AnHardMetJetPtCut, mhtjetetacut)
	mHardMetVec-=AcmeVector
	mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()

	mNJets = countJets(rebalancedJets,AnHardMetJetPtCut)
	mBTags = countBJets(rebalancedJets,AnHardMetJetPtCut)###

	hope = (fitsucceed and mHardMetPt<rebalancedMetCut)# mHardMetPt>min(mSt/2,180):# was 160

	redoneMET = redoMET(MetVec,recojets,rebalancedJets)
	mMetPt,mMetPhi = redoneMET.Pt(), redoneMET.Phi()
	#mindphi = 4
	#for jet in rebalancedJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))	
	fv = [mSt,mHardMetPt,mNJets,mBTags,dphi, int(recophotons.size())]

	for regionkey in regionCuts:      
		for ivar, varname in enumerate(varlist):
			hname = regionkey+'_'+varname
			if selectionFeatureVector(fv,regionkey,varname,''): 
				fillth1(histoStructDict[hname].Rebalanced, fv[ivar],weight)

	fillth1(hTotFit, fv[3], weight)
 
	if hope: fillth1(hPassFit, fv[3], weight)

	nsmears = 3*bootupfactor
	weight = xsec_times_lumi_over_nevents / nsmears

	for i in range(nsmears):

		if not hope: break
		RplusSJets = smearJets(rebalancedJets,99+_Templates_.nparams)
		rpsSt = getHt(RplusSJets,AnHardMetJetPtCut)
		for obj in acme_objects: rpsSt+=obj.Pt()
		rHardMetVec = getHardMet(RplusSJets,AnHardMetJetPtCut, mhtjetetacut)
		rHardMetVec-=AcmeVector
		rpsHardMet, rpsHardMetPhi = rHardMetVec.Pt(), rHardMetVec.Phi()
		rpsNJets = countJets(RplusSJets,AnHardMetJetPtCut)
		rpsBTags = countBJets(RplusSJets,AnHardMetJetPtCut)
		rpsHardMetVec = mkmet(rpsHardMet,rpsHardMetPhi)
		redoneMET = redoMET(MetVec, recojets, RplusSJets)
		rpsMetPt, rpsMetPhi = redoneMET.Pt(), redoneMET.Phi()
		
		#mindphi = 4
		#for jet in RplusSJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(rpsHardMetVec)))	
		fv = [rpsSt,rpsHardMet,rpsNJets,rpsBTags,dphi, int(recophotons.size())]
		for regionkey in regionCuts:     
			for ivar, varname in enumerate(varlist):
				hname = regionkey+'_'+varname
				if selectionFeatureVector(fv,regionkey,varname,''):
					fillth1(histoStructDict[hname].RplusS, fv[ivar],weight)
		if mktree and 'JetHT' in physicsProcess and fv[1]>=150 and fv[1]<200 and fv[0]>=500 and fv[2]>=4: 
			growTree(littletree, fv, jetPhis, weight)


	if isdata: continue    
	genMetVec = mkmet(branchGenMissingET[0].MET, branchGenMissingET[0].Phi)
	weight = xsec_times_lumi_over_nevents



	weight = xsec_times_lumi_over_nevents
	gSt = getHt(genjets,AnHardMetJetPtCut)
	for obj in acme_objects: gSt+=obj.Pt()
	
	gHardMetVec = getHardMet(genjets,AnHardMetJetPtCut, mhtjetetacut)
	gHardMetVec-=AcmeVector
	gHardMetPt, gHardMetPhi = gHardMetVec.Pt(), gHardMetVec.Phi()
	gNJets = countJets(genjets,AnHardMetJetPtCut)
	gBTags = countBJets(genjets,AnHardMetJetPtCut)
	
	#mindphi = 4
	#for jet in genjets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(gHardMetVec)))	
	fv = [gSt,gHardMetPt,gNJets,gBTags,dphi, int(recophotons.size())]
	
	for regionkey in regionCuts:
		for ivar, varname in enumerate(varlist):
			hname = regionkey+'_'+varname
			if selectionFeatureVector(fv,regionkey,varname,''): 
				fillth1(histoStructDict[hname].Gen, fv[ivar],weight)

	#Gen-smearing
	nsmears = 3
	if isdata: weight = 1.0*prescaleweight/nsmears    
	else: weight = xsec_times_lumi_over_nevents / nsmears
	for i in range(nsmears):
		if not (gHardMetPt<rebalancedMetCut): break
		smearedJets = smearJets(genjets,9999)
		#smearedJets,csvSmearedJets = smearJets(genjets,matchedCsvVec,_Templates_.ResponseFunctions,_Templates_.hEtaTemplate,_Templates_.hPtTemplate,999)
		mSt = getHt(smearedJets,AnHardMetJetPtCut)
		for obj in acme_objects: mSt+=obj.Pt()
		mHardMetVec = getHardMet(smearedJets,AnHardMetJetPtCut, mhtjetetacut)
		mHardMetVec-=AcmeVector
		mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()
		mNJets = countJets(smearedJets,AnHardMetJetPtCut)
		mBTags = countBJets(smearedJets,AnHardMetJetPtCut)
		redoneMET = redoMET(genMetVec, genjets, smearedJets)
		mMetPt, mMetPhi = redoneMET.Pt(), redoneMET.Phi()
		#mindphi = 4
		#for jet in smearedJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))	
		fv = [mSt,mHardMetPt,mNJets,mBTags,dphi, int(recophotons.size())]

		for regionkey in regionCuts:
			for ivar, varname in enumerate(varlist):
				hname = regionkey+'_'+varname
				if selectionFeatureVector(fv,regionkey,varname,''): 
					fillth1(histoStructDict[hname].GenSmeared, fv[ivar],weight)

fnew.cd()
writeHistoStruct(histoStructDict)
hSt.Write()
hStWeighted.Write()

hPassFit.Write()
hTotFit.Write()
if mktree:
	treefile.cd()
	littletree.Write()
	treefile.Close()
print 'just created', fnew.GetName()
fnew.Close()






