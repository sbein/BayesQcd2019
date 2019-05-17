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
AnHardMetJetPtCut = 25.0
cutoff = 15.0
isdata = False
rebalancedMetCut = 150
cleanrecluster = False
useMediumPho = False

debugmode = False
sayalot = False
nametag = {'Nom':'', 'Up': 'JerUp'}

##load in UsefulJet class, which the rebalance and smear code uses
gROOT.ProcessLine(open('src/UsefulJet.cc').read())
exec('from ROOT import *')
gROOT.ProcessLine(open('src/BayesRandS.cc').read())
exec('from ROOT import *')


##read in command line arguments
defaultInfile_ = "Summer16v3.GJets_DR-0p4_HT"
#T2qqGG.root
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbosity", type=int, default=1,help="analyzer script to batch")
parser.add_argument("-printevery", "--printevery", type=int, default=100,help="short run")
parser.add_argument("-fin", "--fnamekeyword", type=str,default=defaultInfile_,help="file")
parser.add_argument("-bootstrap", "--bootstrap", type=str, default='0',help="boot strapping (0,1of5,2of5,3of5,...)")
parser.add_argument("-jersf", "--JerUpDown", type=str, default='Nom',help="JER scale factor (JerNom, JerUp, ...)")
parser.add_argument("-forcetemplates", "--forcetemplates", type=str, default='',help="you can use this to override the template choice")
parser.add_argument("-quickrun", "--quickrun", type=bool, default=False,help="short run")
parser.add_argument("-debugmode", "--debugmode", type=bool, default=False,help="short run")
parser.add_argument("-sayalot", "--sayalot", type=bool, default=False,help="short run")
args = parser.parse_args()
fnamekeyword = args.fnamekeyword
inputFiles = glob(fnamekeyword)
bootstrap = args.bootstrap
JerUpDown = args.JerUpDown
forcetemplates = args.forcetemplates
verbosity = args.verbosity
printevery = args.printevery
debugmode = args.debugmode
printevery = args.printevery
quickrun = args.quickrun
sayalot = args.sayalot
if quickrun: 
    n2process = 10000
    if 'T2' in fnamekeyword: n2process = 20000
else: n2process = 9999999999999
print 'will analyze', n2process, 'events'
llhdHardMetThresh = 15
mktree = False

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


ntupleV = '17'
if 'Summer16' in fnamekeyword or 'Fall17' in fnamekeyword or 'Autumn18' in fnamekeyword:
    isdata___ = False
else: 
    isdata___ = True

is2017 = False
is2016 = False
is2018 = False

if 'Run2016' in fnamekeyword or 'Summer16' in fnamekeyword: 
    BTAG_deepCSV = 0.6324
    is2016 = True
if 'Run2017' in fnamekeyword or 'Fall17' in fnamekeyword: 
    BTAG_deepCSV = 0.4941
    is2017 = True
if 'Run2018' in fnamekeyword or 'Autumn18' in fnamekeyword: 
    BTAG_deepCSV = 0.4184#0.4941####
    is2018 = True

BTag_Cut = BTAG_deepCSV
#could try to veto events that have a non-isolated clean photon, say.

#Dictionary list of signal regions
regionCuts = {}
pi = 3.14159
Inf = 9999
#varlist =                            ['Ht',    'HardMet','NJets','BTags','MinDPhi','NPhotons', 'MetSignificance']
regionCuts['Baseline1Pho']          = [[0,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,Inf],    [1,1],    [-1,Inf]]
regionCuts['Baseline2Pho']          = [[0,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,Inf],    [2,Inf],  [-1,Inf]]
regionCuts['LdpBaseline1Pho']          = [[0,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,0.3],    [1,1],    [-1,Inf]]
regionCuts['LdpBaseline2Pho']          = [[0,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,0.3],    [2,Inf],  [-1,Inf]]
regionCuts['HighHtBaseline1Pho']          = [[200,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,Inf],    [1,1],    [-1,Inf]]
regionCuts['HighHtBaseline2Pho']          = [[200,Inf],[150,Inf],[0,Inf],[0,Inf],[0.0,Inf],    [2,Inf],  [-1,Inf]]

#################
# Load in chain #
#################
fnamefile = open('usefulthings/filelistV'+ntupleV+'.txt')
lines = fnamefile.readlines()
shortfname = fnamekeyword.strip()
print 'going to check for ', shortfname
fnamefile.close()
c = TChain('TreeMaker2/PreSelection')
filelist = []
for line in lines:
    if not shortfname in line: continue
    fname = '/eos/uscms//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV'+ntupleV+'/'+line#RelValQCD_FlatPt
    fname = fname.strip().replace('/eos/uscms/','root://cmseos.fnal.gov//')
    print 'adding', fname
    c.Add(fname)
    filelist.append(fname)
    if quickrun: break
n2process = c.GetEntries()
if quickrun: n2process = min(10000,n2process)
c.Show(0)
print 'n(entries) =', n2process


varlist = ['Ht','HardMet','NJets','BTags','MinDPhiHardMetJets','NPhotons', 'MetSignificance']
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

if ('Summer16' in fnamekeyword or 'Run2016' in fnamekeyword): 
    templateFileName = 'usefulthings/ResponseFunctionsMC16AllFilters'+nametag[JerUpDown]+'_deepCsv.root'
if ('Fall17' in fnamekeyword or 'Run2017' in fnamekeyword): 
    templateFileName = 'usefulthings/ResponseFunctionsMC17'+nametag[JerUpDown]+'_deepCsv.root'
if ('Autumn18' in fnamekeyword or 'Run2018' in fnamekeyword): 
    templateFileName = 'usefulthings/ResponseFunctionsMC18'+nametag[JerUpDown]+'_deepCsv.root'
if not forcetemplates=='': templateFileName = forcetemplates

print 'using templates from',templateFileName
ftemplate = TFile(templateFileName)
ftemplate.ls()
hPtTemplate = ftemplate.Get('hPtTemplate')
templatePtAxis = hPtTemplate.GetXaxis()
hEtaTemplate = ftemplate.Get('hEtaTemplate')
templateEtaAxis = hEtaTemplate.GetXaxis()

hHtTemplate = ftemplate.Get('hHtTemplate')
templateHtAxis = hHtTemplate.GetXaxis()
'''#option for using FullSim-based prior
priorFileName = templateFileName
priorFileName = 'usefulthings/ResponseFunctionsMC17AllFilters_deepCsv.root'
fprior = TFile(priorFileName)
hHtTemplate = fprior.Get('hHtTemplate')
templateHtAxis = hHtTemplate.GetXaxis()
'''

##Create output file
infileID = fnamekeyword.split('/')[-1].replace('.root','')
newname = 'posterior-'+infileID+'.root'
fnew = TFile(newname, 'recreate')
print 'creating', newname
if mktree:
    treefile = TFile('littletreeLowHardMet'+fnamekeyword+'.root','recreate')
    littletree = TTree('littletree','littletree')
    prepareLittleTree(littletree)

hHt = TH1F('hHt','hHt',120,0,2500)
hHt.Sumw2()
hHtWeighted = TH1F('hHtWeighted','hHtWeighted',120,0,2500)
hHtWeighted.Sumw2()
hGenMetGenHardMetRatio = TH1F('hGenMetGenHardMetRatio','hGenMetGenHardMetRatio',50,0,5)
hGenMetGenHardMetRatio.Sumw2()
hPassFit = TH1F('hPassFit','hPassFit',5,0,5)
hPassFit.Sumw2()
hTotFit = TH1F('hTotFit','hTotFit',5,0,5)
hTotFit.Sumw2()


GleanTemplatesFromFile(ftemplate)#, fprior)

histoStructDict = {}
for region in regionCuts:
    for var in varlist:
        histname = region+'_'+var
        histoStructDict[histname] = mkHistoStruct(histname, binning)

t0 = time.time()
for ientry in range(n2process):


    if ientry%printevery==0:
        print "processing event", ientry, '/', n2process, 'time', time.time()-t0

    if debugmode:
        #if not ientry>122000: continue
        if ientry in [298]: continue
        #if not ientry in [1839, 5348]: continue#,30548,49502]: continue
        a = 2
    c.GetEntry(ientry)

    weight = c.CrossSection

    #if not c.MET>100: continue
    if not passesUniversalSelection(c): continue

    acme_objects = vector('TLorentzVector')()
    recophotons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    if not len(c.Photons)>0: continue

    #idea: use HT to reference prior instead of ST

    for ipho, pho in enumerate(c.Photons):
        if not pho.Pt()>75: continue #trigger is pho 70
        if not bool(c.Photons_fullID[ipho]): continue ##need to BOOL this?
        if not abs(pho.Eta())<2.4: continue
        tlvpho = TLorentzVector()
        tlvpho.SetPtEtaPhiE(pho.Pt(), pho.Eta(), pho.Phi(), pho.E())        
        acme_objects.push_back(tlvpho)
        if useMediumPho: 
          if abs(pho.Eta())<1.48:   
            if not (c.Photons_hadTowOverEM[ipho]<0.0396  and c.Photons_sigmaIetaIeta[ipho]<0.01022 and c.Photons_pfChargedIsoRhoCorr[ipho]<0.441 and c.Photons_pfNeutralIsoRhoCorr[ipho]<(2.725+0.0148*c.Photons[ipho].Pt()+0.000017*pow(c.Photons[ipho].Pt(),2)) and c.Photons_pfGammaIsoRhoCorr[ipho]<(2.571+0.0047*c.Photons[ipho].Pt())): continue
          if abs(pho.Eta())>=1.48:
            if not (c.Photons_hadTowOverEM[ipho]<0.0219 and c.Photons_sigmaIetaIeta[ipho]<0.03001  and c.Photons_pfChargedIsoRhoCorr[ipho]< 0.442  and c.Photons_pfNeutralIsoRhoCorr[ipho]<(1.715+0.0163*c.Photons[ipho].Pt()+0.000014*pow(c.Photons[ipho].Pt(),2)) and c.Photons_pfGammaIsoRhoCorr[ipho]<(3.863+0.0034*c.Photons[ipho].Pt())): continue
        if sayalot:
            print ientry, 'acme photon', pho.Pt(), pho.Eta(), pho.Phi()
            print 'Photons_genMatched', c.Photons_genMatched[ipho]
            print 'Photons_nonPrompt', bool(c.Photons_nonPrompt[ipho])
            print 'Photons_pfGammaIsoRhoCorr', c.Photons_pfGammaIsoRhoCorr[ipho]

        #if not c.Photons_genMatched[ipho]: continue
        #if bool(c.Photons_nonPrompt[ipho]): continue
        if bool(c.Photons_hasPixelSeed[ipho]): continue                	
        recophotons.push_back(tlvpho)

    if not len(c.Photons)==recophotons.size():
        #print 'this is important'
        continue
    else: 
        a = 2
        #print 'not important'

    '''
    for ipho, pho in enumerate(c.GenParticles):
        if not pho.Pt()>100: continue #trigger is pho 70
        if not c.GenParticles_PdgId[ipho]==22: continue
        if not c.GenParticles_Htatus[ipho]==1: continue
        if not abs(pho.Eta())<2.4: continue        
        acme_objects.push_back(pho)				
        recophotons.push_back(pho)
    '''
    
    if not int(recophotons.size())>0: continue
    #print ientry, 'doing this'
    #c.Show(ientry)
    #if c.NJets>8: exit
    

    recoelectrons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for iel, el in enumerate(c.Electrons):
        if not el.Pt()>10: continue
        if not c.Electrons_mediumID[iel]: continue
        if not c.Electrons_passIso[iel]: continue        
        tlvel = TLorentzVector()
        tlvel.SetPtEtaPhiE(el.Pt(), el.Eta(), el.Phi(), el.Pt()*TMath.CosH(el.Eta()))
        if debugmode:
            print ientry, 'acme electron', el.Pt()		
        #acme_objects.push_back(tlvel)		
        if not abs(el.Eta())<2.4: continue		
        recoelectrons.push_back(tlvel)
    if not len(recoelectrons)==0: continue    

    recomuons = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using
    for imu, mu in enumerate(c.Muons):
        if not mu.Pt()>10: continue
        if not c.Muons_mediumID[imu]: continue
        if not c.Muons_passIso[imu]: continue
        tlvmu = TLorentzVector()
        tlvmu.SetPtEtaPhiE(mu.Pt(), mu.Eta(), mu.Phi(), mu.Pt()*TMath.CosH(mu.Eta()))
        if debugmode:
            print ientry, 'acme muon', mu.Pt()				
        #acme_objects.push_back(tlvmu)			
        if not abs(mu.Eta())<2.4: continue		
        recomuons.push_back(tlvmu)		
    if not len(recomuons)==0: continue        		

    AcmeVector = TLorentzVector()
    AcmeVector.SetPxPyPzE(0,0,0,0)
    for obj in acme_objects: AcmeVector+=obj		

    _Templates_.AcmeVector = AcmeVector	


    if cleanrecluster: recojets_ = CreateUsefulJetVector(c.Jetsclean, c.Jetsclean_bJetTagDeepCSVBvsAll)
    else: recojets_ = CreateUsefulJetVector(c.Jets, c.Jets_bJetTagDeepCSVBvsAll)



    recojets_.clear()
    for ijet, jet in enumerate(c.Jets):
        if not (jet.Pt()>2 and abs(jet.Eta())<5.0): continue
        recojets_.push_back(UsefulJet(jet, c.Jets_bJetTagDeepCSVBvsAll[ijet], float(int(bool(c.Jets_ID[ijet])))))
        if sayalot and jet.Pt()>AnHardMetJetPtCut:
            print 'ijet original', ijet, 'pt, eta, phi', jet.Pt(), jet.Eta(), jet.Phi(), 'Jets_bJetTagDeepCSVBvsAll', c.Jets_bJetTagDeepCSVBvsAll[ijet], 'Jets_bDiscriminatorCSV', c.Jets_bDiscriminatorCSV[ijet], 'Jets_chargedEmEnergyFraction', c.Jets_chargedEmEnergyFraction[ijet], 'Jets_chargedHadronEnergyFraction', c.Jets_chargedHadronEnergyFraction[ijet], 'Jets_chargedMultiplicity', c.Jets_chargedMultiplicity[ijet], ', Jets_multiplicity', c.Jets_multiplicity[ijet], ', Jets_partonFlavor', c.Jets_partonFlavor[ijet] , ', Jets_photonEnergyFraction', c.Jets_photonEnergyFraction[ijet] , ', Jets_photonMultiplicity', c.Jets_photonMultiplicity[ijet], 'Jets_hadronFlavor', c.Jets_hadronFlavor[ijet], 'Jets_hadronFlavor', c.Jets_hadronFlavor[ijet], 'Jets_ID', bool(c.Jets_ID[ijet])

            
    if is2017: # ecal noise treatment
        recojets.clear()
        for ijet, jet in enumerate(c.Jets):
            if not (jet.Pt()>2 and abs(jet.Eta())<5.0): continue
            if abs(jet.Eta())>2.65 and abs(jet.Eta()) < 3.139 and jet.Pt()/c.Jets_jecFactor[ijet]<50: continue #/c.Jets_jerFactor[ijet]
            recojets.push_back(UsefulJet(jet, c.Jets_bJetTagDeepCSVBvsAll[ijet], jet.Pt()))  


    ##declare empty vector of UsefulJets (in c++, std::vector<UsefulJet>):
    recojets = vector('UsefulJet')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using

    nMatchedAcmeOuterPairs = 0
    nMatchedAcmeInnerPairs = 0
    passesJetId = True
    for ijet, jet in enumerate(recojets_):
        if not jet.Pt()>15: continue
        if not abs(jet.Eta())<5: continue
        closestAcme = getClosestObject(acme_objects, jet, 0.1)
        if jet.DeltaR(closestAcme)<0.6:
            nMatchedAcmeOuterPairs+=1
            if jet.DeltaR(closestAcme)<0.2:
                nMatchedAcmeInnerPairs+=1
                if sayalot: print 'skipping reco jet with pT, eta', jet.Pt(), jet.Eta(), jet.DeltaR(acme_objects[0])
                continue
        #    #jet.tlv-=closestAcme
        if jet.Pt()>AnHardMetJetPtCut and jet.OriginalPt()<0.5: 
            passesJetId = False
            print 'failed jet id'
            break
        ujet = UsefulJet(jet.tlv, jet.btagscore)        
        recojets.push_back(ujet)

    if not passesJetId: continue
    if cleanrecluster:
        if not nMatchedAcmeJetParis==0: continue
    else:
        if not nMatchedAcmeInnerPairs==len(acme_objects): continue
        if not nMatchedAcmeOuterPairs==nMatchedAcmeInnerPairs: continue

    genjets_ = vector('TLorentzVector')()
    #build up the vector of jets using TLorentzVectors; 
    #this is where you have to interface with the input format you're using (treemaker)
    for ijet, jet in enumerate(c.GenJets):
        #if not jet.Pt()>15: continue
        #if not abs(jet.Eta())<5: continue
        tlvjet = TLorentzVector()
        tlvjet.SetPtEtaPhiE(jet.Pt(), jet.Eta(), jet.Phi(), jet.Pt()*TMath.CosH(jet.Eta()))
        closestAcme = getClosestObject(acme_objects, tlvjet, 0.1)
        if tlvjet.DeltaR(closestAcme)<0.1: 
            if sayalot: print 'ueberspringen jet with pT, eta', tlvjet.Pt(), tlvjet.Eta(), '(',tlvjet.DeltaR(acme_objects[0]),')'
            #tlvjet-=closestAcme
            continue
        genjets_.push_back(tlvjet)
    gHt = getHt(genjets_,AnHardMetJetPtCut)
    gHt = gHt
    #for obj in acme_objects: gHt+=obj.Pt()
    fillth1(hHt, gHt,1)

    matchedBtagscoreVec = createMatchedBtagscoreVector(genjets_, recojets)
    genjets = CreateUsefulJetVector(genjets_, matchedBtagscoreVec)

    ##a few global objects
    MetVec = TLorentzVector()
    MetVec.SetPtEtaPhiE(c.MET,0,c.METPhi,c.MET)

    CleanMetVec = TLorentzVector()
    CleanMetVec.SetPtEtaPhiE(c.METclean,0,c.METPhiclean,c.METclean)
    redirtiedMetVec = CleanMetVec.Clone()
    redirtiedMetVec-=AcmeVector 

    ##observed histogram
    tHt = getHt(recojets,AnHardMetJetPtCut)
    tHt = tHt
    #for obj in acme_objects: tHt+=obj.Pt()
    tHardMhtVec = getHardMet(recojets,AnHardMetJetPtCut, mhtjetetacut)
    tHardMetVec = tHardMhtVec.Clone()
    tHardMetVec-=AcmeVector
    tHardMetPt, tHardMetPhi = tHardMetVec.Pt(), tHardMetVec.Phi()


    if  not abs(MetVec.Pt()-tHardMetVec.Pt())<60:
            print ientry, 'met/mht consistency', abs(MetVec.Pt()-tHardMetVec.Pt())
            continue
    #tHardMhtPt, tHardMhtPhi = tHardMhtVec.Pt(), tHardMhtVec.Phi()


    if sayalot:
        print 'going to play with this photon'
        pho = recophotons[0]    
        print 'pho pt, eta, phi', pho.Pt(), pho.Eta(), pho.Phi()
        print 'MHT = ', c.MHT
        mht = TLorentzVector()
        mht.SetPtEtaPhiE(0,0,0,0)
        for jet in c.Jets: 
            if jet.Pt()>30 and abs(jet.Eta())<5.0: mht-=jet
        print 'MHT-homegrown', mht.Pt()
        print 'MHT-clean', c.MHTclean
        mhtagain = mkmet(c.MHTclean, c.MHTPhiclean)
        mhtagain-=pho
        print 'MHT-clean-dirtied', mhtagain.Pt()
        print 'len(c.Photons)', len(c.Photons)


    if cleanrecluster:
        mhtagain = mkmet(c.MHTclean, c.MHTPhiclean)
        mhtagain-=AcmeVector
        tHardMetPt, tHardMetPhi = mhtagain.Pt(), mhtagain.Phi()
        if c.NJets-c.NJetsclean==len(acme_objects): tHardMetPt = mhtagain.Clone()
        tHardMetPt, tHardMetPhi = tHardMetPt.Pt(), tHardMetPt.Phi()
    
    mindphi = 4
    for jet in recojets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(tHardMetVec)))


    #met_consistency = abs(c.MET-tHardMetPt)/tHardMetPt
    #if not met_consistency<0.5: 
    #    print 'throwing away this event on grounds of suspicion', tHardMhtPt
    #    continue
    #print 'accepting', tHardMhtPt

    if tHt>0: tMetSignificance = tHardMetPt/TMath.Sqrt(tHt)
    else: tMetSignificance = 99
    tNJets = countJets(recojets,AnHardMetJetPtCut)
    tBTags = countBJets(recojets,AnHardMetJetPtCut)

    if debugmode:
        if not tHardMetPt>100 and tHardMetPt<130: continue

    fv = [tHt,tHardMetPt,tNJets,tBTags,mindphi, int(recophotons.size()),tMetSignificance]

    if tHardMetPt>100: 
        print ientry, 'fv', fv
            
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
    #else:
    #	fitsucceed = True
    #	rebalancedJets = recojets
    mHt = getHt(rebalancedJets,AnHardMetJetPtCut)
    mHt = mHt
    #for obj in acme_objects: mHt+=obj.Pt()
    mHardMetVec = getHardMet(rebalancedJets,AnHardMetJetPtCut, mhtjetetacut)
    mHardMetVec-=AcmeVector
    mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()
    if mHt>0: mMetSignificance = mHardMetPt/TMath.Sqrt(mHt)
    else: mMetSignificance = 8	

    mNJets = countJets(rebalancedJets,AnHardMetJetPtCut)
    mBTags = countBJets(rebalancedJets,AnHardMetJetPtCut)###

    #hope = (fitsucceed and mHardMetPt<rebalancedMetCut)# mHardMetPt>min(mHt/2,180):# was 160
    #hope = (fitsucceed and mSignificance<rebalancedSignificanceCut)# mHardMetPt>min(mHt/2,180):# was 160
    hope = (fitsucceed and mHardMetPt<rebalancedMetCut)# mHardMetPt>min(mHt/2,180):# was 160	

    redoneMET = redoMET(MetVec,recojets,rebalancedJets)
    mMetPt,mMetPhi = redoneMET.Pt(), redoneMET.Phi()
    mindphi = 4
    for jet in rebalancedJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))
    fv = [mHt,mHardMetPt,mNJets,mBTags,mindphi, int(recophotons.size()),mMetSignificance]

    for regionkey in regionCuts:      
        for ivar, varname in enumerate(varlist):
            hname = regionkey+'_'+varname
            if selectionFeatureVector(fv,regionkey,varname,''): 
                fillth1(histoStructDict[hname].Rebalanced, fv[ivar],weight)

    fillth1(hTotFit, fv[3], weight)

    if hope: fillth1(hPassFit, fv[3], weight)

    nsmears = 100*bootupfactor
    weight = c.puWeight * c.CrossSection / nsmears

    for i in range(nsmears):

        if not hope: break
        RplusSJets = smearJets(rebalancedJets,99+_Templates_.nparams)
        rpsHt = getHt(RplusSJets,AnHardMetJetPtCut)
        rpsHt = rpsHt
        #for obj in acme_objects: rpsHt+=obj.Pt()
        rpsHardMetVec = getHardMet(RplusSJets,AnHardMetJetPtCut, mhtjetetacut)
        rpsHardMetVec-=AcmeVector
        rpsHardMetPt, rpsHardMetPhi = rpsHardMetVec.Pt(), rpsHardMetVec.Phi()
        if rpsHt>0: rpsMetSignificance = rpsHardMetPt/TMath.Sqrt(rpsHt)
        else: rpsMetSignificance = 8			
        rpsNJets = countJets(RplusSJets,AnHardMetJetPtCut)
        rpsBTags = countBJets(RplusSJets,AnHardMetJetPtCut)
        mindphi = 4
        for jet in RplusSJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(rpsHardMetVec)))
        fv = [rpsHt,rpsHardMetPt,rpsNJets,rpsBTags,mindphi, int(recophotons.size()),rpsMetSignificance]
        #print i, 'of', nsmears, fv
        for regionkey in regionCuts:     
            for ivar, varname in enumerate(varlist):
                hname = regionkey+'_'+varname
                if selectionFeatureVector(fv,regionkey,varname,''):
                    fillth1(histoStructDict[hname].RplusS, fv[ivar],weight)
        if mktree and 'JetHT' in physicsProcess and fv[1]>=150 and fv[1]<200 and fv[0]>=500 and fv[2]>=4: 
            growTree(littletree, fv, jetPhis, weight)


    if isdata: continue    
    genMetVec = mkmet(c.GenMET, c.GenMETPhi)
    weight = c.CrossSection

    gHt = getHt(genjets,AnHardMetJetPtCut)
    gHt = gHt
    #for obj in acme_objects: gHt+=obj.Pt()

    gHardMetVec = getHardMet(genjets,AnHardMetJetPtCut, mhtjetetacut)
    tmpgmet = gHardMetVec.Pt()
    gHardMetVec-=AcmeVector
    if sayalot: print 'gen mht after/before', gHardMetVec.Pt()/tmpgmet
    gHardMetPt, gHardMetPhi = gHardMetVec.Pt(), gHardMetVec.Phi()
    if gHt>0: gMetSignificance = gHardMetPt/TMath.Sqrt(gHt)
    else: gMetSignificance = 8	

    ###Delphes filter, because I think Delphes is mis-computing its own gen MET
    if gHardMetPt>80 and False: 
        fillth1(hGenMetGenHardMetRatio,abs(gHardMetPt-genMetVec.Pt())/gHardMetPt)
        if abs(gHardMetPt-genMetVec.Pt())/gHardMetPt>0.5: 
            print 'skipping', ientry, gHardMetPt, genMetVec.Pt()
            continue
    ###Delphes

    mindphi = 4
    for jet in genjets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(gHardMetVec)))
        
    gNJets = countJets(genjets,AnHardMetJetPtCut)
    gBTags = countBJets(genjets,AnHardMetJetPtCut)


    if debugmode:
        #c.Show(ientry)
        print 'missingET', c.MET
        print 'gen missingET', genMetVec.Pt()
        print 'genHardMet', gHardMetVec.Pt(), gHardMetVec.Phi()		
        print 'genjets'
        for ijet, jet in enumerate(genjets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()
            matchedreco = calcMinDr(recojets, jet, 0.01)
            print '....matched reco pt', matchedreco.Pt(), matchedreco.Eta(), matchedreco.Phi(), 'dr', matchedreco.DeltaR(jet.tlv)            
        print 'recoHardMet', tHardMetVec.Pt(), tHardMetVec.Phi()
        print 'recoHardMetUppy', (tHardMetVec+acme_objects[0]).Pt()
        print 'recoHardMetDowny', (tHardMetVec-acme_objects[0]).Pt()        
        print 'recojets:'
        for ijet, jet in enumerate(recojets):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi(), 'btag=', jet.btagscore
            matchedgen = calcMinDr(genjets, jet, 0.01)
            print '....matchedgenjet pt', matchedgen.Pt(), matchedgen.Eta(), matchedgen.Phi(), 'dr', matchedgen.DeltaR(jet.tlv)
            matchedacme = calcMinDr(acme_objects, jet, 0.01)
            print '....matchedacme dr',  matchedacme

        print 'acme_objects'
        for ijet, jet in enumerate(acme_objects):
            print ijet, jet.Pt(), jet.Eta(), jet.Phi()


    fv = [gHt,gHardMetPt,gNJets,gBTags,mindphi, int(recophotons.size()),gMetSignificance]	
    for regionkey in regionCuts:
        for ivar, varname in enumerate(varlist):
            hname = regionkey+'_'+varname
            if selectionFeatureVector(fv,regionkey,varname,''): 
                fillth1(histoStructDict[hname].Gen, fv[ivar],weight)

    #Gen-smearing
    nsmears = 3
    if isdata: weight = 1.0*prescaleweight/nsmears    
    else: weight = c.CrossSection / nsmears
    for i in range(nsmears):
        if not (gHardMetPt<rebalancedMetCut): break
        smearedJets = smearJets(genjets,9999)
        #smearedJets,csvSmearedJets = smearJets(genjets,matchedCsvVec,_Templates_.ResponseFunctions,_Templates_.hEtaTemplate,_Templates_.hPtTemplate,999)
        mHt = getHt(smearedJets,AnHardMetJetPtCut)
        mHt = mHt
        #for obj in acme_objects: mHt+=obj.Pt()
        mHardMetVec = getHardMet(smearedJets,AnHardMetJetPtCut, mhtjetetacut)
        mHardMetVec-=AcmeVector
        mHardMetPt, mHardMetPhi = mHardMetVec.Pt(), mHardMetVec.Phi()
        if mHt>0: mMetSignificance = mHardMetPt/TMath.Sqrt(mHt)
        else: mMetSignificance = 8		
        mNJets = countJets(smearedJets,AnHardMetJetPtCut)
        mBTags = countBJets(smearedJets,AnHardMetJetPtCut)
        redoneMET = redoMET(genMetVec, genjets, smearedJets)
        mMetPt, mMetPhi = redoneMET.Pt(), redoneMET.Phi()
        mindphi = 4
        for jet in smearedJets[:4]: mindphi = min(mindphi, abs(jet.DeltaPhi(mHardMetVec)))	
        fv = [mHt,mHardMetPt,mNJets,mBTags,mindphi, int(recophotons.size()),mMetSignificance]

        for regionkey in regionCuts:
            for ivar, varname in enumerate(varlist):
                hname = regionkey+'_'+varname
                if selectionFeatureVector(fv,regionkey,varname,''): 
                    fillth1(histoStructDict[hname].GenSmeared, fv[ivar],weight)

fnew.cd()
writeHistoStruct(histoStructDict)
hGenMetGenHardMetRatio.Write()
hHt.Write()
hHtWeighted.Write()

hPassFit.Write()
hTotFit.Write()
if mktree:
    treefile.cd()
    littletree.Write()
    treefile.Close()
print 'just created', fnew.GetName()
fnew.Close()



