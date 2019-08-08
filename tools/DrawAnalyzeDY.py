from ROOT import *
from utils import *
gROOT.SetBatch(1)
gStyle.SetOptStat(0)
chain = TChain('TreeMaker2/PreSelection')
chain.Add('/eos/uscms//store/user/sbein/RebalanceAndSmear/Run2ProductionV17/*Summer16v3.DYJetsToLL_M-50_HT*9*.root')
chain.Show(0)
print 'nevents =', chain.GetEntries()

universalconstraint = ' abs(MET-HardMETPt)<90 && @Muons.size()==0 && @Electrons.size()==2 && Electrons[1].Pt()>100 &&  @GenTaus.size()==0'
#universalconstraint = ' HardMETPt/MET<1.01 && (@Muons.size()==2 || @Electrons.size()==2) && @GenTaus.size()==0'

plotBundle = {}
plotBundle['Baseline_HardMet'] = ['HardMETPt>>hadc(3,130,530)','NJets>1 && abs(ZCandidates.M()-91)<1500', True]
#plotBundle['Baseline_MET'] = ['MET>>hadc(20,130,530)','NJets>1', True]
plotBundle['Baseline_NJets'] = ['NJets>>hadc(10,0,10)','HardMETPt>210', True]
plotBundle['Baseline_HT'] = ['HT>>hadc(24,0,2200)','HardMETPt>210', True]
plotBundle['Baseline_ZMass'] = ['ZCandidates.M()>>hadc(100,70,570)','HardMETPt>210', True]

plotBundle['Baseline_ZPt'] = ['ZCandidates.Pt()>>hadc(30,0,600)','HardMETPt>210', True]

fnew = TFile('ezplots_dyjets.root', 'recreate')
c1 = mkcanvas()


for key in plotBundle:
    drawarg, constraint, usernsvalue = plotBundle[key]
    obsweight = '1*('+constraint + ' && '+ universalconstraint + ' && IsUniqueSeed==1)'#puWeight * CrossSection
    print 'drawing', drawarg, ', with constraint:', obsweight


    
    chain.Draw(drawarg,obsweight)
    hobs = chain.GetHistogram().Clone(key+'_obs')
    hobs.GetYaxis().SetRangeUser(0.01,10000*hobs.GetMaximum())
    if not 'Vs' in key:
        c1.cd()
        if usernsvalue: 
            drawarg = drawarg.replace('METPt','METPtRandS').replace('NJets','NJetsRandS').replace('BTags','BTagsRandS').replace('HT','HTRandS').replace('MET>','HardMETPt>')
            randsconstraint = constraint.replace('METPt','METPtRandS').replace('NJets','NJetsRandS').replace('BTags','BTagsRandS').replace('HT','HTRandS')
        methweight = '1/NSmearsPerEvent*('+ randsconstraint + ' && '+universalconstraint+')'#puWeight * CrossSection
        print 'drawing', drawarg, ', with constraint:', methweight
        chain.Draw(drawarg, methweight)
        hrands = chain.GetHistogram().Clone(key+'_rands') 
        hrands.GetYaxis().SetRangeUser(0.01,10000*hrands.GetMaximum())
        histoStyler(hrands, kAzure-8)
        hrands.SetFillColor(hrands.GetLineColor())
        hrands.SetFillStyle(1001)

        
        leg = mklegend(x1=.45, y1=.57, x2=.95, y2=.74, color=kWhite)
        hobs.SetTitle('Summer16 QCD observed')
        hobs.GetXaxis().SetTitle(key.split('_')[-1])
        hrands.GetXaxis().SetTitle(key.split('_')[-1])
        hrands.SetTitle('QCD sim. rebalance and smeared')
        hratio, hmethodsyst = FabDrawSystyRatio(c1,leg,hobs,[hrands],datamc='MC',lumi='n/a', title = '', LinearScale=False, fractionthing='truth / method')
        c1.Update()
        c1.Write('c_'+key)        
    else:
        c2.cd()
        hobs.Draw('colz')
        c2.Update()
        hobs.Write()
        c2.Write('c_'+key)

        
'''
for key in plotBundle:
    drawarg, constraint, usernsvalue = plotBundle[key]
    obsweight = '1 * 1*('+constraint + ' && '+ universalconstraint + ' && IsUniqueSeed==1)'#puWeight
    print 'drawing', drawarg, ', with constraint:', obsweight
    chain.Draw(drawarg,obsweight)
    hobs = chain.GetHistogram().Clone(key+'_obs')
    hobs.GetYaxis().SetRangeUser(0.01,10000*hobs.GetMaximum())

    if usernsvalue: drawarg = drawarg.replace('>>','RandS>>')
    methweight = '1 * 1/NSmearsPerEvent*('+constraint.replace('>','RandS>') + ' && '+universalconstraint+')'#puWeight
    print 'drawing', drawarg, ', with constraint:', methweight
    chain.Draw(drawarg, methweight)
    hrands = chain.GetHistogram().Clone(key+'_rands') 
    hrands.GetYaxis().SetRangeUser(0.01,10000*hrands.GetMaximum())
    histoStyler(hrands, kAzure-8)
    hrands.SetFillColor(hrands.GetLineColor())
    hrands.SetFillStyle(1001)
    hobs.SetTitle('DY+Jets MC')
    hrands.SetTitle('Rebalance and smear DY+Jets')    
    leg = mklegend(x1=.42, y1=.52, x2=.92, y2=.71, color=kWhite)
    hratio, hmethodsyst = FabDrawSystyRatio(c1,leg,hobs,[hrands],datamc='MC',lumi='n/a', title = '', LinearScale=False, fractionthing='truth / method')

    c1.Update()
    c1.Write('c_'+key)
'''

'''
c2 = mkcanvas('c2')
#drawarg = "Jets[JetsRebalanced_origIdx[0]].Pt()/GenJets[0].Pt()>>hadc(100,0,3)"
#constraint = "IsUniqueSeed==1 && fabs(JetsRebalanced[0].Eta()-GenJets[0].Eta())<0.4"

drawarg = "Jets[JetsRebalanced_origIdx[0]].Pt()/GenJets[0].Pt()>>hadc(100,0,3)"
constraint = "IsUniqueSeed==1 && HardMETPt>100 && fabs(JetsRebalanced[0].Eta()-GenJets[0].Eta())<0.4 && fabs(JetsRebalanced[0].Phi()-GenJets[0].Phi())<0.4"

print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histJetResponse = chain.GetHistogram().Clone('histJet1Response')
histoStyler(histJetResponse, kRed+2)
histJetResponse.SetFillStyle(3012)
histJetResponse.SetFillColor(kRed+1)
print 'RMS', histJetResponse.GetRMS()

drawarg = "JetsRebalanced[0].Pt()/GenJets[0].Pt()>>hadc(100,0,3)"
constraint = "IsUniqueSeed==1 && HardMETPt>100 && fabs(JetsRebalanced[0].Eta()-GenJets[0].Eta())<0.4 && fabs(JetsRebalanced[0].Phi()-GenJets[0].Phi())<0.4"
print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histRebalanced = chain.GetHistogram().Clone('hRebalancedJet1Response')
histoStyler(histRebalanced, kBlack)
histRebalanced.SetFillStyle(3005)
histRebalanced.SetFillColor(kBlack)
histRebalanced.SetLineWidth(3)
print 'RMS', histRebalanced.GetRMS()

histRebalanced.GetXaxis().SetTitle('pT(RECO)/pT(GEN)')
histRebalanced.SetTitle('')
histRebalanced.GetYaxis().SetRangeUser(0,1.2*histRebalanced.GetMaximum())
histRebalanced.Draw('hist')
histJetResponse.Draw('hist same')

leg = mklegend(x1=.47, y1=.76, x2=.93, y2=.88, color=kWhite)
leg.AddEntry(histRebalanced, 'rebalanced jet')
leg.AddEntry(histJetResponse, 'original jet')
leg.Draw()
c2.Update()

c2.Write('responseJet1')
histJetResponse.Write()
histRebalanced.Write()

c3 = mkcanvas('c3')
drawarg = "Jets[JetsRebalanced_origIdx[1]].Pt()/GenJets[1].Pt()>>hadc(100,0,3)"
constraint = "IsUniqueSeed==1 && HardMETPt>100 && fabs(JetsRebalanced[1].Eta()-GenJets[1].Eta())<0.4"
print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histJetResponse = chain.GetHistogram().Clone('histJet2Response')
histoStyler(histJetResponse, kRed+2)
histJetResponse.SetFillStyle(3012)
histJetResponse.SetFillColor(kRed+1)
print 'RMS', histJetResponse.GetRMS()

drawarg = "JetsRebalanced[1].Pt()/GenJets[1].Pt()>>hadc(100,0,3"
constraint = "IsUniqueSeed==1 && fabs(JetsRebalanced[1].Eta()-GenJets[1].Eta())<0.4"
print 'now on to this', drawarg
print 'with constraint', constraint
chain.Draw(drawarg, constraint)

histRebalanced = chain.GetHistogram().Clone('hRebalancedJet2Response')
histoStyler(histRebalanced, kBlack)
histRebalanced.SetFillStyle(3005)
histRebalanced.SetFillColor(kBlack)
histRebalanced.SetLineWidth(3)
print 'RMS', histRebalanced.GetRMS()

histRebalanced.GetXaxis().SetTitle('pT(RECO)/pT(GEN)')
histRebalanced.SetTitle('')
histRebalanced.GetYaxis().SetRangeUser(0,1.2*histRebalanced.GetMaximum())
histRebalanced.Draw('hist')
histJetResponse.Draw('hist same')

leg = mklegend(x1=.47, y1=.76, x2=.93, y2=.88, color=kWhite)
leg.AddEntry(histRebalanced, 'rebalanced jet')
leg.AddEntry(histJetResponse, 'original jet')
leg.Draw()
c3.Update()

c3.Write('responseJet2')
histJetResponse.Write()
histRebalanced.Write()
'''

print 'just created', fnew.GetName()
fnew.Close()
