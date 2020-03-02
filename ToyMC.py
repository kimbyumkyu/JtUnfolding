from ROOT import gROOT
from numpy import number
gROOT.ProcessLine( 'gSystem->Load("/Users/tuomas/OneDrive/work/032.JTAnalysis/Unfolding/Root6/RooUnfold/libRooUnfold");')
from ROOT import gRandom, TH1, TH1D, cout, TF1, TH2D, TH2
from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import TMath
from ROOT import TVector3
from ROOT import TRandom3
from ROOT import TNtuple
from ROOT import TFile
import matplotlib
import rootpy.ROOT as ROOT
# from ROOT import RooUnfoldSvd
# from ROOT import RooUnfoldTUnfold
# from ROOT import RooUnfoldIds
import rootpy
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
from rootpy.plotting import Hist,Hist2D,Hist3D
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt
import math
from rootpy.io import root_open
import sys
from drawing import *
from defs import *
from datetime import datetime
from datetime import timedelta 



def main():
  if(len(sys.argv) < 3):
    print("Usage: ToyMC [numberEvents] [randomSeed]")
    return
  numberEvents = int(sys.argv[1])
  seed = int(sys.argv[2])
  print ("==================================== TRAIN ====================================")

  f = root_open("legotrain_350_20161117-2106_LHCb4_fix_CF_pPb_MC_ptHardMerged.root", 'read')
  hJetPt = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(2))
  hZ = f.Get('AliJJetJtTask/AliJJetJtHistManager/Z/ZNFin{:02d}'.format(2))

  FillFakes = False
  dummy_variable = True
  weight = True

  NBINS=50;
  LimL=0.1
  LimH=500
  logBW = (TMath.Log(LimH)-TMath.Log(LimL))/NBINS
  LogBinsX = [LimL*math.exp(ij*logBW) for ij in range(0,NBINS+1)]   
    
    
  hJetPtMeas = Hist(LogBinsX)
  hJetPtTrue = Hist(LogBinsX)

  myRandom = TRandom3(seed);
  fEff = TF1("fEff","1-0.5*exp(-x)")
  jetBinBorders = [5,10,20,30,40,60,80,100,150,500]
  hJetPtMeasCoarse = Hist(jetBinBorders)
  hJetPtTrueCoarse = Hist(jetBinBorders)

  NBINSJt = 64
  low = 0.01
  high = 10
  BinW = (TMath.Log(high)-TMath.Log(low))/NBINSJt
  LogBinsJt = [low*math.exp(i*BinW) for i in range(NBINSJt+1)]
  hJtTrue = Hist(LogBinsJt)
  hJtMeas = Hist(LogBinsJt)
  hJtFake = Hist(LogBinsJt)
  LogBinsPt = jetBinBorders
  jetPtBins = [(a,b) for a,b in zip(jetBinBorders,jetBinBorders[1:])]

  hJtTrue2D = Hist2D(LogBinsJt,LogBinsPt)
  hJtMeas2D = Hist2D(LogBinsJt,LogBinsPt)
  hJtFake2D = Hist2D(LogBinsJt,LogBinsPt)
  hJtMeasBin = [Hist(LogBinsJt) for i in jetBinBorders]
  hJtTrueBin = [Hist(LogBinsJt) for i in jetBinBorders]
  
  response= RooUnfoldResponse (hJtMeas,hJtTrue)
  response2D= RooUnfoldResponse (hJtMeas2D,hJtTrue2D)
  responseBin = [RooUnfoldResponse(hJtMeas,hJtTrue) for i in jetBinBorders]
  responseJetPt = RooUnfoldResponse(hJetPtMeas,hJetPtTrue)
  responseJetPtCoarse = RooUnfoldResponse(hJetPtMeasCoarse,hJetPtTrueCoarse)
  
  #Histogram index is jet pT index, Bin 0 is 5-10 GeV
  #Histogram X axis is observed jT, Bin 0 is underflow
  #Histogram Y axis is observed jet Pt, Bin 0 is underflow
  #Histogram Z axis is True jT, Bin 0 is underflow
  responses = [Hist3D(LogBinsJt,LogBinsPt,LogBinsJt) for i in jetPtBins]
  misses = Hist2D(LogBinsJt,LogBinsPt)
  fakes2D = Hist2D(LogBinsJt,LogBinsPt)
  outFile = TFile("tuple.root","recreate")
  responseTuple = TNtuple("responseTuple","responseTuple",'jtObs:ptObs:jtTrue:ptTrue')
  
  hMultiTrue = Hist(50,0,50)
  hMultiMeas = Hist(50,0,50)
  hZMeas = Hist(50,0,1)
  hZTrue = Hist(50,0,1)
  hZFake = Hist(50,0,1)
  responseMatrix = Hist2D(LogBinsJt,LogBinsJt)
  numberJets = 0
  numberFakes = 0
  numberJetsMeasBin = [0 for i in jetBinBorders]
  numberJetsTrueBin = [0 for i in jetBinBorders]
  numberFakesBin = [0 for i in jetBinBorders]
  ieout = numberEvents/10
  if(ieout > 10000):
    ieout = 10000
  fakeRate = 1
  start_time = datetime.now()
  print("Processing Training Events")
  for ievt in range(numberEvents):
    tracksTrue = []
    tracksMeas = [0 for x in range(100)]
    if ievt%ieout == 0 and ievt > 0:
      time_elapsed = datetime.now() - start_time
      time_left = timedelta(seconds = time_elapsed.total_seconds() * 1.0*(numberEvents-ievt)/ievt)
      print('Event {} [{:.2f}%] Time Elapsed: {} ETA: {}'.format(ievt,100.0*ievt/numberEvents,fmtDelta(time_elapsed),fmtDelta(time_left)))
    jetTrue = TVector3(0,0,0)
    jetMeas = TVector3(0,0,0)
    jetPt = hJetPt.GetRandom()
    remainder = jetPt
    if jetPt < 5:
        continue
    nt = 0
    nt_meas = 0
    while remainder > 0:
        trackPt = hZ.GetRandom()*jetPt
        if trackPt < remainder:
            track = TVector3()
            remainder = remainder - trackPt
        else:
            trackPt = remainder
            remainder = -1
        if(trackPt > 0.15):
            track.SetPtEtaPhi(trackPt,myRandom.Gaus(0,0.1),myRandom.Gaus(math.pi,0.2))
            tracksTrue.append(track)
            jetTrue += track
            if fEff.Eval(trackPt) > myRandom.Uniform(0,1):
                tracksMeas[nt] = 1
                jetMeas += track
                nt_meas += 1
            else:
                tracksMeas[nt] = 0
            nt +=1
    fakes = []
    for it in range(fakeRate*100):
      if(myRandom.Uniform(0,1) > 0.99):
        fake = TVector3()
        fake.SetPtEtaPhi(myRandom.Uniform(0.15,1),myRandom.Gaus(0,0.1),myRandom.Gaus(math.pi,0.2))
        fakes.append(fake)
        jetMeas += fake
    
    hJetPtMeas.Fill(jetMeas.Pt())
    hJetPtTrue.Fill(jetTrue.Pt())
    responseJetPt.Fill(jetMeas.Pt(),jetTrue.Pt())
    responseJetPtCoarse.Fill(jetMeas.Pt(),jetTrue.Pt())
    hMultiTrue.Fill(nt)
    hMultiMeas.Fill(nt_meas)
    ij_meas = GetBin(jetBinBorders,jetMeas.Pt())
    ij_true = GetBin(jetBinBorders,jetTrue.Pt())
    if nt < 5 or nt_meas < 5:
        continue
    numberJets += 1
    if ij_meas >= 0:
      numberJetsMeasBin[ij_meas] += 1
      hJetPtMeasCoarse.Fill(jetMeas.Pt())
    if ij_true >= 0:
      numberJetsTrueBin[ij_true] += 1
      hJetPtTrueCoarse.Fill(jetTrue.Pt())
    for track,it in zip(tracksTrue,range(100)):
      zTrue = (track*jetTrue.Unit())/jetTrue.Mag()
      jtTrue = (track - scaleJet(jetTrue,zTrue)).Mag()
      hZTrue.Fill(zTrue)
      if ij_true >= 0:
        if weight:
          hJtTrue.Fill(jtTrue,1.0/jtTrue)
          hJtTrueBin[ij_true].Fill(jtTrue,1.0/jtTrue)
          hJtTrue2D.Fill(jtTrue,jetTrue.Pt(),1.0/jtTrue)
        else:
          hJtTrue.Fill(jtTrue)
          hJtTrueBin[ij_true].Fill(jtTrue)
          hJtTrue2D.Fill(jtTrue,jetTrue.Pt())
      if ij_meas >= 0:
        if tracksMeas[it] == 1:
          zMeas= (track*jetMeas.Unit())/jetMeas.Mag()
          jtMeas = (track - scaleJet(jetMeas,zMeas)).Mag()
          hZMeas.Fill(zMeas)
          if weight:
            hJtMeasBin[ij_meas].Fill(jtMeas,1.0/jtMeas)
            hJtMeas.Fill(jtMeas,1.0/jtMeas)
            hJtMeas2D.Fill(jtMeas,jetMeas.Pt(),1.0/jtMeas)
          else: 
            hJtMeas.Fill(jtMeas)
            hJtMeasBin[ij_meas].Fill(jtMeas)
            hJtMeas2D.Fill(jtMeas,jetMeas.Pt())
          response.Fill(jtMeas,jtTrue)
          responseBin[ij_true].Fill(jtMeas,jtTrue)
          response2D.Fill(jtMeas,jetMeas.Pt(),jtTrue,jetTrue.Pt())
          responseMatrix.Fill(jtMeas,jtTrue)
          responses[ij_true].Fill(jtMeas,jetMeas.Pt(),jtTrue)
          responseTuple.Fill(jtMeas,jetMeas.Pt(),jtTrue,jetTrue.Pt())
        else:
          response.Miss(jtTrue)
          responseBin[ij_true].Miss(jtTrue)
          response2D.Miss(jtTrue,jetTrue.Pt())
          misses.Fill(jtTrue,jetTrue.Pt())
          responseTuple.Fill(-1,-1,jtTrue,jetTrue.Pt())
    if ij_meas >= 0:
      for fake in fakes:
        zFake = (fake*jetMeas.Unit())/jetMeas.Mag()
        jtFake = (fake-scaleJet(jetMeas,zFake)).Mag()
        hZMeas.Fill(zFake)
        hZFake.Fill(zFake)
        if weight:
          hJtMeas.Fill(jtFake,1.0/jtFake)
          hJtMeasBin[ij_meas].Fill(jtFake,1.0/jtFake)
          hJtMeas2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
          hJtFake2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
          hJtFake.Fill(jtFake,1.0/jtFake)
        else:
          hJtMeas.Fill(jtFake)
          hJtMeasBin[ij_meas].Fill(jtFake)
          hJtMeas2D.Fill(jtFake,jetMeas.Pt())
          hJtFake2D.Fill(jtFake,jetMeas.Pt())
          hJtFake.Fill(jtFake)
        if FillFakes:
          response.Fake(jtFake)
          responseBin[ij_true].Fake(jtFake)
          response2D.Fake(jtFake,jetMeas.Pt())
          fakes2D.Fill(jtFake,jetMeas.Pt())
          responseTuple.Fill(jtFake,jetMeas.Pt(),-1,-1)
          numberFakes += 1
          numberFakesBin[ij_true] += 1

  response2Dtest = make2Dresponse(responses, jetPtBins, hJtMeas2D,hJtTrue2D,misses = misses, fakes = fakes2D)

  if dummy_variable:
    hJetPtMeas.Reset()
    hJetPtTrue.Reset()
    hMultiTrue.Reset()
    hMultiMeas.Reset()
    hJetPtMeasCoarse.Reset()
    hJetPtTrueCoarse.Reset()
    hZTrue.Reset()
    hZMeas.Reset()
    hJtTrue.Reset()
    hJtTrue2D.Reset()
    hJtMeas.Reset()
    hJtMeas2D.Reset()
    hJtFake.Reset()
    hJtFake2D.Reset()
    for h,h2 in zip(hJtTrueBin,hJtMeasBin):
      h.Reset()
      h2.Reset()
    numberJetsMeasBin = [0 for i in jetBinBorders]
    numberJetsTrueBin = [0 for i in jetBinBorders]    
    numberJets = 0
    print("Create testing data")
    start_time = datetime.now()
    numberEvents = numberEvents/2
    for ievt in range(numberEvents):
      tracksTrue = []
      tracksMeas = [0 for x in range(100)]
      if ievt%ieout == 0 and ievt > 0:
        time_elapsed = datetime.now() - start_time
        time_left = timedelta(seconds = time_elapsed.total_seconds() * 1.0*(numberEvents-ievt)/ievt)
        print('Event {} [{:.2f}%] Time Elapsed: {} ETA: {}'.format(ievt,100.0*ievt/numberEvents,fmtDelta(time_elapsed),fmtDelta(time_left)))
      jetTrue = TVector3(0,0,0)
      jetMeas = TVector3(0,0,0)
      jetPt = hJetPt.GetRandom()
      remainder = jetPt
      if jetPt < 5:
          continue
      nt = 0
      nt_meas = 0
      while remainder > 0:
          trackPt = hZ.GetRandom()*jetPt
          if trackPt < remainder:
              track = TVector3()
              remainder = remainder - trackPt
          else:
              trackPt = remainder
              remainder = -1
          if(trackPt > 0.15):
              track.SetPtEtaPhi(trackPt,myRandom.Gaus(0,0.1),myRandom.Gaus(math.pi,0.2))
              tracksTrue.append(track)
              jetTrue += track
              if fEff.Eval(trackPt) > myRandom.Uniform(0,1):
                  tracksMeas[nt] = 1
                  jetMeas += track
                  nt_meas += 1
              else:
                  tracksMeas[nt] = 0
              nt +=1
      fakes = []
      for it in range(fakeRate*100):
        if(myRandom.Uniform(0,1) > 0.99):
          fake = TVector3()
          fake.SetPtEtaPhi(myRandom.Uniform(0.15,1),myRandom.Gaus(0,0.1),myRandom.Gaus(math.pi,0.2))
          fakes.append(fake)
          jetMeas += fake      
      hJetPtMeas.Fill(jetMeas.Pt())
      hJetPtTrue.Fill(jetTrue.Pt())
      hMultiTrue.Fill(nt)
      hMultiMeas.Fill(nt_meas)
      ij_meas = GetBin(jetBinBorders,jetMeas.Pt())
      ij_true = GetBin(jetBinBorders,jetTrue.Pt())
      if nt < 5 or nt_meas < 5:
          continue
      numberJets += 1
      if ij_meas >= 0:
        numberJetsMeasBin[ij_meas] += 1
        hJetPtMeasCoarse.Fill(jetMeas.Pt())
      if ij_true >= 0:
        numberJetsTrueBin[ij_true] += 1
        hJetPtTrueCoarse.Fill(jetTrue.Pt())
      for track,it in zip(tracksTrue,range(100)):
        zTrue = (track*jetTrue.Unit())/jetTrue.Mag()
        jtTrue = (track - scaleJet(jetTrue,zTrue)).Mag()
        hZTrue.Fill(zTrue)
        if ij_true >= 0:
          if weight:
            hJtTrue.Fill(jtTrue,1.0/jtTrue)
            hJtTrueBin[ij_true].Fill(jtTrue,1.0/jtTrue)
            hJtTrue2D.Fill(jtTrue,jetTrue.Pt(),1.0/jtTrue)
          else:
            hJtTrue.Fill(jtTrue)
            hJtTrueBin[ij_true].Fill(jtTrue)
            hJtTrue2D.Fill(jtTrue,jetTrue.Pt())
        if ij_meas >= 0:
          if tracksMeas[it] == 1:
            zMeas= (track*jetMeas.Unit())/jetMeas.Mag()
            jtMeas = (track - scaleJet(jetMeas,zMeas)).Mag()
            hZMeas.Fill(zMeas)
            if weight:
              hJtMeasBin[ij_meas].Fill(jtMeas,1.0/jtMeas)
              hJtMeas.Fill(jtMeas,1.0/jtMeas)
              hJtMeas2D.Fill(jtMeas,jetMeas.Pt(),1.0/jtMeas)
            else: 
              hJtMeas.Fill(jtMeas)
              hJtMeasBin[ij_meas].Fill(jtMeas)
              hJtMeas2D.Fill(jtMeas,jetMeas.Pt())
      if ij_meas >= 0:
        for fake in fakes:
          zFake = (fake*jetMeas.Unit())/jetMeas.Mag()
          jtFake = (fake-scaleJet(jetMeas,zFake)).Mag()
          hZMeas.Fill(zFake)
          hZFake.Fill(zFake)
          if weight:
            hJtMeas.Fill(jtFake,1.0/jtFake)
            hJtMeasBin[ij_meas].Fill(jtFake,1.0/jtFake)
            hJtMeas2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
            hJtFake2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
            hJtFake.Fill(jtFake,1.0/jtFake)
          else:
            hJtMeas.Fill(jtFake)
            hJtMeasBin[ij_meas].Fill(jtFake)
            hJtMeas2D.Fill(jtFake,jetMeas.Pt())
            hJtFake2D.Fill(jtFake,jetMeas.Pt())
            hJtFake.Fill(jtFake)
            
  time_elapsed = datetime.now() - start_time
  print('Event {} [{:.2f}%] Time Elapsed: {}'.format(numberEvents,100.0,fmtDelta(time_elapsed)))          
  if not FillFakes:
    hJtMeas.Add(hJtFake,-1)
    hJtMeas2D.Add(hJtFake2D,-1)  
  responseTuple.Print()
  outFile.Write()
#  printTuple(responseTuple)
  
  hJtMeasProjBin = [makeHist(hJtMeas2D.ProjectionX("histMeas{}".format(i), i,i),bins= LogBinsJt) for i in range(1,len(jetBinBorders))]
  hJtMeasProj = makeHist(hJtMeas2D.ProjectionX("histMeas"),bins=LogBinsJt)
  hJtTrueProjBin = [makeHist(hJtTrue2D.ProjectionX("histTrue{}".format(i), i,i),bins = LogBinsJt) for i in range(1,len(jetBinBorders))]
  hJtTrueProj = makeHist(hJtTrue2D.ProjectionX("histTrue"),bins=LogBinsJt)
  hJtFakeProjBin = [makeHist(hJtFake2D.ProjectionX("histFake{}".format(i),i,i), bins = LogBinsJt) for i in range(1,len(jetBinBorders))]
  
  if not FillFakes:
    for h,h2 in zip(hJtMeasBin,hJtFakeProjBin):
      h.Add(h2,-1)
  
  for h in (hJtMeasProj,hJtTrueProj,hJtMeas,hJtTrue,hJtFake,hZFake,hZMeas,hZTrue):
    h.Scale(1.0/numberJets,"width")
  for meas,true,n_meas,n_true in zip(hJtMeasBin,hJtTrueBin,numberJetsMeasBin,numberJetsTrueBin):
      if n_meas > 0:
          meas.Scale(1.0/n_meas,"width")
      if n_true > 0:
          true.Scale(1.0/n_true,"width")
  
  numberJetsMeasFromHist = [hJetPtMeasCoarse.GetBinContent(i) for i in range(1,hJetPtMeasCoarse.GetNbinsX()+1)]
  numberJetsTrueFromHist = [hJetPtTrueCoarse.GetBinContent(i) for i in range(1,hJetPtTrueCoarse.GetNbinsX()+1)]
  print("Total number of jets: {}".format(numberJets))
  print("Total number of fakes: {}".format(numberFakes))
  print("Measured jets by bin")
  print(numberJetsMeasBin)
  print(numberJetsMeasFromHist)
  print("True jets by bin")
  print(numberJetsTrueBin)
  print(numberJetsTrueFromHist)
  hRecoJetPtCoarse = unfoldJetPt(hJetPtMeasCoarse,responseJetPtCoarse,jetBinBorders)
  numberJetsFromReco = [hRecoJetPtCoarse.GetBinContent(i) for i in range(1,hRecoJetPtCoarse.GetNbinsX())]
  print("Unfolded jet numbers by bin:")
  print(numberJetsFromReco)
  
  print("Fakes by bin")
  print(numberFakesBin)

  print ("==================================== UNFOLD ===================================")
  unfold= RooUnfoldBayes     (response, hJtMeas, 4)    #  OR
  unfoldSVD= RooUnfoldSvd     (response, hJtMeas, 20)     #  OR
  unfold2D = RooUnfoldBayes (response2D, hJtMeas2D,4)
  for u in (unfold,unfoldSVD,unfold2D):
    u.SetVerbose(0)
  #response2Dtest = makeResponseFromTuple(responseTuple,hJtMeas2D,hJtTrue2D)


  unfold2Dtest = RooUnfoldBayes (response2Dtest,hJtMeas2D,4)
  
  unfoldBin = [RooUnfoldBayes (responseBin[i],hJtMeasBin[i]) for i in range(len(jetBinBorders))]
  for u in unfoldBin:
    u.SetVerbose(0)
  hRecoBayes= makeHist(unfold.Hreco(),bins=LogBinsJt)
  hRecoSVD= makeHist(unfoldSVD.Hreco(),bins=LogBinsJt)
  hRecoBin = [makeHist(unfoldBin[i].Hreco(),bins=LogBinsJt) for i in range(len(jetBinBorders))]
  hReco2D = make2DHist(unfold2D.Hreco(),xbins=LogBinsJt,ybins=LogBinsPt)
  hReco2Dtest = make2DHist(unfold2Dtest.Hreco(),xbins=LogBinsJt,ybins=LogBinsPt)
  hRecoJetPt = unfoldJetPt(hJetPtMeas,responseJetPt,LogBinsX)



  hReco2DProjBin = [makeHist(hReco2D.ProjectionX("histReco{}".format(i),i,i), bins = LogBinsJt) for i in range(1,len(jetBinBorders))]
  hReco2DTestProjBin = [makeHist(hReco2Dtest.ProjectionX("histRecoTest{}".format(i),i,i), bins=LogBinsJt) for i in range(1,len(jetBinBorders))]
  
  hReco2DProj = makeHist(hReco2D.ProjectionX("histReco"),bins = LogBinsJt)
  hReco2DProj.Scale(1.0/numberJets,"width")
  for h,h2, n in zip(hReco2DProjBin,hReco2DTestProjBin,numberJetsFromReco):
      if n > 0:
          h.Scale(1.0/n,"width")
          h2.Scale(1.0/n,"width")
  #unfold.PrintTable (cout, hJtTrue)
  for h,h2,nj in zip(hJtMeasProjBin,hJtFakeProjBin,numberJetsMeasBin):
      if(nj > 0):
          h.Scale(1.0/nj,"width")
          h2.Scale(1.0/nj,"width")
      #else:
      #    print("nj is 0 for {}".format(h.GetName()))
  for h,nj in zip(hJtTrueProjBin,numberJetsTrueBin):
      if(nj > 0):
          h.Scale(1.0/nj,"width")
          
  #draw8grid(hJtMeasBin[1:],hJtTrueBin[1:],jetPtBins[1:],xlog = True,ylog = True,name="newfile.pdf",proj = hJtMeasProjBin[2:], unf2d = hReco2DProjBin[2:], unf=hRecoBin[1:])
  if(numberEvents > 1000):
    if(numberEvents > 1000000):
      filename = 'ToyMC_{}M_events.pdf'.format(numberEvents/1000000)
    else:
      filename = 'ToyMC_{}k_events.pdf'.format(numberEvents/1000)
  else:
    filename='ToyMC_{}_events.pdf'.format(numberEvents)
  draw8gridcomparison(hJtMeasBin,hJtTrueBin,jetPtBins,xlog = True,ylog = True,name=filename,proj = None, unf2d = hReco2DProjBin, unf2dtest=hReco2DTestProjBin, unf=hRecoBin,fake=hJtFakeProjBin, start=1, stride=1)
  drawQA(hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2DProj,hZ,hZTrue,hZMeas,hZFake,hMultiMeas,hMultiTrue,hJetPt,hJetPtTrue,hJetPtMeas,hRecoJetPt,responseMatrix)
  outFile.Close()




if __name__ == "__main__": main()

