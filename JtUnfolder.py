import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
from ROOT import gROOT
from pickle import BINSTRING
gROOT.ProcessLine( 'gSystem->Load("/Users/tuomas/OneDrive/work/032.JTAnalysis/Unfolding/Root6/RooUnfold/libRooUnfold");')
from rootpy.io import root_open
from ROOT import TMath, TRandom3, TVector3
from ROOT import TF1
from ROOT import TString
from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBayes
import math
from rootpy.plotting import Hist,Hist2D,Hist3D
import matplotlib.pyplot as plt
from datetime import datetime
from datetime import timedelta
from defs import *
import drawing
import sys

class JtUnfolder(object):
  def __init__(self,name,**kwargs):
    self._name = name
    print("Create Unfolder")
    self._Njets = kwargs.get('Njets',0)
    self._fEff = None
    self._jetBinBorders = kwargs.get('jetBinBorders',None)
    if self._jetBinBorders is not None:
      self._jetPtBins = [(a,b) for a,b in zip(self._jetBinBorders,self._jetBinBorders[1:])]
    self._fakeRate = kwargs.get('fakeRate',1)
    self._randomSeed = kwargs.get('randomSeed',123)
    self._weight = kwargs.get('weight',True)
    self._fillFakes = kwargs.get('fillFakes',False)
    self._NBINSJt = kwargs.get('NBINSJt',64)
    self._NBINS=kwargs.get('NBINS',64)
    self._responseJetPt = None
    self._responseJetPtCoarse = None
    self._IsData = kwargs.get('Data',False) 
    self._hJtTrue2D = None
    self._Niter = kwargs.get('Iterations',4)



    
  def setTrackMatch(self,hists):
    self._matching = hists  
    
  def setJtBins(self,bins):
    self._LogBinsJt = bins

  def setPtBins(self,bins):
    self._LogBinsX = bins
        
  def setNbinsJt(self,nbins):
    self._NBINSJt = nbins
    
  def setNbinsPt(self,nbins):
    self._NBINS = nbins  
       
  def setRandomSeed(self,seed):
    self._randomSeed = seed
        
  def setJtTestMeas2D(self,hJtMeas):
    self._hJtTestMeas2D = hJtMeas

  def setJtTestTrue2D(self,hJtTrue):
    self._hJtTestTrue2D = hJtTrue

  def setJtMeas2D(self,hJtMeas):
    self._hJtMeas2D = hJtMeas
  
  def setJtTrue2D(self,hJtTrue):
    self._hJtTrue2D = hJtTrue
  
  def setJetPtMeas(self,hPtMeas):
    self._hJetPtMeas = hPtMeas
  
  def setJetPtTrue(self,hPtTrue):
    self._hJetPtTrue = hPtTrue

  def setJetPtMeasCoarse(self,hPtMeas):
    self._hJetPtMeasCoarse = hPtMeas
  
  def setJetPtTrueCoarse(self,hPtTrue):
    self._hJetPtTrueCoarse = hPtTrue
 
  def setFakeRate(self,rate):
    self._fakeRate = rate
  
  def setWeight(self,weight):
    self._weight = weight
    
  def setMisses2D(self,misses):
    self._misses2D = misses
  
  def setFakes2D(self,fakes):
    self._hJtFake2D = fakes

  def setJetPtResponseHisto(self,hresponse):
    self._hresponseJetPt = hresponse

  def setJetPtResponseCoarseHisto(self,hresponse):
    self._hresponseJetPtCoarse = hresponse

  def setJetPtResponse(self,response):
    self._responseJetPt = response

  def setJetPtResponseCoarse(self,response):
    self._responseJetPtCoarse = response

  def setNumberJetsMeas(self,n):
    self._numberJetsMeasBin = n
  
  def setNumberJetsTestMeas(self,n):
    self._numberJetsTestMeasBin = n

  def setNumberJetsTrue(self,n):
    self._numberJetsTrueBin = n

  def setNumberJetsTestTrue(self,n):
    self._numberJetsTestTrueBin = n
    
  def setNumberJetsMeasTrain(self,n):
    self._numberJetsMeasTrain = n
    
  def set2Dresponse(self,responses):
    self._2Dresponses = responses
  
  def setWeighting(self,wei):
    self._weighting = wei
    
  def drawTrackMatch(self,name,option):
    drawing.drawMatchHisto(self._matching, self._jetPtBins,name,option)  
    
  def createToyTraining(self,rootFile,numberEvents):
    self._f = root_open(rootFile, 'read')
    self._hJetPtIn = self._f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(2))
    self._hZIn = self._f.Get('AliJJetJtTask/AliJJetJtHistManager/Z/ZNFin{:02d}'.format(2))
    LimL=0.1
    LimH=500
    logBW = (TMath.Log(LimH)-TMath.Log(LimL))/self._NBINS
    self._LogBinsX = [LimL*math.exp(ij*logBW) for ij in range(0,self._NBINS+1)]  

    self._hJetPtMeas = Hist(self._LogBinsX)
    self._hJetPtTrue = Hist(self._LogBinsX)  

    self._myRandom = TRandom3(self._randomSeed);
    if self._fEff is None:
      self._fEff = TF1("fEff","1-0.5*exp(-x)")
    if self._jetBinBorders is None:
      self._jetBinBorders = [5,10,20,30,40,60,80,100,150,500]
    
    self._hJetPtMeasCoarse = Hist(self._jetBinBorders)
    self._hJetPtTrueCoarse = Hist(self._jetBinBorders)
     
    low = 0.01
    high = 10
    BinW = (TMath.Log(high)-TMath.Log(low))/self._NBINSJt
    self._LogBinsJt = [low*math.exp(i*BinW) for i in range(self._NBINSJt+1)]
    self._jetBinBorders = self._jetBinBorders
    self._jetPtBins = [(a,b) for a,b in zip(self._jetBinBorders,self._jetBinBorders[1:])]

    self._hJtTrue2D = Hist2D(self._LogBinsJt,self._jetBinBorders)
    self._hJtMeas2D = Hist2D(self._LogBinsJt,self._jetBinBorders)
    self._hJtFake2D = Hist2D(self._LogBinsJt,self._jetBinBorders)
    
    #Histogram to store jT with respect to the leading hadron
    self._hJtTestMeas2D = Hist2D(self._LogBinsJt,self._jetBinBorders) #Needs a better name
    self._hJtTestTrue2D = Hist2D(self._LogBinsJt,self._jetBinBorders) #Needs a better name

    self._responseJetPt = RooUnfoldResponse(self._hJetPtMeas,self._hJetPtTrue)
    self._responseJetPtCoarse = RooUnfoldResponse(self._hJetPtMeasCoarse,self._hJetPtTrueCoarse)
    self._hresponseJetPt = Hist2D(self._jetBinBorders,self._jetBinBorders)
    self._hresponseJetPtCoarse = Hist2D(self._jetBinBorders,self._jetBinBorders)
    #Histogram index is jet pT index, Bin 0 is 5-10 GeV
    #Histogram X axis is observed jT, Bin 0 is underflow
    #Histogram Y axis is observed jet Pt, Bin 0 is underflow
    #Histogram Z axis is True jT, Bin 0 is underflow
    self._2Dresponses = [Hist3D(self._LogBinsJt,self._jetBinBorders,self._LogBinsJt) for i in self._jetPtBins]
    self._misses2D = Hist2D(self._LogBinsJt,self._jetBinBorders)
    self._fakes2D = Hist2D(self._LogBinsJt,self._jetBinBorders)
    

    self._numberJetsMeasBin = [0 for i in self._jetBinBorders]
    self._numberJetsTrueBin = [0 for i in self._jetBinBorders]
    self._numberJetsTestMeasBin = [0 for i in self._jetBinBorders]
    self._numberJetsTestTrueBin = [0 for i in self._jetBinBorders]
    self._numberFakesBin = [0 for i in self._jetBinBorders]
    ieout = numberEvents/10
    if(ieout > 20000):
      ieout = 20000
    start_time = datetime.now()
    print("Processing MC Training Events")
    for ievt in range(numberEvents):
      tracksTrue = []
      tracksMeas = [0 for x in range(100)]
      if ievt%ieout == 0 and ievt > 0:
        time_elapsed = datetime.now() - start_time
        time_left = timedelta(seconds = time_elapsed.total_seconds() * 1.0*(numberEvents-ievt)/ievt)
        print('Event {} [{:.2f}%] Time Elapsed: {} ETA: {}'.format(ievt,100.0*ievt/numberEvents,fmtDelta(time_elapsed),fmtDelta(time_left)))
      jetTrue = TVector3(0,0,0)
      jetMeas = TVector3(0,0,0)
      jetPt = self._hJetPtIn.GetRandom()
      remainder = jetPt
      if jetPt < 5:
          continue
      nt = 0
      nt_meas = 0
      while remainder > 0:
          trackPt = self._hZIn.GetRandom()*jetPt
          if trackPt < remainder:
              track = TVector3()
              remainder = remainder - trackPt
          else:
              trackPt = remainder
              remainder = -1
          if(trackPt > 0.15):
              track.SetPtEtaPhi(trackPt,self._myRandom.Gaus(0,0.1),self._myRandom.Gaus(math.pi,0.2))
              tracksTrue.append(track)
              jetTrue += track
              if self._fEff.Eval(trackPt) > self._myRandom.Uniform(0,1):
                  tracksMeas[nt] = 1
                  jetMeas += track
                  nt_meas += 1
              else:
                  tracksMeas[nt] = 0
              nt +=1
      fakes = []
      for it in range(self._fakeRate*100):
        if(self._myRandom.Uniform(0,1) > 0.99):
          fake = TVector3()
          fake.SetPtEtaPhi(self._myRandom.Uniform(0.15,1),self._myRandom.Gaus(0,0.1),self._myRandom.Gaus(math.pi,0.2))
          fakes.append(fake)
          jetMeas += fake
      
      self._responseJetPt.Fill(jetMeas.Pt(),jetTrue.Pt())
      self._responseJetPtCoarse.Fill(jetMeas.Pt(),jetTrue.Pt())
      self._hresponseJetPt.Fill(jetMeas.Pt(),jetTrue.Pt())
      self._hresponseJetPtCoarse.Fill(jetMeas.Pt(),jetTrue.Pt())
      ij_meas = GetBin(self._jetBinBorders,jetMeas.Pt())
      ij_true = GetBin(self._jetBinBorders,jetTrue.Pt())
      if nt < 5 or nt_meas < 5:
          continue
      if ij_meas >= 0:
        self._numberJetsMeasBin[ij_meas] += 1
      if ij_true >= 0:
        self._numberJetsTrueBin[ij_true] += 1
      for track,it in zip(tracksTrue,range(100)):
        zTrue = (track*jetTrue.Unit())/jetTrue.Mag()
        jtTrue = (track - scaleJet(jetTrue,zTrue)).Mag()
        if ij_meas >= 0:
          if tracksMeas[it] == 1:
            zMeas= (track*jetMeas.Unit())/jetMeas.Mag()
            jtMeas = (track - scaleJet(jetMeas,zMeas)).Mag()
            self._2Dresponses[ij_true].Fill(jtMeas,jetMeas.Pt(),jtTrue)
          else:
            self._misses2D.Fill(jtTrue,jetTrue.Pt())
      if ij_meas >= 0:
        for fake in fakes:
          zFake = (fake*jetMeas.Unit())/jetMeas.Mag()
          jtFake = (fake-scaleJet(jetMeas,zFake)).Mag()
          if self._weight:
            self._hJtFake2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
          else:
            self._hJtFake2D.Fill(jtFake,jetMeas.Pt())
          if self._fillFakes:
            self._fakes2D.Fill(jtFake,jetMeas.Pt())
    print('Event {} [{:.2f}%] Time Elapsed: {}'.format(numberEvents,100.0,fmtDelta(time_elapsed)))      
    self._numberJetsMeasTrain = sum(self._numberJetsMeasBin)
    
  def createToyData(self,rootFile,numberEvents):
    self._f = root_open(rootFile, 'read')
    self._hJetPtIn = self._f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(2))
    self._hZIn = self._f.Get('AliJJetJtTask/AliJJetJtHistManager/Z/ZNFin{:02d}'.format(2))
    if self._fEff is None:
      self._fEff = TF1("fEff","1-0.5*exp(-x)")
    if self._jetBinBorders is None:
      self._jetBinBorders = [5,10,20,30,40,60,80,100,150,500]        
    self._hMultiTrue = Hist(50,0,50)
    self._hMultiMeas = Hist(50,0,50)
    self._hZMeas = Hist(50,0,1)
    self._hZTrue = Hist(50,0,1)
    self._hZFake = Hist(50,0,1)
    self._numberJetsMeasBin = [0 for i in self._jetBinBorders]
    self._numberJetsTrueBin = [0 for i in self._jetBinBorders]    
    print("Create testing data")
    start_time = datetime.now()
    numberEvents = numberEvents
    ieout = numberEvents/10
    if(ieout > 20000):
      ieout = 20000
    for ievt in range(numberEvents):
      tracksTrue = []
      tracksMeas = [0 for x in range(100)]
      if ievt%ieout == 0 and ievt > 0:
        time_elapsed = datetime.now() - start_time
        time_left = timedelta(seconds = time_elapsed.total_seconds() * 1.0*(numberEvents-ievt)/ievt)
        print('Event {} [{:.2f}%] Time Elapsed: {} ETA: {}'.format(ievt,100.0*ievt/numberEvents,fmtDelta(time_elapsed),fmtDelta(time_left)))
      jetTrue = TVector3(0,0,0)
      jetMeas = TVector3(0,0,0)
      jetPt = self._hJetPtIn.GetRandom()
      remainder = jetPt
      if jetPt < 5:
          continue
      nt = 0
      nt_meas = 0
      leading = None
      leadingPt = 0
      while remainder > 0:
          trackPt = self._hZIn.GetRandom()*jetPt
          if trackPt < remainder:
              track = TVector3()
              remainder = remainder - trackPt
          else:
              trackPt = remainder
              remainder = -1
          if(trackPt > 0.15):
              track.SetPtEtaPhi(trackPt,self._myRandom.Gaus(0,0.1),self._myRandom.Gaus(math.pi,0.2))
              if track.Pt() > leadingPt:
                leading = track
                leadingPt = leading.Pt()
              tracksTrue.append(track)
              jetTrue += track
              if self._fEff.Eval(trackPt) > self._myRandom.Uniform(0,1):
                  tracksMeas[nt] = 1
                  jetMeas += track
                  nt_meas += 1
              else:
                  tracksMeas[nt] = 0
              nt +=1
      fakes = []
      if(leadingPt > 0.25*jetPt):
        doLeading = True
      else:
        doLeading = False
      for it in range(self._fakeRate*100):
        if(self._myRandom.Uniform(0,1) > 0.99):
          fake = TVector3()
          fake.SetPtEtaPhi(self._myRandom.Uniform(0.15,1),self._myRandom.Gaus(0,0.1),self._myRandom.Gaus(math.pi,0.2))
          fakes.append(fake)
          jetMeas += fake      
      self._hJetPtMeas.Fill(jetMeas.Pt())
      self._hJetPtTrue.Fill(jetTrue.Pt())
      self._hMultiTrue.Fill(nt)
      self._hMultiMeas.Fill(nt_meas)
      ij_meas = GetBin(self._jetBinBorders,jetMeas.Pt())
      ij_true = GetBin(self._jetBinBorders,jetTrue.Pt())
      if nt < 5 or nt_meas < 5:
          continue
      if ij_meas >= 0:
        self._numberJetsMeasBin[ij_meas] += 1
        self._hJetPtMeasCoarse.Fill(jetMeas.Pt())
        if(doLeading):
          self._numberJetsTestMeasBin[ij_meas] += 1          
      if ij_true >= 0:
        self._numberJetsTrueBin[ij_true] += 1
        self._hJetPtTrueCoarse.Fill(jetTrue.Pt())
        if(doLeading):
          self._numberJetsTestTrueBin[ij_true] += 1
      for track,it in zip(tracksTrue,range(100)):
        zTrue = (track*jetTrue.Unit())/jetTrue.Mag()
        jtTrue = (track - scaleJet(jetTrue,zTrue)).Mag()
        self._hZTrue.Fill(zTrue)
        if(track.Pt() < 0.95*leadingPt and doLeading):
          zLeadingTrue = (track*leading.Unit())/leading.Mag()
          jtLeadingTrue = (track - scaleJet(leading,zLeadingTrue)).Mag()
        if ij_true >= 0:
          if self._weight:
            self._hJtTrue2D.Fill(jtTrue,jetTrue.Pt(),1.0/jtTrue)
            if((track.Pt() < 0.95*leading.Pt()) and doLeading): self._hJtTestTrue2D.Fill(jtLeadingTrue,jetTrue.Pt(),1.0/jtLeadingTrue)
          else:
            self._hJtTrue2D.Fill(jtTrue,jetTrue.Pt())
            if(track.Pt() < 0.95*leading.Pt() and doLeading): self._hJtTestTrue2D.Fill(jtLeadingTrue,jetTrue.Pt())
        if ij_meas >= 0:
          if tracksMeas[it] == 1:
            zMeas= (track*jetMeas.Unit())/jetMeas.Mag()
            jtMeas = (track - scaleJet(jetMeas,zMeas)).Mag()
            if(track.Pt() < 0.95*leading.Pt()):
              zLeadingMeas = (track*leading.Unit())/leading.Mag()
              jtLeadingMeas = (track - scaleJet(leading,zLeadingMeas)).Mag()
            self._hZMeas.Fill(zMeas)
            if self._weight:
              self._hJtMeas2D.Fill(jtMeas,jetMeas.Pt(),1.0/jtMeas)
              if(track.Pt() < 0.95*leadingPt and doLeading): self._hJtTestMeas2D.Fill(jtLeadingMeas,jetMeas.Pt(),1.0/jtLeadingMeas)
            else: 
              self._hJtMeas2D.Fill(jtMeas,jetMeas.Pt())
              if(track.Pt() < 0.95*leadingPt and doLeading): self._hJtTestMeas2D.Fill(jtLeadingMeas,jetMeas.Pt())
      if ij_meas >= 0:
        for fake in fakes:
          zFake = (fake*jetMeas.Unit())/jetMeas.Mag()
          jtFake = (fake-scaleJet(jetMeas,zFake)).Mag()
          zLeadingFake = (fake*leading.Unit())/leading.Mag()
          jtLeadingFake = (fake - scaleJet(leading,zLeadingFake)).Mag()
          self._hZMeas.Fill(zFake)
          self._hZFake.Fill(zFake)
          if self._weight:
            self._hJtMeas2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
            self._hJtTestMeas2D.Fill(jtLeadingFake,leadingPt,1.0/jtLeadingFake)
            #self._hJtFake2D.Fill(jtFake,jetMeas.Pt(),1.0/jtFake)
          else:
            self._hJtMeas2D.Fill(jtFake,jetMeas.Pt())
            self._hJtTestMeas2D.Fill(jtLeadingFake,leadingPt)
            #self._hJtFake2D.Fill(jtFake,jetMeas.Pt())
            
    time_elapsed = datetime.now() - start_time
    print('Event {} [{:.2f}%] Time Elapsed: {}'.format(numberEvents,100.0,fmtDelta(time_elapsed)))      
    
  def unfold(self):
    if self._fillFakes:
      fakes = self._hJtFake2D
    else:
      fakes = None
    if(self._hJtTrue2D):
      self._response2D = make2Dresponse(self._2Dresponses, self._jetPtBins, self._hJtMeas2D, self._hJtTrue2D, self._weighting, misses = self._misses2D, fakes = fakes)
    else: 
      self._response2D = make2Dresponse(self._2Dresponses, self._jetPtBins, self._hJtMeas2D, self._hJtMeas2D, self._weighting, misses = self._misses2D, fakes = fakes)
    if not self._fillFakes:
      print("Scale hJtFake2D by {}".format(1.0*sum(self._numberJetsMeasBin)/self._numberJetsMeasTrain))
      self._hJtFake2D.Scale(1.0*sum(self._numberJetsMeasBin)/self._numberJetsMeasTrain)
      self._hJtMeas2D.Add(self._hJtFake2D,-1)
    
    self._hJtFakeProjBin = [makeHist(self._hJtFake2D.ProjectionX("histFake{}".format(i),i,i), bins = self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
    
    for h,nj in zip(self._hJtFakeProjBin,self._numberJetsMeasBin):
      if(nj > 0):
        h.Scale(1.0/nj,"width")
    if self._responseJetPtCoarse is None:
      print(self._responseJetPtCoarse)
      self._responseJetPtCoarse = createResponse(self._hJetPtMeasCoarse,self._hresponseJetPtCoarse)
    if self._responseJetPt is None:
      self._responseJetPt = createResponse(self._hJetPtMeas, self._hresponseJetPt)
    self._hJetPtRecoCoarse = unfoldJetPt(self._hJetPtMeasCoarse,self._responseJetPtCoarse,self._jetBinBorders)
    self._hJetPtReco = unfoldJetPt(self._hJetPtMeas,self._responseJetPt,self._LogBinsX)

    self._numberJetsMeasFromHist = [self._hJetPtMeasCoarse.GetBinContent(i) for i in range(1,self._hJetPtMeasCoarse.GetNbinsX()+1)]
    if(not self._IsData):
      self._numberJetsTrueFromHist = [self._hJetPtTrueCoarse.GetBinContent(i) for i in range(1,self._hJetPtTrueCoarse.GetNbinsX()+1)]
    self._numberJetsFromReco = [self._hJetPtRecoCoarse.GetBinContent(i) for i in range(1,self._hJetPtRecoCoarse.GetNbinsX())]
    self.printJetNumbers()
    
    unfold2D = RooUnfoldBayes (self._response2D,self._hJtMeas2D,self._Niter)
    self._hJtReco2D = make2DHist(unfold2D.Hreco(),xbins=self._LogBinsJt,ybins=self._jetBinBorders)
    self._hJtMeasBin = [makeHist(self._hJtMeas2D.ProjectionX("histMeas{}".format(i), i,i),bins= self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
    self._hJtTestMeasBin = [makeHist(self._hJtTestMeas2D.ProjectionX("histTestMeas{}".format(i),i,i), bins=self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
    if(not self._IsData):
      self._hJtTrueBin = [makeHist(self._hJtTrue2D.ProjectionX("histMeas{}".format(i), i,i),bins= self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
      self._hJtTestTrueBin = [makeHist(self._hJtTestTrue2D.ProjectionX("histMeas{}".format(i), i,i),bins= self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
    self._hJtRecoBin = [makeHist(self._hJtReco2D.ProjectionX("histMeas{}".format(i), i,i),bins= self._LogBinsJt) for i in range(1,len(self._jetBinBorders))]
    for h,h2,n,n2 in zip(self._hJtMeasBin,self._hJtTestMeasBin,self._numberJetsMeasFromHist,self._numberJetsTestMeasBin):
      if n> 0:
        h.Scale(1.0/n,"width") 
      if n2 > 0:
        h2.Scale(1.0/n2,"width")  
    if(not self._IsData):
      for h,h2,n,n2 in zip(self._hJtTrueBin,self._hJtTestTrueBin,self._numberJetsTrueFromHist,self._numberJetsTestTrueBin):
        if n > 0:
          h.Scale(1.0/n,"width")
        if n2 > 0:
          h2.Scale(1.0/n2,"width")
    for h,n in zip(self._hJtRecoBin,self._numberJetsFromReco):
      if n > 0:
        h.Scale(1.0/n,"width")
        
  def plotJetPt(self):
    if(self._IsData):
      drawing.drawJetPt(self._hJetPtMeas,None,self._hJetPtReco,filename='JetPtUnfoldedData.pdf')
      drawing.drawJetPt(self._hJetPtMeasCoarse,None,self._hJetPtRecoCoarse,filename='JetPtCoarseUnfoldedData.pdf')
    else:
      drawing.drawJetPt(self._hJetPtMeas,self._hJetPtTrue,self._hJetPtReco,filename='JetPtUnfolded.pdf')
      drawing.drawJetPt(self._hJetPtMeasCoarse,self._hJetPtTrueCoarse,self._hJetPtRecoCoarse,filename='JetPtCoarseUnfolded.pdf')
      
  def plotJt(self,filename,**kwargs):
    nRebin = kwargs.get('Rebin',1)
    histlist = [self._hJtMeasBin,self._hJtRecoBin,self._hJtFakeProjBin]
    if(not self._IsData):
      histlist.append(self._hJtTrueBin)
    if nRebin > 0:
      for hists in histlist:
        for h in hists:
          h.Rebin(nRebin)
          h.Scale(1.0/nRebin)
    if(self._IsData):
      drawing.draw8gridcomparison(self._hJtMeasBin,None,self._jetPtBins,xlog = True,ylog = True,name='{}.pdf'.format(filename),unf2d = self._hJtRecoBin, start=4, stride=1)
      drawing.draw8gridcomparison(self._hJtMeasBin,None,self._jetPtBins,xlog = True,ylog = True,name='{}_Extra.pdf'.format(filename),unf2d = self._hJtRecoBin, start=0, stride=1)      
    else:
      drawing.draw8gridcomparison(self._hJtMeasBin,self._hJtTrueBin,self._jetPtBins,xlog = True,ylog = True,name='{}.pdf'.format(filename),unf2d = self._hJtRecoBin,fake= self._hJtFakeProjBin, start=4, stride=1)
      drawing.draw8gridcomparison(self._hJtMeasBin,self._hJtTrueBin,self._jetPtBins,xlog = True,ylog = True,name='{}_Extra.pdf'.format(filename),unf2d = self._hJtRecoBin,fake= self._hJtFakeProjBin, start=0, stride=1)
  
  def plotLeadingJtComparison(self,filename,**kwargs):
    nRebin = kwargs.get('Rebin',1)
    histlist = [self._hJtTestTrueBin,self._hJtTrueBin,self._hJtMeasBin]
    if nRebin > 1:
      for hists in histlist:
        for h in hists:
          h.Rebin(nRebin)
          h.Scale(1.0/nRebin)

    drawing.draw8gridcomparison(self._hJtMeasBin,self._hJtTrueBin,self._jetPtBins,xlog = True,ylog = True,name='{}.pdf'.format(filename), leadingJt=self._hJtTestTrueBin,start=2, stride=1)
    drawing.draw8gridcomparison(self._hJtMeasBin,self._hJtTrueBin,self._jetPtBins,xlog = True,ylog = True,name='{}_Extra.pdf'.format(filename), leadingJt=self._hJtTestTrueBin,start=0, stride=1)
  
  
  def plotResponse(self):
    fig,axs = plt.subplots(2,1,figsize=(10,10))
    axs = axs.reshape(2)
    for ax,h in zip(axs,(self._responseJetPt,self._response2D)):
      hResponse = h.Hresponse()
      xbins = [hResponse.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hResponse.GetNbinsX()+2)]
      ybins = [hResponse.GetYaxis().GetBinLowEdge(iBin) for iBin in range(1,hResponse.GetNbinsY()+2)]
      drawing.drawResponse(ax,make2DHist(hResponse,xbins=xbins,ybins=ybins))
    plt.show() #Draw figure on screen
    
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    hResponse=self._responseJetPt.Hresponse()
    xbins = [hResponse.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hResponse.GetNbinsX()+2)]
    ybins = [hResponse.GetYaxis().GetBinLowEdge(iBin) for iBin in range(1,hResponse.GetNbinsY()+2)]    
    drawing.drawPtResponse(ax,make2DHist(hResponse,xbins=xbins,ybins=ybins))
    plt.savefig("JetPtResponse",format='pdf') #Save figure
    plt.show()

  def printJetNumbers(self):
    print("Measured jets by bin")
    print(self._numberJetsMeasBin)
    print(self._numberJetsMeasFromHist)
    if(not self._IsData):
      print("True jets by bin")
      print(self._numberJetsTrueBin)
      print(self._numberJetsTrueFromHist)
    print("Unfolded jet numbers by bin:")
    print(self._numberJetsFromReco)
    print("Jet bin centers:")
    print([self._hJetPtRecoCoarse.GetBinCenter(i) for i in range(1,self._hJetPtRecoCoarse.GetNbinsX())])

  #TODO
  def writeFiles(self,file):
    print("{}: write results to file".format("JtUnfolder"))
    root_open(file, 'UPDATE')
    #self._hJetPtMeas.Write()
    for h,n in zip(self._hJtRecoBin, range(10)):
      h.SetName("hist{:02d}".format(n))
      h.Write()

def main():
  if(len(sys.argv) > 1):
    numberEvents = int(sys.argv[1])
    if(len(sys.argv) > 2):
      randomSeed = int(sys.argv[2])
    else: randomSeed = 123
  else:
    numberEvents = 5000
  rootFile = "legotrain_350_20161117-2106_LHCb4_fix_CF_pPb_MC_ptHardMerged.root"
  rootFile = 'AnalysisResults.root'
  rootFile = 'legotrain_512_20180523-2331_LHCb4_fix_CF_pPb_MC_ptHardMerged.root'
  unfolderToy = JtUnfolder("ToyMC",NBINSJt=32,NBINS=64,randomSeed=randomSeed)
  unfolderToy.createToyTraining(rootFile, numberEvents)
  unfolderToy.createToyData(rootFile, numberEvents/2)
  unfolderToy.unfold()
  if(numberEvents >= 1000):
    if(numberEvents >= 1000000):
      eString = "{}M_events".format(numberEvents/1000000)
    else:
      eString = "{}k_events".format(numberEvents/1000)
  else:
    eString = "{}_events".format(numberEvents)
  unfolderToy.plotLeadingJtComparison("ToyMCLeadingJtTest_{}".format(eString))
  return
  unfolderToy.plotResponse()
  unfolderToy.plotJetPt()
  unfolderToy.plotJt('ToyMCUnfolder_{}'.format(eString))
  
if __name__ == "__main__": main()
      
