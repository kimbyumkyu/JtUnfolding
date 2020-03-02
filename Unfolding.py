import logging
mpl_logger = logging.getLogger('matplotlib') 
mpl_logger.setLevel(logging.WARNING) 
from ROOT import gROOT
#gROOT.ProcessLine( 'gSystem->Load("/Users/tuomas/OneDrive/work/032.JTAnalysis/Unfolding/Root6/RooUnfold/libRooUnfold");')
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
from ROOT import TString
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
import time
import JtUnfolder
from ctypes import c_int

def main():
  print('Number of arguments: ', len(sys.argv), 'arguments.')
  print('Argument list:',str(sys.argv))
  filename = sys.argv[1]
  print("Input file: {}".format(filename))
  ff = root_open("weighting.root","read")
  weighting = ff.Get("HoverP")
  
  jetBinBorders = [5,10,20,30,40,60,80,100,150,500]
  jetPtBins = [(a,b) for a,b in zip(jetBinBorders,jetBinBorders[1:])]
  JetPtCenter = [7.5,15,25,35,50,70,90,125,325]
  JetPtError = [2.5,5,5,5,10,10,10,25,175]
  Njets = len(jetBinBorders)-1
  Njets = 9
  f = root_open(filename, 'read')
  if(len(sys.argv) > 2):
    filename_test = sys.argv[2]
    f_test = root_open(filename_test)
    print("Open {} for test data".format(filename_test))
  else:
    f_test = root_open(filename)
  if(len(sys.argv) > 3):
    filename_data = sys.argv[3]
    f_data=root_open(filename_data)
  else:
    f_data = None
  
  nR = 3
  iFinder = 0
  iMCFinder = iFinder + nR*2

  if(f_data):
    numberJetsDataMB = [f_data.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)).Integral() for iJet in range(Njets)]
    hTrackJtDataMB = [f_data.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)] 
    hTrackJt2DDataMB = f_data.Get('AliJJetJtTask/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iFinder))
    hJetPtDataMBCoarse = Hist(jetBinBorders)
    for n,i in zip(numberJetsDataMB,range(1,Njets+1)):
      hJetPtDataMBCoarse.SetBinContent(i,n)
    hJetPtDataMB = f_data.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iFinder))
 
    numberJetsDataTriggered = [f_data.Get('AliJJetJtTask_kEMCEJE/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)).Integral() for iJet in range(Njets)]
    hTrackJtDataTriggered = [f_data.Get('AliJJetJtTask_kEMCEJE/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)] 
    hTrackJt2DDataTriggered = f_data.Get('AliJJetJtTask_kEMCEJE/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iFinder))
    hJetPtDataTriggeredCoarse = Hist(jetBinBorders)
    for n,i in zip(numberJetsDataTriggered,range(1,Njets+1)):
      hJetPtDataTriggeredCoarse.SetBinContent(i,n)
    hJetPtDataTriggered = f_data.Get('AliJJetJtTask_kEMCEJE/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iFinder))    
  
  numberJetsMeas = [f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)).Integral() for iJet in range(Njets)]
  numberJetsTrue = [f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iMCFinder,iJet)).Integral() for iJet in range(Njets)]
  hTrackJtMeas = [f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)] 
  hTrackJtTrue = [f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtWeightBin/JetConeJtWeightBinNFin{:02d}JetPt{:02d}'.format(iMCFinder,iJet)) for iJet in range(Njets)] 
  numberJetsMeasTrain = [f.Get('AliJJetJtTask/AliJJetJtHistManager/JetPtBin/JetPtBinNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)).Integral() for iJet in range(Njets)]
  

  
  print("Measured Jets:")
  print(numberJetsMeas)
  print("True Jets:")
  print(numberJetsTrue)
  hTrackJtCorr2D = [f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackJtCorr2D/TrackJtCorr2DNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)]
  hTrackJtMeas2D = f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iFinder))
  hTrackJtTrue2D = f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JtWeight2D/JtWeight2DNFin{:02d}'.format(iMCFinder))
  hTrackJtMisses2D = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackJtMisses2D/TrackJtMisses2DNFin{:02d}'.format(iFinder))
  hTrackJtFakes2D = f.Get('AliJJetJtTask/AliJJetJtHistManager/JetConeJtUnfBg2D/JetConeJtUnfBg2DNFin{:02d}'.format(iFinder))
  
  print(hTrackJtMeas2D.GetTitle())
  print(hTrackJtTrue2D.GetTitle())
  
  hTrackMatchSuccess = [f.Get('AliJJetJtTask/AliJJetJtMCHistManager/TrackMatchSuccess/TrackMatchSuccessNFin{:02d}JetPt{:02d}'.format(iFinder,iJet)) for iJet in range(Njets)]
  
#   for h in hTrackJtCorr2D:
#     print("Bin 25702 content : {}".format(h.GetBinContent(25702)))
#     ybins = [h.GetYaxis().GetBinLowEdge(iBin) for iBin in range(1,h.GetNbinsY()+2)]
#     print(ybins)
#     ix = c_int()
#     iy = c_int()
#     iz = c_int()
#     h.GetBinXYZ(25702,ix,iy,iz)
#     print("{},{},{}".format(ix,iy,iz))
#     print("{}, Entries: {}".format(h.GetName(),h.GetEntries()))
#     for ibx in range(0,h.GetNbinsX()+1):
#       jtobs = h.GetXaxis().GetBinCenter(ibx)
#       for iby in range(0,h.GetNbinsY()+1):
#         ptobs = h.GetYaxis().GetBinCenter(iby)
#         for ibz in range(0,h.GetNbinsZ()+1):
#           jttrue = h.GetZaxis().GetBinCenter(ibz)
#           ib = h.GetBin(ibx,iby,ibz)
#           content = h.GetBinContent(ib)
#           #print("Bin {},{},{} is number {} (jtobs = {}, ptobs = {}, jttrue = {}, content = {})".format(ibx,iby,ibz,ib,jtobs,ptobs,jttrue,content))
#           if ib == 25702:
#             print("MOI")
#             return
#           if content > 0:
#             print("bin {},{},{} has content {}, (jtobs = {}, ptobs = {}, jttrue = {})".format(ibx,iby,ibz,content,jtobs,ptobs,jttrue))
#   return
  #hTrackJtMeas2D.Add(hTrackJtFakes2D,-1)
  
  LogBinsJt = [hTrackJtMeas2D.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hTrackJtMeas2D.GetNbinsX()+2)]
  print(LogBinsJt)

  for h in hTrackJtCorr2D:
    print("{}".format(h.GetTitle()))
  hJetPtMeas = f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iFinder))
  hJetPtMeasCoarse = Hist(jetBinBorders)
  for n,i in zip(numberJetsMeas,range(1,Njets+1)):
    hJetPtMeasCoarse.SetBinContent(i,n)
    
  hJetPtTrue = f_test.Get('AliJJetJtTask/AliJJetJtHistManager/JetPt/JetPtNFin{:02d}'.format(iMCFinder))
  hJetPtTrueCoarse = Hist(jetBinBorders)
  for n,i in zip(numberJetsTrue,range(1,Njets+1)):
    hJetPtTrueCoarse.SetBinContent(i,n)
  LogBinsPt = [hJetPtTrue.GetXaxis().GetBinLowEdge(iBin) for iBin in range(1,hJetPtTrue.GetNbinsX()+2)]
  hJetPtResponse = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorr/JetPtCorrNFin{:02d}'.format(iFinder))
  hJetPtResponseCoarse = f.Get('AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorrCoarse/JetPtCorrCoarseNFin{:02d}'.format(iFinder))

  
  TrackJtUnfolder = JtUnfolder.JtUnfolder('TrackJtUnfolder',jetBinBorders=jetBinBorders, Njets=Njets,Iterations = 5)
  TrackJtUnfolder.setTrackMatch(hTrackMatchSuccess)
  #TrackJtUnfolder.drawTrackMatch("TrackMatch",'single')
  TrackJtUnfolder.setPtBins(LogBinsPt)
  TrackJtUnfolder.setJtBins(LogBinsJt)
  TrackJtUnfolder.setJtMeas2D(hTrackJtMeas2D)
  TrackJtUnfolder.setJtTestMeas2D(hTrackJtMeas2D)
  TrackJtUnfolder.setJtTrue2D(hTrackJtTrue2D)
  TrackJtUnfolder.setJtTestTrue2D(hTrackJtTrue2D)
  TrackJtUnfolder.setMisses2D(hTrackJtMisses2D)
  TrackJtUnfolder.setFakes2D(hTrackJtFakes2D)
  TrackJtUnfolder.setJetPtMeas(hJetPtMeas)
  TrackJtUnfolder.setJetPtTrue(hJetPtTrue)
  TrackJtUnfolder.setJetPtMeasCoarse(hJetPtMeasCoarse)
  TrackJtUnfolder.setJetPtTrueCoarse(hJetPtTrueCoarse)
  TrackJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas, hJetPtResponse))
  TrackJtUnfolder.setJetPtResponseCoarse(createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse))
  TrackJtUnfolder.setNumberJetsMeas(numberJetsMeas)
  TrackJtUnfolder.setNumberJetsTrue(numberJetsTrue)
  TrackJtUnfolder.setNumberJetsTestMeas(numberJetsMeas)
  TrackJtUnfolder.setNumberJetsTestTrue(numberJetsTrue)
  TrackJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
  TrackJtUnfolder.set2Dresponse(hTrackJtCorr2D)
  TrackJtUnfolder.setWeighting(weighting)
  TrackJtUnfolder.unfold()
  TrackJtUnfolder.plotResponse()
  TrackJtUnfolder.plotJetPt()
  TrackJtUnfolder.plotJt("PythiaTest",Rebin=4)
  TrackJtUnfolder.writeFiles("result.root")
  return
  
#   MBDataJtUnfolder = JtUnfolder.JtUnfolder('MBDataUnfolder',jetBinBorders=jetBinBorders,Njets=Njets,Data=True,Iterations=5)
#   MBDataJtUnfolder.setPtBins(LogBinsPt)
#   MBDataJtUnfolder.setJtBins(LogBinsJt)
#   MBDataJtUnfolder.setJtMeas2D(hTrackJt2DDataMB)
#   MBDataJtUnfolder.setFakes2D(hTrackJtFakes2D)
#   MBDataJtUnfolder.setMisses2D(hTrackJtMisses2D)
#   MBDataJtUnfolder.setJetPtMeas(hJetPtDataMB)
#   MBDataJtUnfolder.setJetPtMeasCoarse(hJetPtDataMBCoarse)
#   MBDataJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas, hJetPtResponse))
#   MBDataJtUnfolder.setJetPtResponseCoarse(createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse))
#   MBDataJtUnfolder.setNumberJetsMeas(numberJetsDataMB)
#   MBDataJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
#   MBDataJtUnfolder.set2Dresponse(hTrackJtCorr2D)
#   MBDataJtUnfolder.unfold()
#   MBDataJtUnfolder.plotJt("MBDataUnfolded",Rebin=4)

  
 # TriggeredDataJtUnfolder = JtUnfolder.JtUnfolder('TriggeredDataUnfolder',jetBinBorders=jetBinBorders,Njets=Njets,Data=True,Iterations=5)
 # TriggeredDataJtUnfolder.setPtBins(LogBinsPt)
 # TriggeredDataJtUnfolder.setJtBins(LogBinsJt)
 # TriggeredDataJtUnfolder.setJtMeas2D(hTrackJt2DDataTriggered)
 # TriggeredDataJtUnfolder.setFakes2D(hTrackJtFakes2D)
 # TriggeredDataJtUnfolder.setMisses2D(hTrackJtMisses2D)
 # TriggeredDataJtUnfolder.setJetPtMeas(hJetPtDataTriggered)
 # TriggeredDataJtUnfolder.setJetPtMeasCoarse(hJetPtDataTriggeredCoarse)
 # TriggeredDataJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas, hJetPtResponse))
 # TriggeredDataJtUnfolder.setJetPtResponseCoarse(createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse))
 # TriggeredDataJtUnfolder.setNumberJetsMeas(numberJetsDataTriggered)
 # TriggeredDataJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
 # TriggeredDataJtUnfolder.set2Dresponse(hTrackJtCorr2D)
 # TriggeredDataJtUnfolder.unfold()
 # TriggeredDataJtUnfolder.plotJetPt()
 # TriggeredDataJtUnfolder.plotJt("TriggeredDataUnfolded",Rebin=4)
  
if __name__ == "__main__": main()
  
