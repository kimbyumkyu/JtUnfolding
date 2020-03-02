from rootpy.plotting import Hist,Hist2D
import rootpy
import math
from ROOT import TMath
from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBayes
from ROOT import TFile
from rootpy.io import root_open
import time

def fmtDelta(delta):
  hours, remainder = divmod(delta.seconds, 3600)
  minutes, seconds = divmod(remainder, 60)
  return '{:02d}:{:02d}:{:02d}'.format(hours,minutes,seconds)

def printTuple(Ntuple):
  for entry in Ntuple:
    jtObs = entry.jtObs
    ptObs = entry.ptObs
    jtTrue = entry.jtTrue
    ptTrue = entry.ptTrue
    print('{},{},{},{}'.format(jtObs,ptObs,jtTrue,ptTrue))
  
def unfoldJetPt(meas,response,xbins):
  unfold = RooUnfoldBayes(response,meas)
  unfold.SetVerbose(0)
  hReco = makeHist(unfold.Hreco(),bins=xbins)
  return hReco

def scaleJets(hist,numberJets):
    h = hist.Clone()
    Njets = h.GetNbinsY()
    Nx = h.GetNbinsX()
    for iy,nj in zip(range(1,Njets+1),numberJets):
        if(nj > 0):
            for ix in range(1,Nx+1):
                iBin = h.GetBin(ix,iy)
                h.SetBinContent(iBin,h.GetBinContent(iBin)/nj)
    return h

def GetBin(array, val):
  for i,a in zip(range(12),array):
    if a > val:
        return i-1
  return -1

def getDiffR(phi1, phi2, eta1, eta2):
  diffPhi = TMath.Abs(phi1-phi2)
  if diffPhi > math.pi: 
    diffPhi = 2*math.pi - diffPhi
  
  return math.sqrt(diffPhi**2+(eta1-eta2)**2)

def makeHist(h,**kwargs):
    if 'bins' in kwargs:
        h2 = Hist(kwargs.get('bins')) 
    else:
        h2 = Hist(h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX()+1),name="_{}".format(h.GetName()))
    for i in range(h.GetNbinsX()):
          h2.SetBinContent(i,h.GetBinContent(i))
    return h2
  
def make2DHist(h,**kwargs):
    if 'xbins' in kwargs and 'ybins' in kwargs:
        h2 = Hist2D(kwargs.get('xbins'),kwargs.get('ybins'))
    else:
        return
    nx = h.GetNbinsX()
    ny = h.GetNbinsY()
    for iy in range(ny):
        for ix in range(nx):
            ib = h.GetBin(ix,iy)
            h2.SetBinContent(ib,h.GetBinContent(ib))
    return h2

def scaleJet(jet,scale):
    jet2 = jet.Clone()
    jet2.SetMag(jet.Mag()*scale)
    return jet2

#Histogram index is jet pT index, Bin 0 is 0-5 GeV
#Histogram X axis is observed jT, Bin 0 is underflow
#Histogram Y axis is observed jet Pt, Bin 0 is underflow
#Histogram Z axis is True jT, Bin 0 is underflow

def make2Dresponse(responses,jetPt,meas,true,weight,misses=None,fakes=None):
  print("++++++++ Create Response from 3D histograms ++++++++")
  response2D = RooUnfoldResponse(meas,true)
  for r,pT,i in zip(responses,jetPt,range(len(responses))):
    nx = r.GetNbinsX()
    ny = r.GetNbinsY()
    nz = r.GetNbinsZ()
    pttrue = (pT[0]+pT[1])/2.0
    print("{:.0f}%".format(100.0*i/len(jetPt)))
    for ix in range(0,nx):
      for iy in range(0,ny):
        for iz in range(1,nz):
          ib = r.GetBin(ix,iy,iz)
          c = r.GetBinContent(ib)
          jttrue = r.GetZaxis().GetBinCenter(iz)
          #if iy == 0 and ix == 0:
            #response2D.Miss(jttrue,pttrue,c)
          if ix > 0 and iy > 0:
            jtobs = r.GetXaxis().GetBinCenter(ix)
            jtobs_meas = meas.GetXaxis().GetBinCenter(ix)
            ptobs_meas = meas.GetYaxis().GetBinCenter(iy)
            if TMath.Abs(jtobs - jtobs_meas) > 0.01:
              print ("jtobs: {}, jtobs_meas: {}".format(jtobs,jtobs_meas))
              raise ValueError("Incorrect binning in make2Dresponse")             
            ptobs = r.GetYaxis().GetBinCenter(iy)
            if TMath.Abs(ptobs - ptobs_meas) > 0.01:
              print ("ptobs: {}, ptobs_meas: {}".format(ptobs,ptobs_meas))
              raise ValueError("Incorrect binning in make2Dresponse") 
            jttrue = r.GetZaxis().GetBinCenter(iz)
            wxbin = weight.GetXaxis().FindBin(jttrue)
            wybin = weight.GetYaxis().FindBin(pttrue)
            wbincontent = weight.GetBinContent(wxbin,wybin)
            if wbincontent > 0.1 and wbincontent<10:
               response2D.Fill(jtobs,ptobs,jttrue,pttrue,c*wbincontent)
            else:
               response2D.Fill(jtobs,ptobs,jttrue,pttrue,c)

            #response2D.Fill(jtobs,ptobs,jttrue,pttrue,c)
  print("{:.0f}%".format(100))

  if misses != None:
    nx = misses.GetNbinsX()
    ny = misses.GetNbinsY()
    for ix in range(1,nx):
      for iy in range(1,ny):
        ib = misses.GetBin(ix,iy)
        c = misses.GetBinContent(ib)
        jttrue = misses.GetXaxis().GetBinCenter(ix)
        pttrue = misses.GetYaxis().GetBinCenter(iy)
        #print("jtTrue: {}, ptTrue: {}, Misses: {}".format(jttrue,pttrue,c))
        response2D.Miss(jttrue,pttrue,c)
  if fakes != None:
    nx = fakes.GetNbinsX()
    ny = fakes.GetNbinsY()
    for ix in range(1,nx):
      for iy in range(1,ny):
        ib = fakes.GetBin(ix,iy)
        c = fakes.GetBinContent(ib)
        jtobs = fakes.GetXaxis().GetBinCenter(ix)
        ptobs = fakes.GetYaxis().GetBinCenter(iy)
        #print("jtObs: {}, ptObs: {}, Fakes: {}".format(jtobs,ptobs,c))
        response2D.Fake(jtobs,ptobs,c)      
  return response2D

def makeResponseFromTuple(Ntuple,meas,true):
  print("Start makeResponseFromTuple")
  start = time.time()
  response2D = RooUnfoldResponse(meas,true)
  for entry in Ntuple:
    jtObs = entry.jtObs
    ptObs = entry.ptObs
    jtTrue = entry.jtTrue
    ptTrue = entry.ptTrue
    if ptObs > 0 and ptTrue > 0:
      response2D.Fill(jtObs,ptObs,jtTrue,ptTrue)
    elif ptObs < 0:
      response2D.Miss(jtTrue,ptTrue)
    elif ptTrue < 0:
      response2D.Fake(jtObs,ptObs)
    else:
      print("ERROR")
  
  end = time.time()
  print("Finished in {}s".format(end - start))
  return response2D

  
def createResponseInverse(hMeas,hResponse):
  response = RooUnfoldResponse(hMeas,hMeas)
  for ibx in range(1,hResponse.GetNbinsX()+1):
    jt = hResponse.GetXaxis().GetBinCenter(ibx)
    ib = hResponse.GetBin(ibx,0)
    N  = hResponse.GetBinContent(ib)
    response.Miss(jt,N)
    for iby in range(1,hResponse.GetNbinsY()+1):
      jtobs = hResponse.GetYaxis().GetBinCenter(iby)
      ib = hResponse.GetBin(ibx,iby)
      N = hResponse.GetBinContent(ib)
      response.Fill(jtobs,jt,N)
  return response

def createResponse(hMeas,hResponse):
  response = RooUnfoldResponse(hMeas,hMeas)
  for iby in range(1,hResponse.GetNbinsY()+1):
    jt = hResponse.GetYaxis().GetBinCenter(iby)
    ib = hResponse.GetBin(iby,0)
    N  = hResponse.GetBinContent(ib)
    response.Miss(jt,N)
    for ibx in range(1,hResponse.GetNbinsX()+1):
      jtobs = hResponse.GetXaxis().GetBinCenter(ibx)
      ib = hResponse.GetBin(ibx,iby)
      N = hResponse.GetBinContent(ib)
      response.Fill(jtobs,jt,N)
  return response
