import matplotlib
import rootpy.ROOT as ROOT
import rootpy
from matplotlib.colors import LogNorm
from matplotlib.colors import PowerNorm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LogLocator
from matplotlib.backends.backend_pdf import PdfPages
from rootpy.plotting import Hist,Hist2D
import rootpy.plotting.root2matplotlib as rplt
import matplotlib.pyplot as plt

def drawJetPt(hMeas,hTrue,hReco,filename='JetPt'): 
  if(hTrue is None):
    divider = hMeas
  else:
    divider = hTrue
    hTrue.SetMarkerColor('red')

  hMeas.SetMarkerColor('blue')
  hReco.SetMarkerColor(3)
  ratioMeas = hMeas.Clone()
  ratioMeas.Divide(divider)
  ratioReco = hReco.Clone()
  ratioReco.Divide(divider)
  
  fig, axs = plt.subplots(2,1,figsize=(5,8),sharey=False,sharex=True)
  axs = axs.reshape(2)
  axs[0].text(45,1e-2,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Full jets\n' r'Anti-$k_T$, R=0.4',fontsize=7)
  axs[0].set_ylabel(r'$\frac{dN}{dp_{T}}$',fontsize=18)
  axs[0].set_xlabel(r'$p_{T}$',fontsize = 18)
  axs[1].set_xlabel(r'$p_{T}$',fontsize = 18)
  axs[1].set_ylabel('Ratio',fontsize=18)

  ax = axs[0]
  ax.set_xscale('log') #Set logarithmic scale
  ax.set_yscale('log')
  ax.set_xlim([5,500]) #Set x-axis limits
  ax.set_ylim([1,500]) #Set y-axis limits
  
  rplt.errorbar(hMeas,xerr=False,emptybins=False,axes=ax,label='Measured',fmt='go') #Plot measured pT
  rplt.errorbar(hReco,xerr=False,emptybins=False,axes=ax,label='Unfolded',fmt='go') #Plot unfolded pT
  if(hTrue is not None):
    rplt.errorbar(hTrue,xerr=False,emptybins=False,axes=ax,label='True',fmt='go') #Plot True pT
  
  maxContent = hMeas.GetBinContent(hMeas.GetMaximumBin())
  ax.set_ylim([3e-10*maxContent,10*maxContent])
  ax.legend(loc ='lower left')
  ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
  ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
  #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 

  ax.grid(True) #Draw grid 

  ax= axs[1]
  rplt.errorbar(ratioMeas,xerr=False,emptybins=False,axes=ax,fmt='o')
  rplt.errorbar(ratioReco,xerr=False,emptybins=False,axes=ax,fmt='o')
  ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
  ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
  ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
  #ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
  ax.set_xscale('log') #Set logarithmic scale
  #ax.set_yscale('log')
  ax.set_xlim([5,500]) #Set x-axis limits
  ax.set_ylim([0,2]) #Set y-axis limits
  ax.grid(True) #Draw grid  

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  plt.savefig(filename,format='pdf') #Save figure
  plt.show() #Draw figure on screen

def draw2D(hist,title):
  fig, axs = plt.subplots(1,1,figsize = (10,10), sharex = True, sharey = True)
  #axs = axs.reshape(4)
  ax = axs
  ax.set_xlabel(r'$j_T$')
  ax.set_ylabel(r'$p_{T,jet}$')
  ax.set_xscale('log')
  #ax.set_yscale('log')
  #rplt.hist2d(hist,label="Hist",norm=PowerNorm(0.8),colorbar=True)
  rplt.hist2d(hist, label = "Hist",norm=LogNorm(),colorbar=True)
  ax.set_xlim([0.01,10])
  ax.set_ylim([5,100])
  ax.text(1,7,title,fontsize = 10,bbox=dict(boxstyle="round",
                 ec=(1., 0.5, 0.5),
                 fc=(1., 0.8, 0.8),
                 alpha = 0.8
                 ))
  plt.show()

def draw2DComparison(hists,titles):
  fig, axs = plt.subplots(1,2,figsize = (10,10))
  axs = axs.reshape(2)
  for ax,h,title in zip(axs,hists,titles):
      ax.set_xlabel(r'$j_T$')
      ax.set_ylabel(r'$p_{T,jet}$')
      ax.set_xscale('log')
      #ax.set_yscale('log')
      #rplt.hist2d(hist,label="Hist",norm=PowerNorm(0.8),colorbar=True)
      rplt.hist2d(h,axes = ax, label = title,norm=LogNorm(),colorbar=True)
      ax.set_xlim([0.01,10])
      ax.set_ylim([5,100])
      ax.text(1,7,title,fontsize = 10,bbox=dict(boxstyle="round",
                     ec=(1., 0.5, 0.5),
                     fc=(1., 0.8, 0.8),
                     alpha = 0.8
                     ))
  plt.show()

def drawQA(hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2D,hZ,hZTrue,hZMeas,hZFake,hMultiMeas,hMultiTrue,hJetPt,hJetPtTrue,hJetPtMeas,hJetPtReco,responseM):
  fig, axs = plt.subplots(2,3,figsize = (15,10))
  axs = axs.reshape(6)
  axs[2].yaxis.set_ticks_position('right')
  axs[2].yaxis.set_ticks_position('right')
  axs[5].yaxis.set_label_position('right')
  axs[5].yaxis.set_label_position('right')
  for ax in axs:
    ax.grid(True)
  drawJt(axs[0],hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2D)
  hZInput = hZ.Clone()
  hZInput.Scale(hZTrue.Integral()/hZ.Integral())
  drawJtRatio(axs[3],hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2D)
  drawZ(axs[1],hZMeas,hZTrue,hZInput,hZFake)
  drawMulti(axs[2],hMultiMeas,hMultiTrue)
  hJetPtInput = hJetPt.Clone()
  hJetPtInput.Scale(hJetPtTrue.Integral()/hJetPt.Integral())
  drawJetPtQA(axs[4],hJetPtMeas,hJetPtTrue,hJetPtInput,hJetPtReco)
  drawResponse(axs[5],responseM)
  plt.show() #Draw figure on screen

def drawResponse(ax,hist):
  ax.set_xlabel(r'$j_{T,obs}\left[GeV\right]$') #Add x-axis labels for bottom row
  ax.set_ylabel(r'$j_{T,true}\left[GeV\right]$') #Add x-axis labels for bottom row
  #rplt.hist2d(hist, label = "Hist",norm=LogNorm(),colorbar=True)
  rplt.hist2d(hist, label = "Hist")
  ax.set_xlim([0.01,10])
  ax.set_ylim([0.01,10])
  #ax.set_xscale('log')
  #ax.set_yscale('log')

def drawPtResponse(ax,hist,filename=""):
  ax.set_xlabel(r'$p_{T,obs}\left[GeV\right]$') #Add x-axis labels for bottom row
  ax.set_ylabel(r'$p_{T,true}\left[GeV\right]$') #Add x-axis labels for bottom row
  #rplt.hist2d(hist, label = "Hist",norm=LogNorm(),colorbar=True)
  rplt.hist2d(hist, label = "Hist",norm=LogNorm())
  ax.set_xlim([5,150])
  ax.set_ylim([5,150])
  #ax.set_xscale('log')
  #ax.set_yscale('log')
    

def drawJtRatio(ax,hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2D):

  ax.set_xlabel(r'$j_T$')
  ax.set_ylabel('Ratio to Truth',fontsize=12)
  ax.set_xscale('log')
  
  hJtMeas.linecolor = 'blue'
  hJtTrue.linecolor = 'red'
  hRecoBayes.linecolor = 'green'
  hRecoSVD.linecolor = 'black'
  hReco2D.linecolor = 'cyan'
  hJtFake.linecolor = 'yellow'
  
  h1 = hJtMeas.Clone()
  h2 = hJtTrue.Clone()
  h3 = hRecoBayes.Clone()
  h4 = hRecoSVD.Clone()
  h5 = hReco2D.Clone()
  h6 = hJtFake.Clone()
  for h in (h1,h2,h3,h4,h5,h6):
    h.Divide(hJtTrue)
    
  rplt.hist(h1, axes = ax, label = 'Measured')
  rplt.hist(h2, axes = ax,label='True') #Plot jT histogram,
  rplt.hist(h3, axes = ax, label = 'Bayes')
  rplt.hist(h4, axes = ax, label = 'SVD')
  rplt.hist(h5, axes = ax, label = '2D Unfolded')
  rplt.hist(h6, axes = ax, label = 'Fake')
  ax.set_xlim([0.01,10])
  ax.set_ylim([0,1.5])


def drawJt(ax,hJtMeas,hJtTrue,hJtFake,hRecoBayes,hRecoSVD,hReco2D):
  ax.set_xlim([0.01,10])
  ax.set_ylim([1e-5,1000])
  ax.set_xlabel(r'$j_T$')
  ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18)
  ax.set_yscale('log')
  ax.set_xscale('log')

  hJtMeas.linecolor = 'blue'
  hJtTrue.linecolor = 'red'
  hRecoBayes.linecolor = 'green'
  hRecoSVD.linecolor = 'black'
  hReco2D.linecolor = 'cyan'
  hJtFake.linecolor = 'yellow'
  rplt.hist(hJtMeas, axes = ax, label = 'Measured')
  rplt.hist(hJtTrue, axes = ax,label='True') #Plot jT histogram,
  rplt.hist(hRecoBayes, axes = ax, label = 'Bayes')
  rplt.hist(hRecoSVD, axes = ax, label = 'SVD')
  rplt.hist(hReco2D, axes = ax, label = '2D')
  rplt.hist(hJtFake, axes = ax, label = 'Fake')
  ax.legend(loc = 'lower left')

def drawZ(ax,hMeas,hTrue,hInput,hFake):
  #ax.set_xlim([0,10])
  #ax.set_ylim([1e-5,1000])
  ax.set_xlabel(r'$Z$')
  ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{dZ}$',fontsize=18)
  ax.set_yscale('log')
  #ax.set_xscale('log')

  hMeas.linecolor = 'blue'
  hTrue.linecolor = 'red'
  hInput.linecolor = 'pink'
  hFake.linecolor = 'yellow'
  rplt.hist(hInput,axes = ax, label = 'Input')
  rplt.hist(hMeas,axes = ax, label = 'Measured')
  rplt.hist(hTrue,axes = ax, label ='True') #Plot jT histogram,
  rplt.hist(hFake,axes = ax, label ='Fake')

  #ax.legend(loc = 'lower left')

def drawMulti(ax,hMeas,hTrue):
  #ax.set_ylim([1e-5,1000])
  ax.set_xlabel(r'$N_{tracks}$')
  ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{dN_{tracks}}$',fontsize=18)
  ax.set_yscale('log')
  #ax.set_xscale('log')

  hMeas.linecolor = 'blue'
  hTrue.linecolor = 'red'
  rplt.hist(hMeas,axes = ax, label = 'Measured')
  rplt.hist(hTrue,axes = ax, label ='True') #Plot jT histogram,
  ax.set_xlim([0,30])


def drawJetPtQA(ax,hMeas,hTrue,hInput,hReco):
  #ax.set_ylim([1e-5,1000])
  ax.set_xlabel(r'$p_{T}$')
  ax.set_ylabel(r'$\frac{dN}{dp_{jet}}$',fontsize=18)
  ax.set_yscale('log')
  ax.set_xscale('log')

  hMeas.linecolor = 'blue'
  hTrue.linecolor = 'red'
  hInput.linecolor = 'pink'
  hReco.linecolor = 'green'
  rplt.hist(hInput,axes = ax, label = 'Input')
  rplt.hist(hMeas,axes = ax, label = 'Measured')
  rplt.hist(hTrue,axes = ax, label ='True')
  rplt.hist(hReco, axes = ax, label = 'Unfolded')
  ax.set_xlim([5,150])


def drawMatchHisto(hists,jetPt,name,option='grid'):
  """Create an 4 by 2 grid of subfigures with shared axes and plots jT with background in jet pT bins
  Args:
    measJt: List of jT histograms
    jetPt: List of jet Pt bins in tuples (low border, high border)
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    name: Name of output file
  """
  if(option == 'grid'):
    fig, axs = plt.subplots(2,4,figsize=(10,5),sharey=True,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
    axs = axs.reshape(8) #Because the figures is in a 2x4 layout axs is a 2 dimensional array with 2x4 elements, this makes it a 1 dimensional array with 8 elements
    #axs[1].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
    for ax in [axs[0],axs[3],axs[4],axs[7]]: 
      ax.set_ylabel('Normalized',fontsize=12) #Add y-axis labels to left- and righmost subfigures
    for ax in axs[4:]:
      #ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
      ax.xaxis.set_ticks((0.5,1.5,2.5,3.5,4.5,5.5))
      ax.set_xticklabels(('No match','One match','More than 2 matches','No measured jT','Two matches','Fake tracks'),rotation = 'vertical')
    for ax in [axs[3],axs[7]]: 
      ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures
  
    for (h,ax,i,pT) in zip (hists,axs,range(0,8),jetPt):
      h.Scale(1.0/h.Integral(0,5))
      rplt.hist(h,xerr=False,emptybins=False,axes=ax)
      #rplt.errorbar(true,xerr=False,emptybins=False,axes=ax,label='True jT',fmt='go') 
      #if i == 0: #For the first subfigure add a legend to bottom left corner
      #  ax.legend(loc ='lower left')
      ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
      ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
      ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
      ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
      ax.text(3,0.5,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
      ax.set_xlim([0,6]) #Set x-axis limits
      ax.set_ylim([0,1]) #Set y-axis limits
      ax.grid(False) #Draw grid 
      ax.xaxis.set_ticks((0.5,1.5,2.5,3.5,4.5,5.5))
  
  
    plt.tight_layout()
    plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
    plt.savefig(name,format='pdf') #Save figure
    plt.show() #Draw figure on screen
  if(option == 'single'):
    ax2 = plt.gca()
    ax2.set_ylabel('Normalized')
    ax2.set_ylim([0,1])
    ax2.set_xlim([0,5])
    ax2.set_xticklabels(('No match','One match','More than 2 matches','No measured jT','Two matches','Fake tracks'),rotation = 'vertical')
    ax2.xaxis.set_ticks((0.5,1.5,2.5,3.5,4.5,5.5))
    ax2.yaxis.set_ticks_position('both')
    for h,pT,i,c in zip(hists,jetPt,range(0,5),('red','blue','green','yellow','cyan','orange','black','magenta')):
      h.color = c
      h.fillstyle = 'hollow'
      h.Scale(1.0/h.Integral(0,5))
      rplt.hist(h,xerr=False,emptybins = False,axes=ax2,histtype='bar',label = r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1]))
    ax2.legend(loc ='upper right')

    plt.savefig(name,format='pdf')
    plt.show()
    
def draw8grid(measJt,trueJt,jetPt,xlog = True,ylog = True,name="newfile.pdf",proj=None,unf2d = None, unf=None):
  """Create an 4 by 2 grid of subfigures with shared axes and plots jT with background in jet pT bins
  Args:
    measJt: List of jT histograms
    jetPt: List of jet Pt bins in tuples (low border, high border)
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    name: Name of output file
  """
  fig, axs = plt.subplots(2,4,figsize=(10,5),sharey=True,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  axs = axs.reshape(8) #Because the figures is in a 2x4 layout axs is a 2 dimensional array with 2x4 elements, this makes it a 1 dimensional array with 8 elements
  #axs[1].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  for ax in [axs[0],axs[3],axs[4],axs[7]]: 
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
  for ax in [axs[3],axs[7]]: 
    ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures

  for (jT,true,ax,i,pT) in zip (measJt,trueJt,axs,range(0,8),jetPt): 
    if(xlog):
      ax.set_xscale('log') #Set logarithmic scale
    if(ylog):
      ax.set_yscale('log')
    jT.SetMarkerColor('blue')
    true.SetMarkerColor('red')
    rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='Measured jT',fmt='+')
    rplt.errorbar(true,xerr=False,emptybins=False,axes=ax,label='True jT',fmt='go') 
    for h,color,title in zip((proj,unf,unf2d),(1,3,6),("Projected Meas","1D unfolded","2D unfolded")):
        if h is not None:
            h[i].SetMarkerColor(color)
            rplt.errorbar(h[i],xerr=False,emptybins=False,axes=ax,label=title,fmt='go')
    if i == 0: #For the first subfigure add a legend to bottom left corner
      ax.legend(loc ='lower left')
    ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
    ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    ax.set_xlim([0.01,5]) #Set x-axis limits
    ax.set_ylim([5e-4,2e3]) #Set y-axis limits
    ax.grid(True) #Draw grid 

  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  plt.savefig(name,format='pdf') #Save figure
  plt.show() #Draw figure on screen

def draw8gridcomparison(measJt,trueJt,jetPt,xlog = True,ylog = True,name="newfile.pdf",proj=None,unf2d = None, unf=None,fake=None, unf2dtest=None,leadingJt=None, **kwargs):
  """Create an 4 by 2 grid of subfigures with shared axes and plots jT with background in jet pT bins
  Args:
    measJt: List of jT histograms
    jetPt: List of jet Pt bins in tuples (low border, high border)
    xlog: Whether to use logarithmic scale on X-axis, default is True
    ylog: Whether to use logarithmic scale on Y-axis, default is True
    name: Name of output file
  """
  if 'start' in kwargs:
    start = kwargs.get('start')
  else:
    start = 0
  if 'stride' in kwargs:
    stride = kwargs.get('stride')
  else:
    stride = 2
  fig, axs = plt.subplots(2,4,figsize=(14,7),sharey=False,sharex=True) #Create figure with 8 subfigures, axs is a list of subfigures, fig is the whole thing
  axs = axs.reshape(8) #Because the figures is in a 2x4 layout axs is a 2 dimensional array with 2x4 elements, this makes it a 1 dimensional array with 8 elements
  #axs[1].text(0.02,0.005,r'pPb $\sqrt{s_{NN}} = 5.02 \mathrm{TeV}$' '\n Charged jT\n' r'Anti-$k_T$, R=0.4' '\nJet Cone',fontsize=7) #Add text to second subfigure, first parameters are coordinates in the drawn scale/units
  for ax in [axs[0],axs[3]]: 
    ax.set_ylabel(r'$\frac{1}{N_{jets}}\frac{dN}{j_{T}dj_{T}}$',fontsize=18) #Add y-axis labels to left- and righmost subfigures
  for ax in [axs[4],axs[7]]: 
    ax.set_ylabel('Ratio to truth',fontsize=12) #Add y-axis labels to left- and righmost subfigures
  for ax in axs[4:]:
    ax.set_xlabel(r'$j_{T}\left[GeV\right]$') #Add x-axis labels for bottom row
  for ax in [axs[3],axs[7]]: 
    ax.yaxis.set_label_position('right') #Set the y-axis label position to right hand side for the rightmost subfigures

  if(trueJt):
    divider = trueJt
  else:
    divider = measJt
  ratios = []
  for hists in (measJt,proj,unf2d,unf2dtest,unf,fake,leadingJt):
      if hists is not None:
          ratio = []
          for h,true in zip(hists,divider):
              h2 = h.Clone()
              h2.Divide(true)
              ratio.append(h2)
      else:
          ratio = None
      ratios.append(ratio)

  for (jT,ax,i,pT) in zip (measJt[start::stride],axs[0:4],range(start,8,stride),jetPt[start::stride]): 
    if(xlog):
      ax.set_xscale('log') #Set logarithmic scale
    if(ylog):
      ax.set_yscale('log')
    jT.SetMarkerColor('blue')
    #true.SetMarkerColor('red')
    if(leadingJt is None):
      rplt.errorbar(jT,xerr=False,emptybins=False,axes=ax,label='Measured jT',fmt='+')
    #rplt.errorbar(true,xerr=False,emptybins=False,axes=ax,label='True jT',fmt='go') 
    for h,color,title in zip((trueJt,proj,unf2d,unf2dtest,unf,fake,leadingJt),('red',1,3,8,6,7,9),("True jT","Projected Meas","2D unfolded","2D unfolded test","1D unfolded","Fakes","Leading ref.")):
        if h is not None:
            h[i].SetMarkerColor(color)
            rplt.errorbar(h[i],xerr=False,emptybins=False,axes=ax,label=title,fmt='go')
    if i == start: #For the first subfigure add a legend to bottom left corner
      ax.legend(loc ='lower left')
    ax.yaxis.set_ticks_position('both') #Show ticks on left and right side
    ax.xaxis.set_ticks_position('both')  #Show ticks on bottom and top
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.tick_params(which='both',direction='in') #Move ticks from outside to inside
    ax.text(0.3,1e2,r'$p_{{T,\mathrm{{jet}}}}$:''\n'r' {:02d}-{:02d} GeV'.format(pT[0],pT[1])) 
    if(leadingJt is None):
      ax.set_xlim([0.01,5]) #Set x-axis limits
    else:
      ax.set_xlim([0.01,20]) #Set x-axis limits
    maxContent = measJt[start].GetBinContent(measJt[start].GetMaximumBin())
    
    #ax.set_ylim([5e-4,2e4]) #Set y-axis limits
    ax.set_ylim([5e-8*maxContent,5*maxContent]) #Set y-axis limits
    ax.grid(True) #Draw grid 
  
  for ratio,color,title in zip(ratios,('blue',1,3,8,6,7,9),("Measured","Projected Meas","2D unfolded","2D unfolded test","1D unfolded","Fakes","Leading ref.")):
    if ratio is not None:
      for r,ax,i in zip(ratio[start::stride],axs[4:8],range(start,8,stride),):
        if(xlog):
          ax.set_xscale('log')
        r.SetMarkerColor(color)
        if(leadingJt is None or title != "Measured"):
          rplt.errorbar(r,xerr=False,emptybins=False,axes=ax,label=title,fmt='go')
        if(leadingJt is None):
          ax.set_xlim([0.01,5]) #Set x-axis limits
          ax.set_ylim([0,2]) #Set y-axis limits
        else:
          ax.set_xlim([0.01,20]) #Set x-axis limits
          ax.set_ylim([0,2]) #Set y-axis limits
        ax.grid(True)
      
  plt.tight_layout()
  plt.subplots_adjust(wspace =0,hspace=0) #Set space between subfigures to 0
  plt.savefig(name,format='pdf') #Save figure
  plt.show() #Draw figure on screen

