# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 18:48:43 2016

@author: William
"""

import sys
import numpy as np
import scipy as sp
from scipy import linalg
from scipy import special
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

from osgeo import ogr, osr, gdal

from scipy.stats import chi2
import cv2
import time

from skimage.segmentation import random_walker
from skimage import morphology
import maxflow

#private libraries
import GeoFileHandle as gfh


class ChangeDetection:
#==============================================================================
# Change detection Class to perform the IRMAD 
#==============================================================================
    def __init__(self,dirpath,xfn,yfn):
#==============================================================================
# init change detection instance
#==============================================================================
        #define file and path names
        self.dirpath = dirpath
        self.outputpath = dirpath + 'output/'
        self.png_path = dirpath + 'PNG/'
        self.png_fn_tag = '.png'
        self.shape_path = dirpath + 'shape/'
        self.Xfn = dirpath + xfn
        self.Yfn = dirpath + yfn
        
        
        self.ROIfn = 'ROI.tif'
        self.XCVfn = 'XCV.tif'
        self.YCVfn = 'YCV.tif'
        self.MADfn = 'MAD.tif'
        self.whtfn = 'wht.tif'
        self.chi2fn = 'chi2.tif'
        self.chngfn = 'chng.tif'
        
        self.itermax = 1
        
        #open image files
        self.X = gfh.GetImageInfo(self.Xfn)
        self.Y = gfh.GetImageInfo(self.Yfn)
        
        try:
            temp = gfh.GetImageInfo(self.Xfn)
        except:
            print 'Error, file could not be loaded'
            
        self.nrow, self.ncol, self.nvar = gfh.GetImageSize(temp)
        
        self.ROI_toogle = False
        self.ROI = None
        
        self.Ndata = self.nrow*self.ncol

        #output file initialization
        self.XCV = gfh.CreateFileAs(self.outputpath + self.XCVfn, temp)
        self.YCV = gfh.CreateFileAs(self.outputpath + self.YCVfn, temp) 
        self.MAD = gfh.CreateFileAs(self.outputpath + self.MADfn, temp)
        self.chi2 = gfh.CreateFileAs(self.outputpath + self.chi2fn, temp, Nband = 1)

        whtArray = np.ones((self.nrow,self.ncol,1))
        self.wht  = gfh.CreateFileAs(self.outputpath + self.whtfn, temp, Nband = 1,data = whtArray)
        
        self.chng = gfh.CreateFileAs(self.outputpath + self.chngfn, temp, Nband = 1)
        
#==============================================================================
# Path functions 
#==============================================================================
    
    def SetOutputPath(self,filename):
        self.outputpath = filename  
        return 0
    
    def SetPNGPath(self,filename):
        self.png_path = filename
        return 0
    
    def SHPpath(self,filename):
        self.shape_path = filename
        return 0

#==============================================================================
# Create ROI from shape file
#==============================================================================
    def CreateROIfromSHP(self,shpfn):
        binary = gfh.polygon2Binary(shpfn,self.Xfn)
        
        ROIimg = binary[:,:,np.newaxis]
        
        self.ROI = gfh.CreateFileAs(self.outputpath +self.ROIfn,self.X,
                                    Nband = 1,data = ROIimg)
        self.ROI_toogle = True        
        
        return 0

        
    def SetROIfile(self,filepath):
        self.ROI_toogle = True
        self.ROI = gfh.GetImageInfo('ROI.tif')
        return 0
        
#==============================================================================
# IRMAD
#==============================================================================

    def IRMAD(self,itermax,**warg):
        self.itermax = itermax
#------------------------------------------------------------------------------

        rho = np.zeros((itermax,self.nvar))   
        
#------------------------------------------------------------------------------
#   IRMAD iterations
#------------------------------------------------------------------------------
    
        for n in range(0, itermax):
            print "Iteration: ", n+1
            
    #        if (n+1) == itermax:
    #            whtdata = gfh.ReadData(whtfile)
    #            whtdata = (whtdata > 0.5).astype(np.float64)
    #            ROIdata = gfh.ReadData(ROIfile)
    #            ix,iy,iz = np.where(ROIdata == 0)
    #            whtdata[ix,iy] = -9999
    #            gfh.WriteToFile(whtfile,whtdata)
            
            #cumulative cov determ. from imagefiles (filepaths as input)
            D, meanw = self.PrvCovW()
            
            meanX, meanY = np.split(meanw,2)    
            
            rho[n,:],varmads, A1, A2 = self.CanCorrAnalysis(D)
            
            print ["%.4f" % number for number in rho[n,:]]
            
            #loop over rows            
            for irow in range(0,self.nrow):
                Xr = gfh.ReadImgRow(self.X,irow)
                Yr = gfh.ReadImgRow(self.Y,irow)
        
                MADr, chi2r, whtr, XCVr, YCVr = self.Eig2Chi2Weights(Xr,Yr,
                                                A1,A2,meanX,meanY,varmads)
          
                
                if self.ROI_toogle == True:
                    #find indicies in ROI file where not selected
                    ROIRow = gfh.ReadImgRow(self.ROI,irow)
                    _,ind = np.where(ROIRow == 0)
                else:
                    ind = np.ones(1,self.ncol)
                    
                whtr[:,ind] = -9999 #no data value
                gfh.WriteImgRow(self.wht,whtr,irow)
                #Save Row to files
                
                #if last iteration? - write to MAD and chi2 file
                if (n+1) == itermax:
                    
                    MADr[:,ind] = -9999 #no data value
                    gfh.WriteImgRow(self.MAD,MADr,irow)
                    
                    chi2r[:,ind] = -9999 #no data value
                    gfh.WriteImgRow(self.chi2,chi2r,irow)
                    
                    XCVr[:,ind] = -9999 #no data value
                    gfh.WriteImgRow(self.XCV,XCVr,irow)
                    #Save Row to files
                    
                    YCVr[:,ind] = -9999 #no data value
                    gfh.WriteImgRow(self.YCV,YCVr,irow)
                    #Save Row to files                
        
        print "IRMAD complete"
        return rho, A1, A2

    
#==============================================================================
# Provisional mean and covariance
#==============================================================================
    def PrvCovW(self):
        #calculate the mean row by row
        mean = np.zeros((self.nvar*2,1))
        wsum = np.float64(0.0)
        SoS = np.zeros((self.nvar*2,self.nvar*2))
        ROIRow = 0
        for R in range(self.nrow):
            
            X = gfh.ReadImgRow(self.X,R)
            Y = gfh.ReadImgRow(self.Y,R)
            x = np.vstack((X,Y))
            
            WhtRow = gfh.ReadImgRow(self.wht,R)
            
            if self.ROI_toogle == True:
                ROIRow = gfh.ReadImgRow(self.ROI,R)
                _,iROI = np.where(ROIRow == 1)
                WhtRow = WhtRow[:,iROI]
                x = x[:,iROI]
            else:
                ROIRow = np.ones(1,self.ncol)

            if np.size(WhtRow) != 0:
                if np.sum(WhtRow) != 0:
                    d = WhtRow*(x - mean)
                    wsum += np.sum(WhtRow)
                    mean += np.sum(d,axis = 1)[:,np.newaxis]/wsum
                    SoS += np.mat(d)*np.mat((x - mean)).T
            
        D = np.mat(SoS/((self.Ndata-1)*wsum/self.Ndata))
                
        return D, mean
    
#==============================================================================
# Canonical correlation analysis    
#==============================================================================
    def CanCorrAnalysis(self,D):
        s11, s12, s21, s22 = splitDispersion(D)    
        
        #generalized eigenvalueproblem
        Lmbd1, A1 = sp.linalg.eig(s12*s22.I*s21,s11)
       
        #sort eigenvalues
        I1 = np.argsort(Lmbd1)
      
        #remove 0j (imag) from data
        Lmbd1 = np.real(Lmbd1[I1])
       
        #normalize to unit variance
        A1 = UnitVariance(A1,s11)
        
        #sort eigenvectors to CVX
        A1 = A1[:,I1]
        
        #sum of correlations between X and CV(X) positive
        invstderr = np.mat(np.diag(1/np.sqrt(np.diag(s11))))
        sgn1 = np.sign(np.sum(invstderr*s11*A1,axis = 0))
        A1 = A1*np.diag(np.squeeze(np.array(sgn1)))
        
        #determ. eigenvectors to CVY
        A2 = s22.I*s21*A1;
        
        #Unit variance
        A2 = UnitVariance(A2,s22)  
        
        #canonical correlations
        rho = np.sqrt(Lmbd1)
        
        #determ. Variance of the mads
        varmads = 2*(1-rho)[:,np.newaxis]
        
        
        return rho,varmads,A1,A2
    
    
#==============================================================================
# Weight calculations    
#==============================================================================
    def Eig2Chi2Weights(self,Xdata,Ydata,A1,A2,meanX,meanY,varmads):
        #dot image with eigenvector to get the mad images
        XCV = np.dot(A1.T,(Xdata-meanX))
        YCV = np.dot(A2.T,(Ydata-meanY))
        
        nvar = len(meanX)
    
        #diffence images
        MAD = XCV-YCV
    
        chi2 = np.sum(np.square(MAD)/varmads,axis = 0)
        
        wht = 1-sp.special.gammainc(0.5*nvar,0.5*chi2)
        
        return MAD,chi2,wht, XCV, YCV
        
        
        




   

#==============================================================================
# PRINTS and PLOTS begin from here::
#==============================================================================







#==============================================================================
#  Print files to PNG
#==============================================================================
    def printPNGfiles(self,pngfn):
        gfh.MergeAsGrayscale(self.MAD,'MAD_SbS_.png',4)
        gfh.TIFF2PNG(self.MAD,'MAD' + self.png_fn_tag,'mult4add128')
        gfh.TIFF2PNG(self.chi2,'chi' + self.png_fn_tag,'Fit')
        gfh.TIFF2PNG(self.wht,'wht' + self.png_fn_tag,'Fit')
        
        
        #gfh.TIFF2PNG(Xfile,'PNG/xfile.png','noFit')
        #gfh.TIFF2PNG(Yfile,'PNG/yfile.png','noFit')

    
    def plotRhoIterations(self,rho):
        if rho.shape[0] > 1:
            plt.figure(1)
            for i in range(0,self.nvar):
                plt.plot(range(1,self.iter+1),rho[:,i])
            #plt.axis([1, it, 0, 1])
            plt.xlabel('iterations')
            plt.title('Canonical Correlations (rho) development over iterations')
            
            plt.show()
            #rhofn = 'PNG/rho' + pngfn + '.png'
            #plt.savefig('figure.png', format='png')


#==============================================================================
# 
#==============================================================================

    def PlotMADpdf(self):      
        img = gfh.ReadData(self.MAD)
        
        colorList = ('r','g','b','y','c','m')
        color = colorList[0:self.nvar]
            
        #color = ('k','k','k','k')
        hmin = -30  #chi2.ppf(0.0001, df)
        hmax = 30 #chi2.ppf(0.9999, df)
        bins = 100
        
        ##rectancle select
        #rmin = 0
        #rmax = 500
        #cmin = 200
        #cmax = 1000
        #
        #maskrect = np.zeros((img.shape[:2]), np.uint8)[:,:,np.newaxis]
        #maskrect[rmin:rmax, cmin:cmax] = 1
    
    
        wht = gfh.ReadData(self.wht)
        
        #wheight select
        mask  = wht<0.1
        
        if self.ROI_toogle == True:
            maskrect = gfh.ReadData(self.ROI)
            mask = mask*maskrect
            
        mask = mask.astype(np.uint8)
        
        plt.figure(2)
        histrUnNorm = np.zeros((bins,len(color)))
        x = np.linspace(hmin,hmax, bins)
        for i,col in enumerate(color):
            histr = cv2.calcHist([img],[i],mask,[bins],[hmin,hmax])
            histrUnNorm[:,i] = np.squeeze(histr)
            histr /= sp.integrate.trapz(np.squeeze(histr),x = x)
            
            maskdata = mask.reshape(img.shape[0]*img.shape[1],1)
            index = np.where(maskdata == 1)
            imgx = img[:,:,i].reshape(img.shape[0]*img.shape[1],1)
            std = np.std(imgx[index])
            mean = np.mean(imgx[index])
            gauss = sp.stats.norm.pdf(x,mean,std)
            
            
            plt.subplot(2,2,i)
            plt.plot(x,histr,color = col)
            plt.plot(x,gauss,color = col)
            
        plt.show()
        
        plt.figure(3)
        for i,col in enumerate(color):
            plt.plot(x,histrUnNorm[:,i],color = col)
        plt.show()
        
#==============================================================================
# Plot chi2
#==============================================================================
    
    def PlotChi2(self):    
        img = gfh.ReadData(self.chi2)
        hmin = 0#chi2.ppf(0.0001, df)
        hmax = 1000#chi2.ppf(0.9999, df)
        bins = 30000
        
        x = np.linspace(hmin,hmax, bins)
        
        
        chng = gfh.ReadData(self.chng)
        nochng = (1 - chng).astype(np.uint8)
        
        histr = cv2.calcHist([img],[0],nochng,[bins],[hmin,hmax])
        histr /= sp.integrate.trapz(np.squeeze(histr),x = x)
        
        plt.figure(4)
        plt.plot(x, chi2.pdf(x, self.nvar),'r-', lw=5, alpha=0.6, label='chi2 pdf')
        plt.plot(x,histr, label = 'chi2 histogram')
        plt.show()
        plt.legend()
    
#==============================================================================
# create scatterplot
#==============================================================================
    def MakeWeightetHistogram(self):
        
        cmaplist =(plt.cm.Reds,plt.cm.Greens,plt.cm.Blues)
        Titles = ('Red','Green','Blue','DSM','other')
        
        DSMmin = 5
        DSMmax = 20
         
        Hcum = np.zeros((255,255,self.nvar))
        H = Hcum
        for irow in range(0,self.nrow):
            RowX = gfh.ReadImgRow(self.X,irow)
            RowY = gfh.ReadImgRow(self.Y,irow)
            
              
            #find indicies in select file where not selected
            whtRow = gfh.ReadImgRow(self.wht,irow).squeeze()
            
            if self.ROI_toogle == True:
                ROIRow = gfh.ReadImgRow(self.ROI,irow)
                _,ind = np.where(ROIRow == 1)
            else:
                ind = np.ones((1,self.ncol))
            
            if len(ind)>0:
                for ivar in range(self.nvar):
                    if Titles[ivar] == 'DSM':
                        bins = (255,255)
                        Irange = [[DSMmin,DSMmax],[DSMmin,DSMmax]]
                    else:
                        bins = (255,255)
                        Irange = [[0,255],[0,255]]                    
                    H[:,:,ivar],x,y = np.histogram2d(RowX[ivar,ind],RowY[ivar,ind],
                                      bins=bins,range=Irange,weights = whtRow[ind])    
                        
                Hcum = Hcum + H
        
        if self.nvar > 3:
            l = np.ceil(np.sqrt(self.nvar))
            h = np.ceil(np.float(self.nvar)/l)
        else:
            h = 1
            l = self.nvar
    
        plt.figure(num=5, figsize=(6*l, 6*h +1))
    
        for i in range(self.nvar):
            plt.subplot(h,l,i+1)
               
            if i > 2:
                cmap = plt.cm.Greys
            else:
                cmap = cmaplist[i]
            
            if Titles[i] == 'DSM':
                extent = [DSMmin,DSMmax,DSMmin,DSMmax]
                axlabel = 'height'
            else: 
                axlabel = 'pixel intensity'
                extent = [0,255,0,255]
            
            Hcumi = Scale2uint8(Hcum[:,:,i])
            plt.imshow(np.flipud(Hcumi),cmap=cmap, interpolation='none', 
                       extent=extent)
            plt.title(Titles[i])
            plt.xlabel(Titles[i] +' '+ axlabel + ' (Image 2)')
            plt.ylabel(Titles[i] +' '+ axlabel + ' (Image 1)')
        
        plt.suptitle("Weightet no-change Histogram (10 iterations)", size=20)
        #plt.tight_layout()
        plt.show()
        
        plt.savefig('WeightetHistogram/WH_' + str(self.nvar) + 'B.png')
        
        return Hcum
    
#==============================================================================
# MAD correlated with the input paramenters
#==============================================================================
    def MADCorrInput(self):
    
        D1,_ = self.PrvCovW(self.X,self.MAD)
        D2,_ = self.PrvCovW(self.Y,self.MAD)
        df = self.nvar
        L1 = np.diag(D1[0:df,-df:])
        L2 = np.diag(D2[0:df,-df:])
        print L1,L2
        return np.hstack((L1,L2))
    

        
#==============================================================================
# Select raster by size
#==============================================================================
    def FilterChangeByArea(self,area):
        size = np.round(area / (0.4**2)) #devide with square of image resolution
        
        
        
        BinArray = gfh.ReadData(self.chng)
        
        BigChng = np.zeros(np.shape(BinArray))
        
        LabelArray = morphology.label(np.squeeze(BinArray))        
        
        N_chng = np.max(LabelArray)
        if N_chng > 0:
            k = 1
            for i in range(N_chng):
                indx,indy = np.where(LabelArray == i)    
                #larger than x m² and smaller than half the image
                if (len(indx) >= size) & (len(indx) < (self.Ndata/2)): 
                    BigChng[indx,indy] = k
                    k += 1     
            
            print str(k)
            ROI = gfh.ReadData(self.ROI)
            xROI,yROI,_ = np.where(ROI == 0)
            BigChng[xROI,yROI] = -9999
            
            gfh.WriteToFile(self.chng,BigChng)
            
            gfh.TIFF2PNG(self.chng,self.png_path + 'chng_AreaSelect' +
            self.png_fn_tag,'Fit')
        else:
            print 'no changes found, test thresholdvalues!'
        return 0
        
#==============================================================================
# object distance to shape file
#==============================================================================
    def ChangeDistanceToVector(self,vectorfn):
        bufferDist = 1000
        
        shpChng = gfh.RasterToPolygon(self.chng,'ChngPolygon.shp')
        
        dist = gfh.MinDistanceToSHP(shpChng,vectorfn,None)
        
        index = np.array(np.where(dist>1))
        
        #create outputfile
        outputfn = 'shape/ChngSel.shp'
        
        gfh.SelectVectorObjects(index,shpChng,outputfn)
        
#==============================================================================
# Change geometry filter 
#==============================================================================
    def ChangeGeometryfilter(self,AreaMin,distanceTo,vector):
        shpChng = gfh.RasterToPolygon(self.chng,self.shape_path + 'ChngPolygon.shp')
        
        driver = ogr.GetDriverByName('ESRI Shapefile')
        
        #open skel shpfile
        Chngdata = driver.Open(shpChng,1) # 0 means read-only. 1 means writeable.      
        Chnglayer = Chngdata.GetLayer()
        ChngFeatCount = Chnglayer.GetFeatureCount()
        
        #open input vectorfile
        inputdata = driver.Open(vector,0) # 0 means read-only. 1 means writeable.      
        inputlayer = inputdata.GetLayer()            
        inputFeatCount = inputlayer.GetFeatureCount()
        if ChngFeatCount > 0:
            #loop for each Change feature
            for i in range(ChngFeatCount):    
                Chngfeat = Chnglayer.GetFeature(i)
                ChngGeom = Chngfeat.geometry()
                Area = ogr.Geometry.Area(ChngGeom)
                if (Area <= AreaMin) | (Area >= 1000):
                    featId = Chngfeat.GetFID()
                    Chnglayer.DeleteFeature(featId) 
                
                elif not (distanceTo == None): #loop for each Change feature
                    for j in range(inputFeatCount):    
                        inputfeat = inputlayer.GetFeature(j)
                        inputGeom = inputfeat.geometry()
                        
                        dist = ogr.Geometry.Distance(inputGeom,ChngGeom)
                        #save the smallest one        
                        if dist < distanceTo:
                            featId = Chngfeat.GetFID()
                            Chnglayer.DeleteFeature(featId) 
                            break
                    

            
            gfh.TIFF2PNG(self.chng,self.png_path + 'chng_AreaSelect' +
            self.png_fn_tag,'Fit')
        else:
            print 'no changes found, test thresholdvalues!'
            
        
        return 0

#==============================================================================
#  initialize CD interactive
#==============================================================================
    def Chi2ChangeDetectionInteractive(self):
        chi2img = gfh.ReadData(self.chi2).squeeze()
        
        methodfunc = InteractiveFigure(chi2img)
        
        methodfunc.slider1.on_changed(methodfunc.updateImg)
        methodfunc.slider2.on_changed(methodfunc.updateImg)
        methodfunc.radio.on_clicked(methodfunc.ChangeMethod)
        
        plt.show()  

#==============================================================================
# Threshold and create change image
#==============================================================================
    def ChiThreshold(self,threshold):
        chi2data = gfh.ReadData(self.chi2)
        
        Change = thrhold(chi2data,threshold)
        
        gfh.WriteToFile(self.chng,Change)
        
        gfh.TIFF2PNG(self.chng,self.png_path + 'chng_thr' +
        self.png_fn_tag,'Fit')
        return 0
        
    def ChiMRF(self,threshold,beta):
        chi2data = gfh.ReadData(self.chi2)
        
        Change = MRF(chi2data,threshold,beta)
        
        gfh.WriteToFile(self.chng,Change)
        
        gfh.TIFF2PNG(self.chng,self.png_path + 'chng_MRF' +
        self.png_fn_tag,'Fit')
        
    def BinMathMorphOpen(self,threshold,size):
        chi2data = gfh.ReadData(self.chi2)
        
        Change = BinMM(chi2data,threshold,size)
        
        gfh.WriteToFile(self.chng,Change)
        
        gfh.TIFF2PNG(self.chng,self.png_path + 'chng_BMM' +
        self.png_fn_tag,'Fit')
        return 0
        
    def GrayMathMorphOpen(self,threshold,size):
        chi2data = gfh.ReadData(self.chi2)
        
        Change = GrayMM(chi2data,threshold,size)
        
        gfh.WriteToFile(self.chng,Change)
        
        gfh.TIFF2PNG(self.chng,self.png_path + 'chng_GMM' +
        self.png_fn_tag,'Fit')
        return 0

#------------------------------------------------------------------------------       
def thrhold(img,thresh):
    img = img > thresh
    return img
    
def MRF(img,thresh,beta):
    #beta = 10**(self.param2-6)
    E_max = 10
    alpha = 0.0000001

    #intersect of energyfunctions will be 1 independent on alpha
    mu_s = thresh +(1/alpha)**0.5      
    mu_t = thresh -(1/alpha)**0.5  
    
    def SourceEnergy(img,alpha,beta,mu,E_max):
        #source - change
        E_s = alpha*(mu_s - img)**2
        E_s[E_s > E_max] = E_max
        #no penalty with values higher than mu_s
        E_s[img > mu_s] = 0
        return E_s
    def SinkEnergy(img,alpha,beta,mu,E_max):
        #sink no-change     
        E_t = alpha*(mu_t - img)**2
        #maximum penalty
        E_t[E_t > E_max] = E_max
        #no penalty with values lower than mu_t
        E_t[img < mu_t] = 0
        return E_t
    
    E_s = SourceEnergy(img,alpha,beta,mu_s,E_max)
    
    E_t = SinkEnergy(img,alpha,beta,mu_t,E_max)
            
    if 0: #make the plot
        x = np.arange(10000)
        
        E_sx = SourceEnergy(x,alpha,beta,mu_s,E_max)
    
        E_tx = SinkEnergy(x,alpha,beta,mu_t,E_max)
        
        plt.figure(11)
        plt.plot(x,E_tx,x,E_sx)
        plt.show()
        
    
    # Create the graph.
    g = maxflow.Graph[float]()
    # Add the nodes. nodeids has the identifiers of the nodes in the grid.
    nodeids = g.add_grid_nodes(img.shape)
    # Add non-terminal edges with the same capacity.
    g.add_grid_edges(nodeids, beta)
    # Add the terminal edges. The image pixels are the capacities
    # of the edges from the source node. The inverted image pixels
    # are the capacities of the edges to the sink node.
    g.add_grid_tedges(nodeids, E_t, E_s)
    
    # Find the maximum flow.
    g.maxflow()
    # Get the segments of the nodes in the grid.
    sgm = g.get_grid_segments(nodeids)
    
    # The labels should be 1 where sgm is False and 0 otherwise.
    imgout = np.int_(np.logical_not(sgm))
    
    return imgout


def BinMM(img,thresh,size):
    #threshold at sliderval (param1)
    img = thrhold(img,thresh)
    
    structElement = morphology.rectangle(size,size)
    #structElement = morphology.disk(size)
    #structElement = morphology.diamond(size)

    img = morphology.binary_opening(np.squeeze(img), structElement).astype(np.uint8)
    return img[:,:,np.newaxis]

def GrayMM(img,thresh,size):
    structElement = morphology.rectangle(size,size)
  
    img[img == -9999] = 0
    img = img/5000
    img[img < 0] = 0
    img[img > 1] = 1
    
    outdata = morphology.opening(img, structElement)
    #threshold after bin
    imgOut = outdata > 255*thresh/5000
    
    return imgOut
 
#------------------------------------------------------------------------------

def Scale2uint8(I):

    I_out = (((I - I.min()) / (I.max() - I.min())) * 255.9)

    return I_out.astype(np.uint8)

def UnitVariance(A,sii):

    A_unit = np.array(A)*np.array(1/np.sqrt(np.diag(A.T*sii*A)))[np.newaxis,:]
    
    return np.mat(A_unit)
#==============================================================================
# split dispersion into s11 s12 s21 s22    
#==============================================================================
    
def splitDispersion(D):
    s1, s2 = np.vsplit(D,2)
    s11, s12 = np.hsplit(s1,2)
    s21, s22 = np.hsplit(s2,2)

    return s11, s12, s21, s22
    
#==============================================================================
# IRMAD interactive class
#==============================================================================

class InteractiveFigure:
    def __init__(self,image):
        self.Mlist = ('threshold', 'MRF', 'Mathematical Morphology', 'gray MM')
        self.method = self.Mlist[0]
        
        self.img = image
        
        self.nvar = 3
        print 'edit'
        self.param1 = chi2.ppf(0.995, self.nvar)
        img_init = self.img > self.param1
        
        self.param2 = 2
        
        self.fig, ax = plt.subplots(num=20, figsize=(30, 20))
        self.im1 = ax.imshow(img_init, cmap = plt.cm.Greys_r)
        plt.subplots_adjust(left=0.3)
        
        axcolor = 'lightgoldenrodyellow'
        rax = plt.axes([0.05, 0.7, 0.15, 0.15], axisbg=axcolor)
        self.radio = RadioButtons(rax, self.Mlist)
        
        axS1 = self.fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        self.slider1 = Slider(axS1, 'Threshold', 1, 5000, valinit=self.param1)
            
        axS2 = self.fig.add_axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
        self.slider2 = Slider(axS2, 'Smoothnes', 1, 11, valinit=self.param2)
           
        
    def ChangeMethod(self,label):
        self.method = label
        self.updateImg(0)

#------------------------------------------------------------------------------
    
    def updateImg(self,val):
        if val == 0:
            del val
        else:
            self.param1 = self.slider1.val
            self.param2 = self.slider2.val
            
        
        if self.method == self.Mlist[0]:
            img = thrhold(self.img,self.param1)
            
        elif self.method == self.Mlist[1]: 
            img = MRF(self.img,self.param1,self.param2)
            
        elif self.method == self.Mlist[2]:
            img = BinMM(self.img,self.param1,self.param2)    
            
        elif self.method == self.Mlist[3]:
            img = GrayMM(self.img,self.param1,self.param2)
        else:
            print "error"
        
            
        self.im1.set_data(img)
        self.fig.canvas.draw()

 