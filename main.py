# -*- coding: utf-8 -*-
"""
Created on Sat Jan 23 13:10:29 2016

@author: William
"""

import IRMAD as im

#==============================================================================
# GDAL install Error handler
#==============================================================================
#gdal.PushErrorHandler(gfh.gdal_error_handler)
    
#==============================================================================
# Load images
#==============================================================================    
inputpath = 'Input/'

#inputpath = '../Input/mergedDataset/VirtualSpectral07/'

#Xfile = inputpath + 'mergedDataset/VirtualSpectral07/Spectral14LRedit.tif'
#Yfile = inputpath + 'mergedDataset/VirtualSpectral07/Spectral14LR.tif'

#Xfile = inputpath + 'Praestemosen/PM_RGB_14_400mm_GST.tif'
#Yfile = inputpath + 'Praestemosen/PM_RGB_16_400mm_NIRAS.tif'

#Xfile = inputpath + 'Praestemosen/PM_RGB_14_100mm_GST.tif'
#Yfile = inputpath + 'Praestemosen/PM_RGB_16_100mm_NIRAS.tif'

#------------------------------------------------------------------------------
#initialize the combined IRMAD method
#------------------------------------------------------------------------------
Xfile = inputpath + 'Praestemosen/PM_COM_14_400mm.tif'
Yfile = inputpath + 'Praestemosen/PM_COM_16_400mm.tif'

dirpath = 'CD_Combined/'

CDCombined = im.ChangeDetection(dirpath,Xfile,Yfile)

#------------------------------------------------------------------------------
#initialize the binary and/or method
#------------------------------------------------------------------------------
Xfile = inputpath + 'Praestemosen/PM_RGB_14_400mm_GST.tif'
Yfile = inputpath + 'Praestemosen/PM_RGB_16_400mm_NIRAS.tif'

dirpath = 'CD_BinaryAndOr/'

CDBinaryAndOR = im.ChangeDetection(dirpath,Xfile,Yfile)

#==============================================================================
# Create DSM png files
#==============================================================================
#XDSM = inputpath + 'Praestemosen/PM_DSM_14_400mm_GST2.tif'
#YDSM = inputpath + 'Praestemosen/PM_DSM_16_400mm_NIRAS2.tif'

#gfh.TIFF2PNG(XDSM,'PNG/xDSM.png','Fit')
#gfh.TIFF2PNG(YDSM,'PNG/yDSM.png','Fit')

#==============================================================================
# convert ROI polygon to UTM
#==============================================================================
#shapein = inputpath + 'SelectPolygon/Praestemosen2.shp'
#shpfn = inputpath + 'SelectPolygon/Praestemosen2UTM.shp'
#gfh.ShapeProjConvert(shapein,'UTM32',shpfn)

#==============================================================================
# Load ROI polygon and create binary ROI image
#==============================================================================
shpfn = inputpath + 'SelectPolygon/Praestemosen2UTM.shp'
CDCombined.CreateROIfromSHP(shpfn)
CDBinaryAndOR.CreateROIfromSHP(shpfn)
#==============================================================================
# Run IRMAD
#==============================================================================

it = 10

CDCombined.IRMAD(it,converg = 0.01)

CDBinaryAndOR.IRMAD(it, converg = 0.01)

#==============================================================================
#  Generate IRMAD PNG-files
#==============================================================================
#CD_Combined
CDCombined.printPNGfiles()

#CDCombined.plotRhoIterations()

#CDCombined.PlotEigenvectors()

#CD_BinaryAndOr
CDBinaryAndOR.printPNGfiles()

#CDBinaryAndOR.plotRhoIterations()

#CDBinaryAndOR.PlotEigenvectors()

#==============================================================================
# change image using threshold
#==============================================================================
#CD_Combined
threshold = 4000
structsize = 2
beta = 10

CDCombined.ChiThreshold(threshold)
CDCombined.MRFChange(threshold, beta)
CDCombined.GrayMathMorphOpen(threshold,structsize)

CDCombined.Chi2ChangeDetectionInteractive()


#CD_BinaryAndOr
threshold = 16.266236196238129 # chi² (0.999,df = 3)
threshold = 20
structsize = 2
beta = 10

CDBinaryAndOR.ChiThreshold(threshold)
CDBinaryAndOR.MRFChange(threshold, beta)
CDBinaryAndOR.GrayMathMorphOpen(threshold,structsize)

#CDBinaryAndOR.Chi2ChangeDetectionInteractive()


#==============================================================================
# Include DSM in CD_BinaryAndOr
#==============================================================================

DSM14 = inputpath + 'Praestemosen/PM_DSM_14_400mm_GST.tif'
DSM16 = inputpath + 'Praestemosen/PM_DSM_16_400mm_NIRAS.tif'

CDBinaryAndOR.DSMCD(DSM14,DSM16)
CDBinaryAndOR.DSMdifference(DSM14,DSM16)
CDBinaryAndOR.DSMdiffthresh()
#CDBinaryAndOR.Boolean('AND')
#CDBinaryAndOR.Boolean('OR')


#==============================================================================
# Filter objects by area 
#==============================================================================
Area = 50
distanceTo = None
shape = 'combined_CD/Vector_input_ex/skel/skel.shp'

CDCombined.ChangeGeometryfilter(Area,distanceTo,shape)

#==============================================================================
#  """different plots"""
#==============================================================================
#plot MAD Gaussian PDF
#CDCombined.PlotMADpdf()

# chi-squared PDF
#CDCombined.PlotChi2()

#==============================================================================
# create histograms of input and Canonical variates
#==============================================================================
#H = CDCombined.WeightetInputHistogram()

#H = CDBinaryAndOR.WeightetInputHistogram()

#H2 = CDCombined.CanonicalVariatesHistogram()

#==============================================================================
# Correlation MAD with Original vars
#==============================================================================
#Cx = CDCombined.MADCorrInput()


#==============================================================================
# Create PNG spippets for report 
#==============================================================================
#R = CDCombined.Cov2Corr(1)


#==============================================================================
# generate matrices
#==============================================================================
#CDCombined.InputAndCVCorr()

#CDCombined.InputAndMADCorr()

"""

path = '/home/grusinator/Dropbox/DTU/Speciale/Programming/IRMAD/Input/Cutouts/'
nr = '01'
im.TIFF2PNG(path + 'RGB14-' + nr + '.tif',
            path + 'PNG/RGB14-' + nr + '.png','Fit')

im.TIFF2PNG(path + 'RGB16-' + nr + '.tif',
            path + 'PNG/RGB16-' + nr + '.png','Fit')
            
im.TIFF2PNG(path + 'DSM14-' + nr + '.tif',
            path + 'PNG/DSM14-' + nr + '.png','Fit')

im.TIFF2PNG(path + 'DSM16-' + nr + '.tif',
            path + 'PNG/DSM16-' + nr + '.png','Fit')

"""

#------------------------------------------------------------------------------
print 'Main script done!'
#------------------------------------------------------------------------------











