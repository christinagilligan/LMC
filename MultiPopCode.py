import os, sys, re
import csv
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import binned_statistic
from scipy import interpolate
from scipy import spatial
from scipy.interpolate import interp1d
from numpy import zeros
from numpy import extract
from sklearn.mixture import GMM
from sys import argv
from scipy.optimize import leastsq

#set the paramters for running functions over the three filters and color combos
def settingParameters(filter):
    print(str(filter))
    column=0
    color=0
    colorColumn=0
    axis=0
    FieldBlue=0
    FieldRed=0
    MSmag=0
    MSmagend=0
    ClusterReddening=0
    Startmag=0
    Endmag=0
    if filter==filter1:
        column=column1
        colorColumn=column2
        color=str(filter1)+'-'+str(filter2)
        axis=axis1
        FieldBlue=BlueEdge1
        FieldRed=RedEdge1
        MSmag=MSmag1
        MSmagend=MSmagend1
        ClusterReddening=ClusterReddening1*1.281
        Startmag = MSmagend2
        Endmag = MSmag2
        print('Filter1')
        
    if filter==filter2:
        column=column2
        colorColumn=column3
        color=str(filter2)+'-'+str(filter3)
        axis=axis2
        FieldBlue=BlueEdge2
        FieldRed=RedEdge2
        MSmag=MSmag2
        MSmagend=MSmagend2
        ClusterReddening=ClusterReddening2*1.281
        Startmag = MSmagend2
        Endmag = MSmag2
        print('Filter2')


    if filter==filter3:
        column=column1
        colorColumn=column3
        color=str(filter1)+'-'+str(filter3)
        axis=axis3
        FieldBlue=BlueEdge3
        FieldRed=RedEdge3
        MSmag=MSmag1
        MSmagend=MSmagend1
        ClusterReddening=ClusterReddening3*1.281
        Startmag = MSmagend2
        Endmag = MSmag2
        print('Filter3')

    return column, colorColumn, color, axis, FieldBlue, FieldRed, MSmag, MSmagend, ClusterReddening, Startmag, Endmag

#opening the file and getting the data as a numpy array
def openingFile(filename):
    f = open(filename)
    csv_f = csv.reader(f)
    next(csv_f)
    photometry = np.array(list(csv_f),dtype=float)
    return photometry

#remove NAN values in an array and replace them with 0
def removeNAN(array):
    whereNANerror=np.isnan(array)
    array[whereNANerror]=0
    return array

#trying to clean up the data
def cleaningData(filter, plotting):
    print(filter)
    column, colorColumn, color, axis, FieldBlue, FieldRed, MSmag, MSmagend, ClusterReddening, Startmag, Endmag = settingParameters(filter)
    photometry = openingFile(filename)
    print(photometry.shape)
#getting rid of stars that are not measured or have error values (basically ridiculously high in magnitude or error)
    index=[]
    print(column)
    for x in range(0,len(photometry[:,column])):
        if photometry[x,column]>40:
            index.append(x)
        x+=1
    photometry=np.delete(photometry,index,axis=0)
    print(photometry.shape)
    index=[]
    for x in range(0,len(photometry[:,column+1])):
        if photometry[x,column+1]>10:
            index.append(x)
        x+=1

    photometry=np.delete(photometry,index,axis=0)
    print(photometry.shape)
#calculate the std and median for 20 magnitude bins of the photometry
#don't want to preferentially delete stars that are faint since those have higher errors
#only compare stars to those in their bins
    std,bin_edges,binnumber=scipy.stats.binned_statistic(photometry[:,column], photometry[:,column+1], statistic=np.std, bins=20, range=(np.amin(photometry[:,column]),np.amax(photometry[:,column])))
    median,bin_edges,binnumber=scipy.stats.binned_statistic(photometry[:,column], photometry[:,column+1], statistic='median', bins=20, range=(np.amin(photometry[:,column]),np.amax(photometry[:,column])))

    std = removeNAN(std)

    medianarray=np.zeros((1))
    goodindex=np.array([-1])

#finding the stars that are "good"
    for i in range(0,bin_edges.shape[0]-1):
        while True:
            if i==0:
                lowerlimit=0
                upperlimit=bin_edges[i]
            else:
                lowerlimit=bin_edges[i]
                upperlimit=bin_edges[i+1]
            sigmaOkay=medianandstd(column,lowerlimit,photometry,upperlimit,colorColumn,std,i,median,medianarray)
            arrayLength=len(medianarray)
            if np.isclose(medianarray[arrayLength-1],medianarray[arrayLength-2],atol=1e-02)==True:
                i=i+1
                goodindex=np.append(goodindex,sigmaOkay)
                break
    bin_edgesClipped=np.delete(bin_edges,0)
    splinearray=np.vstack((median+std*sigma,bin_edgesClipped-(bin_edges[1]-bin_edges[0])/2))
    splinearray = removeNAN(splinearray)
    fit=interpolate.InterpolatedUnivariateSpline(splinearray[1,:],splinearray[0,:], k=1)
    xnew = np.arange(19, 28, 0.1)
    ynew = fit(xnew)

#python hates me so I have to reopen this file for it to do this correctly
#this section is mapping out where on the chip the "good" and "bad" measurements are
#we want to make sure we're not being biased with respect to position in the cluster
    photometry = openingFile(filename)
    good=np.array([-1])
    bad=np.array([-1])
    
    for i in range(0,len(photometry[:,column])):
        if photometry[i,column+1]<fit(photometry[i,column]):
            good=np.append(good,i)
            i=1+i
        else:
            bad=np.append(bad,i   )
            i=1+i



    if plotting == True:
        #plot of stars' errors with the function used to delete stars
        plt.figure()
        plt.scatter(photometry[:,column],photometry[:,column+1],s=2)
        plt.plot(xnew,ynew)
        plt.axis((14,28,-0.05,0.25))
        plt.xlabel(str(filter),fontsize=18)
        plt.ylabel(str(filter)+' error',fontsize=18)
        plt.show()

        #graph bad stars with respect to their position on the chip
        plt.figure()
        plt.plot(photometry[bad,2],photometry[bad,3],'k.',markersize=2)
        plt.xlabel(r'x',fontsize=18)
        plt.ylabel(r'y',fontsize=18)
        plt.axis((0,5000,0,5000))
        plt.show()

        #Cleaned up CMD
        plt.figure(figsize=(8,16))
        plt.plot(photometry[:,column]-photometry[:,colorColumn],photometry[:,column],'k.',markersize=2,alpha=0.2)
        plt.xlabel(str(color),fontsize=18)
        plt.ylabel(str(filter),fontsize=18)
        plt.axis(axis)
        plt.gca().invert_yaxis()
        plt.show()

    return (goodindex,bin_edges,median,std,bin_edgesClipped,good,bad)

#this function is used in cleaningData
#it's goal is to get the allowed sigma average after we get rid of stars
#we iteratively get rid of stars until it's at the "allowed" sigma defined above with a tolerance of 1%
#I'm sorry about all the parameters here. Using a function within a function really sucks.

def medianandstd(column,lowerlimit,photometry,upperlimit,colorColumn,std,i,median,medianarray):
    magnitudeOkay=np.flatnonzero((lowerlimit<photometry[:,column]) & (photometry[:,column+1]<=upperlimit))
    
    if magnitudeOkay.shape==0:
        magnitudeOkay=np.array([0,1])
    
    subsetMagnitudeBin=photometry[magnitudeOkay,colorColumn]
    subsetErrorBin=photometry[magnitudeOkay,colorColumn+1]
    subsetErrorBin = removeNAN(subsetErrorBin)
    sigmaOkay=np.flatnonzero(subsetErrorBin<=sigma*std[i]+median[i])
    #subsetSubsetMagnitudeBin=subsetMagnitudeBin[sigmaOkay]
    #subsetSubsetErrorBin=subsetErrorBin[sigmaOkay]
    #median[i]=np.median(subsetSubsetErrorBin)
    #medianarray=np.append(medianarray,median[i])
    return sigmaOkay


#this function fits the MS ridgeline
def MSridgelineFitting(filter,finalIntersect, plotting):
    column, colorColumn, color, axis, FieldBlue, FieldRed, MSmag, MSmagend, ClusterReddening, Startmag, Endmag = settingParameters(filter)
    photometry = openingFile(filename)
    
    #this is adapted from Erin O'Malley's code MSfit1.py
    ClusterRed = ClusterReddening*1.281
    BlueEnd = FieldBlue+ClusterRed
    RedEnd = FieldRed+ClusterRed
    MSbins = ((MSmagend-MSmag)/0.1)-1
    MSmedColor = zeros(int(MSbins))
    MSmedMag = zeros(int(MSbins))
    MSmagErr = zeros(int(MSbins))
    MScolorErr = zeros(int(MSbins))
    ColorCluster=photometry[finalIntersect,column]-photometry[finalIntersect,colorColumn]
    ColorClusterErr=np.sqrt(np.power(photometry[finalIntersect,column+1],2)+np.power(photometry[finalIntersect,colorColumn+1],2))
    for b in range(len(MSmedColor)):
        condition2 = (photometry[finalIntersect,column]>MSmag)&(photometry[finalIntersect,column]<(MSmag+0.2))&(ColorCluster[:]>BlueEnd)&(ColorCluster[:]<RedEnd)
        ClusterMag = extract(condition2,photometry[finalIntersect,column])
        ClusterMagErr = extract(condition2,photometry[finalIntersect,column+1])
        ClusterColor = extract(condition2,ColorCluster[:])
        ClusterColorErr = extract(condition2,ColorClusterErr[:])
        ClusterColor = sorted(ClusterColor)
        MSmagErr[b] = sum(ClusterMagErr)/(len(ClusterMagErr)*1.0)
        MScolorErr[b] = sum(ClusterColorErr)/(len(ClusterColorErr)*1.0)
        length = len(ClusterColor)
        MSmedColor[b] = ClusterColor[int(length/2.0)]
        MSmedMag[b] = MSmag+0.1
        MSmag = MSmag + 0.1

    MSridgeline=np.polyfit(MSmedColor,MSmedMag, 1)
    MSridgelinex = np.arange(FieldBlue, FieldRed, 0.01)
    MSridgeliney=MSridgelinex*MSridgeline[0]+MSridgeline[1]

#want to find the nearest (on the chip) 50 stars to every MS star
#find if there is an offset in color for the star and correct for this offset
    MSstars=np.flatnonzero((MSmag<finalIntersectPhotometry[:,column]) & (finalIntersectPhotometry[:,column]<=MSmagend)&(FieldBlue<finalIntersectPhotometry[:,column]-finalIntersectPhotometry[:,colorColumn]) & (finalIntersectPhotometry[:,column]-finalIntersectPhotometry[:,colorColumn]<FieldRed))
    coords=np.column_stack((finalIntersectPhotometry[MSstars,xcoord],finalIntersectPhotometry[MSstars,ycoord]))
    distances=scipy.spatial.distance.cdist(coords,coords, 'euclidean')
    MSphotometry=finalIntersectPhotometry[MSstars,:]
    medianOffsetArray=np.zeros(len(distances))

    for i in range(0,len(distances)):
        sortedDistances=np.argsort(distances[i,:])
        nearest50=sortedDistances[:50]
        goodMagnitude=np.flatnonzero((MSphotometry[nearest50,column]<MSphotometry[i,column]+1)&(MSphotometry[nearest50,column]>MSphotometry[i,column]-1))

        #need distance from median ridgeline in color space
        actualColor=MSphotometry[nearest50[goodMagnitude],column]-MSphotometry[nearest50[goodMagnitude],colorColumn]
        theoreticalColor=(MSphotometry[nearest50[goodMagnitude],column]-MSridgeline[1])/MSridgeline[0]
        medianOffset=np.median(actualColor-theoreticalColor)
        medianOffsetArray[i]=medianOffset

        Color=finalIntersectPhotometry[:,column]-finalIntersectPhotometry[:,colorColumn]
            
        for a in range(0, len(MSstars)):
            index=MSstars[a]
            Color[index]=Color[index]-medianOffsetArray[a]


    if plotting == True:
        plt.figure(figsize=(8,16))
        plt.plot(photometry[:,column]-photometry[:,colorColumn],photometry[:,column],'k.',markersize=2)
        plt.plot(MSridgelinex,MSridgeliney)
        plt.xlabel(color,fontsize=18)
        plt.ylabel(str(filter),fontsize=18)
        plt.axis(axis)
        plt.gca().invert_yaxis()
        #plt.savefig(str(filter1)+'-'+str(filter2)+'Ridgeline.pdf')
        plt.show()

    return(Color, medianOffsetArray, MSstars)


#this function creates a better median ridgeline using rotated histograms
#then we'll straighten out the MS with respect to this line
#finally we can look at histograms of the spread of this function to see if there are two populations
def rotatedHistograms(filter,finalIntersect,Color,distanceError, plotting):
    column, colorColumn, color, axis, FieldBlue, FieldRed, MSmag, MSmagend, ClusterReddening, Startmag, Endmag = settingParameters(filter)
    photometry = openingFile(filename)

#get rid of the crazy outliers
    goodColorIndex=np.flatnonzero((Color>-10)&(Color<10))

#this rotated histogram part is adapted from Erin O'Malley's code MedRidge.py
    finalIntersectPhotometry=photometry[finalIntersect,:]
    VCluster = finalIntersectPhotometry[:,column2]
    ColorCluster = finalIntersectPhotometry[:,column]-finalIntersectPhotometry[:,colorColumn]
    
    Nbins = int(((Startmag-Endmag)/fitsize)-1)
    MedColor = np.zeros(Nbins)
    MedMag = np.zeros(Nbins)
    
    for b in range(len(MedColor)):
        # Set condition to make 0.5 mag bins that overlap at 0.05 mag
        condition2 = (VCluster[:]<Startmag) & (VCluster[:]>(Startmag-0.5)) & (ColorCluster[:]>FieldBlue) & (ColorCluster[:]<FieldRed)
        ClusterMag = extract(condition2,VCluster[:])
        ClusterColor = extract(condition2,ColorCluster[:])
        
        # Now we find the median in each bin by sorting the values and finding the middle
        ClusterColor = sorted(ClusterColor)
        MedColor[b] = ClusterColor[int(len(ClusterColor)/2.0)]
        MedMag[b] = Startmag-fitsize
        Startmag = Startmag-fitsize

    for n in range(3):
        for i in range(len(MedColor)-1):
            m = (MedMag[i]-MedMag[i+1])/(MedColor[i]-MedColor[i+1])
            t = np.arctan(1.0/m)
            
            # Find (x,y) coordinates for bin corners, find (x,y) coordinates for stars in bin, find median y, rotate back to (color,mag)
            xbin = np.cos(t)*MedColor[i]-np.sin(t)*MedMag[i]
            ybin = np.sin(t)*MedColor[i]+np.cos(t)*MedMag[i]
            clusterx = np.cos(t)*ColorCluster[:]-np.sin(t)*VCluster[:]
            clustery = np.sin(t)*ColorCluster[:]+np.cos(t)*VCluster[:]
            bincon = (clusterx[:]>xbin-0.25)&(clusterx[:]<xbin+0.25)&(clustery[:]>ybin-0.05)&(clustery[:]<ybin+0.05)
            x = extract(bincon,clusterx)
            xsort = sorted(x)
            Medx = xsort[int(len(xsort)/2.0)]
            
            # Rotate (Medx,0) back to (color,mag)
            MedColor[i] = np.cos(t)*Medx+np.sin(t)*ybin
            MedMag[i] = -np.sin(t)*Medx+np.cos(t)*ybin

    ColorFit=interp1d(MedMag, MedColor, fill_value='extrapolate')

    condition2 = (VCluster[:]<Startmag) & (VCluster[:]>(Endmag)) & (ColorCluster[:]>FieldBlue) & (ColorCluster[:]<FieldRed)
    a=np.flatnonzero(VCluster<Startmag)
    b=np.flatnonzero(VCluster>(Endmag))
    c=np.flatnonzero(ColorCluster>FieldBlue)
    d=np.flatnonzero(ColorCluster<FieldRed)
    e=np.intersect1d(a,b)
    f=np.intersect1d(c,d)
    selectedIndex=np.intersect1d(e,f)
    selectedStars=VCluster[selectedIndex]
    ClusterColor = ColorCluster[selectedIndex]
    distanceFromFit=np.zeros(ClusterColor.shape[0])
    allDistanceFromFit=np.zeros(ColorCluster.shape[0])
    
    for i in range(ClusterColor.shape[0]-1):
        fitValue=ColorFit(selectedStars[i])
        distanceFromFit[i]=ClusterColor[i]-fitValue

    for i in range(allDistanceFromFit.shape[0]):
        fitValue=ColorFit(VCluster[i])
        allDistanceFromFit[i]=ColorCluster[i]-fitValue

    xerror = np.arange(19, 30, 0.5)
    yzeros=np.zeros(xerror.shape[0])

    RedBlueFunction=interpolate.interp1d(xerror, distanceError,fill_value='extrapolate')
    xerrorRed = np.arange(MSmag2, MSmagend2, 0.5)

#get red and blue stars
    if filter==filter3:
        VCluster=finalIntersectPhotometry[:,6]
        redIndex=[]
        blueIndex=[]
        extremelyRedIndex=[]
        for i in range(len(xerrorRed)-1):
            lowerlimit=xerrorRed[i]
            upperlimit=xerrorRed[i+1]
            a=np.flatnonzero(VCluster>lowerlimit)
            b=np.flatnonzero(VCluster<upperlimit)
            magnitudeIndex=np.intersect1d(a,b)
            #RedBlueLimit=RedBlueConstant*distanceError[i]/2
            RedBlueLimit=RedBlueConstant*RedBlueFunction(VCluster)
            red=np.flatnonzero(allDistanceFromFit>RedBlueLimit)
            blue=np.flatnonzero(allDistanceFromFit<-RedBlueLimit)
            extremelyRed=np.flatnonzero(allDistanceFromFit>3*RedBlueLimit)
            
            #now keep ones in right magnitude bin
            redMagnitude=np.intersect1d(magnitudeIndex,red)
            blueMagnitude=np.intersect1d(magnitudeIndex,blue)
            extremelyRedMagnitude=np.intersect1d(magnitudeIndex,extremelyRed)
            redIndex.append(redMagnitude)
            blueIndex.append(blueMagnitude)
            extremelyRedIndex.append(extremelyRedMagnitude)

        redIndex=np.array(redIndex)
        blueIndex=np.array(blueIndex)
        extremelyRedIndex=np.array(extremelyRedIndex)
        redIndex=np.concatenate(redIndex).ravel()
        blueIndex=np.concatenate(blueIndex).ravel()
        extremelyRedIndex=np.concatenate(extremelyRedIndex).ravel()


        VCluster = finalIntersectPhotometry[selectedIndex,6]
        histogramBin=np.arange(10,30,binsize)
        for j in range(len(histogramBin)-1):
            lowerlimit=histogramBin[j]
            upperlimit=histogramBin[j+1]
            binIndex=np.flatnonzero((lowerlimit<VCluster) & (VCluster<=upperlimit))

#add in some binaries to the fit
        BinaryFit=interp1d(MedMag+0.75, MedColor, fill_value='extrapolate')



    if plotting == True:
        plt.figure(figsize=(8,16))
        plt.scatter(ColorCluster,VCluster,c='k',marker='.')
        plt.plot(MedColor,MedMag,'r')
        plt.axis(axis) #336-606
        plt.gca().invert_yaxis()
        plt.xlabel(color)
        plt.ylabel(str(filter))
        plt.show()

        plt.figure(figsize=(8,16))
        plt.scatter(ColorCluster,VCluster,c='k',marker='.')
        plt.plot(ColorFit(MedMag),MedMag,'r')
        plt.axis(axis) #336-606
        plt.gca().invert_yaxis()
        plt.xlabel(color)
        plt.ylabel(str(filter))
        plt.savefig('fit'+color+'.pdf')
        plt.show()

        plt.figure()
        plt.scatter(allDistanceFromFit,VCluster,c='k',marker='.')
        plt.scatter(allDistanceFromFit[redIndex],VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(allDistanceFromFit[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f'+color)
        plt.ylabel(str(filter))
        plt.show()

        plt.figure()
        plt.scatter(distanceFromFit,selectedStars,c='k',marker='.')
        plt.errorbar(yzeros, xerror, xerr=distanceError/2)
        #plt.axis(axis) #336-606
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f'+color)
        plt.ylabel(str(filter))
        #plt.savefig('distanceFromFit'+color+'.pdf')
        plt.show()

#dm is the horizontal shift from the two means
        if distanceFromFit[binIndex].shape[0] > 10:
            plt.figure()
            plt.hist(distanceFromFit[binIndex], bins='auto')
            title=str(lowerlimit)+str('-')+str(upperlimit)+str('mag ')+str(color)
            plt.gca().set_xlim(left=-0.6,right=0.6)
            #plt.axis((-0.5,0.5,0,300))
            plt.title(title)
            saveTitle=color+'.'+str(title)+str('.pdf')
            #plt.savefig(saveTitle)
            plt.show()

            hist=np.histogram(distanceFromFit[binIndex],bins='auto',normed=True)
            m, dm, sd1, sd2, weight = [0, 0.2, 0.2, 0.2, 0.10]
            p = [m, dm, sd1, sd2, weight] # Initial guesses for leastsq
            oldBins=hist[1]
            bins=np.zeros(len(oldBins)-1)
            for i in range(len(bins)):
                bins[i]=(oldBins[i]+oldBins[i+1])/2
    
            plsq,cov,infodict,mesg,ier = leastsq(res, p, args = (hist[0], bins),full_output = True)
            xrange=np.arange(-1.5, 1.5, 0.001)
            y_est1 = (1-plsq[4])*norm(xrange, plsq[0], plsq[2])
            y_est2 = plsq[4]*norm(xrange, plsq[0] + plsq[1], plsq[3])
            y_est3 = y_est1+y_est2
            
            chi2=(infodict['fvec']**2).sum()/(len(infodict['fvec'])-len(plsq))

            std=np.std(distanceFromFit[binIndex])
    
            leftStars=np.flatnonzero(distanceFromFit[binIndex]<0).shape[0]
            rightStars=np.flatnonzero(distanceFromFit[binIndex]>0).shape[0]
            midStars=np.flatnonzero(distanceFromFit[binIndex]==0).shape[0]

            plt.figure()
            plt.hist(distanceFromFit[binIndex],bins='auto',normed=True)
            #plt.title(title+str(' weight:')+str(round(plsq[4],3))+str(' chi2: ')+str(round(chi2, 3)))
            plt.title(title+str(' weight:')+str(round(plsq[4],3))+str(round(std,5)))
            plt.xlabel('Distance From Fit')
            plt.ylabel('Normalized Counts')
            plt.plot(xrange,y_est1,'r',linewidth=2)
            plt.plot(xrange,y_est2,'y',linewidth=2)
            plt.plot(xrange,y_est3,'k',linewidth=2)
            plt.show()

            #single Gaussian fit
            m, sd1 = [0, 0.05]
            p = [m, sd1]
            plsqsingle=leastsq(ressingle, p, args = (hist[0], bins))
            y_est=norm(bins, plsq[0], plsq[1])
            chi2=np.sum(np.square(hist[0]-y_est)/y_est)
            
            plt.figure()
            plt.hist(distanceFromFit[binIndex],bins='auto',normed=True)
            plt.title(title+str(' chi2: ')+str(chi2))
            plt.plot(bins,y_est,'k',linewidth=2)
            plt.show()
            
            numberPoints=binIndex.shape[0]
            #make a grid of models
            mModel=np.arange(-0.1,0.1,0.025)
            dmModel=np.arange(-0.1,0.5,0.025)
            sd1Model=np.arange(0.025,0.5,0.025)
            sd2Model=np.arange(0.025,0.5,0.025)
            weightModel=np.arange(0,0.5,0.025)

    if filter==filter3:
        return allDistanceFromFit,redIndex,blueIndex,extremelyRedIndex,selectedIndex
    else:
        return allDistanceFromFit,selectedIndex

def sequenceofScripts(filter, plottingvalue):
    print('sequenceofScripts')
    print(filter)
    Data=cleaningData(filter, plottingvalue)
    goodindex=Data[0]
    bin_edges=Data[1]
    median=Data[2]
    std=Data[3]
    bin_edgesClipped=Data[4]
    good=Data[5]
    bad=Data[6]
    print('gothereabc')
#for getting error bars on the straightened plots
    averageerror=np.vstack((median+std,bin_edgesClipped-(bin_edges[1]-bin_edges[0])/2))
    averageerror=np.nan_to_num(averageerror)
    ferror=interpolate.InterpolatedUnivariateSpline(averageerror[1,:],averageerror[0,:], k=1)
    yerror = ferror(xerror)

    return Data, goodindex, bin_edges, median, std, bin_edgesClipped, good, bad, averageerror, ferror, yerror


if __name__ == '__main__':
    #arguments to give which column has what data
    #    script, filename, filter1, column1, filter2, column2, filter3, column3, xcoord, ycoord = argv
    filename, filter1, column1, filter2, column2, filter3, column3, xcoord, ycoord = 'NGC1466.csv', 336, 4, 606, 6, 814, 8, 2, 3
    #convertig inputs to integers
    print(filter1)
    filter1=int(filter1)
    column1=int(column1)
    filter2=int(filter2)
    column2=int(column2)
    filter3=int(filter3)
    column3=int(column3)
    xcoord=int(xcoord)
    ycoord=int(ycoord)

    #color named for plots and saving files
    color1=str(filter1)+'-'+str(filter2)
    color2=str(filter2)+'-'+str(filter3)
    color3=str(filter1)+'-'+str(filter3)

    sigma=100000

    #these values are found for each cluster with trial and error
    RedEdge1=1.75
    RedEdge2=1.0
    RedEdge3=3.0

    BlueEdge1=0.5
    BlueEdge2=0.5
    BlueEdge3=1.0

    MSmag1=24
    MSmag2=24
    MSmag3=24.25

    MSmagend1=27
    MSmagend2=26
    MSmagend3=27.5

    axis1=(0,2.5,19,28)
    axis2=(-0,1.5,16,28)
    axis3=(0,4,19,30)

    ClusterReddening1=0
    ClusterReddening2=0
    ClusterReddening3=0
        


    binaryFraction=0.1

    RedBlueConstant=1

    xerror = np.arange(19, 30, 0.5)

    photometry = openingFile(filename)


    #I keep reading that you should never use global variables, but I don't care. Only way I could get it to work.
    global medianarray
    print('gothere')
    #running the scripts
    print(filter1)
    Data336, goodindex336, bin_edges336, median336, std336, bin_edgesClipped336, good336, bad336, averageerror336, ferror336, yerror336 = sequenceofScripts(filter1, 'FALSE')

    Data606, goodindex606, bin_edges606, median606, std606, bin_edgesClipped606, good606, bad606, averageerror606, ferror606, yerror606 = sequenceofScripts(filter2, 'FALSE')


    Data814, goodindex814, bin_edges814, median814, std814, bin_edgesClipped814, good814, bad814, averageerror814, ferror814, yerror814 = sequenceofScripts(filter3, 'FALSE')
    print('gothere2')
    colorerror336606=np.sqrt(np.square(yerror336)+np.square(yerror606))
    colorerror606814=np.sqrt(np.square(yerror814)+np.square(yerror606))
    colorerror336814=np.sqrt(np.square(yerror336)+np.square(yerror814))

    #stars that are "good" in all three filters
    intersect=np.intersect1d(good336,good606)
    finalIntersect=np.intersect1d(intersect,good814)

    #need to create a radial cut of stars; want them only a certain distance from the center of the cluster

    #finding the center of the cluster

    xCenter=np.median(photometry[finalIntersect,2])
    yCenter=np.median(photometry[finalIntersect,3])


    #want the cluster size to be where the density of stars stays constant

    ClusterSize=100
    density1=0
    density2=1000

    while (np.isclose(density1,density2,atol=2e-05)==False):
        distanceFromCenter=np.sqrt(np.power(xCenter-photometry[finalIntersect,2],2)+np.power(yCenter-photometry[finalIntersect,3],2))
        numberStars1=np.flatnonzero((distanceFromCenter<=ClusterSize)&(distanceFromCenter>ClusterSize-50)).shape[0]
        area1=np.pi*(np.power(ClusterSize,2)-np.power(ClusterSize-50,2))
        density1=numberStars1/area1
        
        numberStars2=np.flatnonzero((distanceFromCenter<=ClusterSize+50)&(distanceFromCenter>ClusterSize)).shape[0]
        area2=np.pi*(np.power(ClusterSize+50,2)-np.power(ClusterSize,2))
        density2=numberStars2/area2
        
        ClusterSize=ClusterSize+50

    ClusterSize=3000 #basically use this if you don't really want to fool around with the cluster size
    print('gothere3')
    distanceFromCenter=np.sqrt(np.power(xCenter-photometry[finalIntersect,2],2)+np.power(yCenter-photometry[finalIntersect,3],2))
    clusterStars=np.flatnonzero(distanceFromCenter<=ClusterSize)

    finalIntersect=finalIntersect[clusterStars]
    print('gothere4')
    fig=plt.figure()
    plt.title('Cluster Size')
    plt.scatter(photometry[:,2], photometry[:,3],c='k',alpha=0.1)
    ax=fig.add_subplot(1,1,1)
    circ=plt.Circle((xCenter,yCenter), radius=ClusterSize, color='g', fill=False)
    ax.add_patch(circ)
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.show()
    
    fig=plt.figure()
    plt.title('Final Selected Stars')
    plt.scatter(photometry[finalIntersect,2], photometry[finalIntersect,3],c='k',alpha=0.1)
    ax=fig.add_subplot(1,1,1)
    circ=plt.Circle((xCenter,yCenter), radius=ClusterSize, color='g', fill=False)
    ax.add_patch(circ)
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.show()
    
    limits=np.arange(19.25, 30.25, 0.5)

    magvalue=np.zeros(limits.shape[0])
    yzeros=np.zeros(magvalue.shape[0])
    for i in range(0, limits.shape[0]-1):
        magvalue[i]=(limits[i]+limits[i+1])/2

    error336=np.zeros(limits.shape[0])
    error336vs606=np.zeros(limits.shape[0])
    error336vs814=np.zeros(limits.shape[0])

    error606=np.zeros(limits.shape[0])
    error606vs336=np.zeros(limits.shape[0])
    error606vs606=np.zeros(limits.shape[0])
    error606vs814=np.zeros(limits.shape[0])

    error814=np.zeros(limits.shape[0])
    error814vs336=np.zeros(limits.shape[0])
    error814vs606=np.zeros(limits.shape[0])

    finalIntersectPhotometry=photometry[finalIntersect,:]
    print('I made it here!')
    for i in range(0, limits.shape[0]-1):
        good814=np.flatnonzero(1>finalIntersectPhotometry[:,9])
        lowerlimit=limits[i]
        upperlimit=limits[i+1]
        starIndex336=np.flatnonzero((lowerlimit<finalIntersectPhotometry[:,4]) & (finalIntersectPhotometry[:,4]<=upperlimit))
        starIndex336=np.intersect1d(starIndex336,good814)
        error336[i]=np.nanmean(finalIntersectPhotometry[starIndex336,5])
        error336vs606[i]=np.nanmean(finalIntersectPhotometry[starIndex336,7])
        error336vs814[i]=np.nanmean(finalIntersectPhotometry[starIndex336,9])
        
        starIndex606=np.flatnonzero((lowerlimit<finalIntersectPhotometry[:,6]) & (finalIntersectPhotometry[:,6]<=upperlimit))
        starIndex606=np.intersect1d(starIndex606,good814)
        error606[i]=np.nanmean(finalIntersectPhotometry[starIndex606,7])
        error606vs336[i]=np.nanmean(finalIntersectPhotometry[starIndex606,5])
        error606vs606[i]=np.nanmean(finalIntersectPhotometry[starIndex606,7])
        error606vs814[i]=np.nanmean(finalIntersectPhotometry[starIndex606,9])
        
        starIndex814=np.flatnonzero((lowerlimit<finalIntersectPhotometry[:,8]) & (finalIntersectPhotometry[:,8]<=upperlimit))
        starIndex814=np.intersect1d(starIndex814,good814)
        error814[i]=np.nanmean(finalIntersectPhotometry[starIndex814,9])
        error814vs336[i]=np.nanmean(finalIntersectPhotometry[starIndex814,5])
        error814vs606[i]=np.nanmean(finalIntersectPhotometry[starIndex814,7])

    distanceError336606814=np.nan_to_num(np.sqrt(np.square(error336vs606)+np.square(error336vs814)))
    distanceError606336814=np.nan_to_num(np.sqrt(np.square(error606vs336)+np.square(error606vs814)))
    distanceError814336606=np.nan_to_num(np.sqrt(np.square(error814vs606)+np.square(error814vs336)))

    distanceError606606814=np.nan_to_num(np.sqrt(np.square(error606vs606)+np.square(error606vs814)))
    distanceError606336606=np.nan_to_num(np.sqrt(np.square(error606vs606)+np.square(error606vs336)))


    #Color Magnitude Diagrams (CMDs)
    print('Whew, made it to the graphs')
    #@@@PLOT@@@
    plt.figure(figsize=(8,16))
    plt.plot(photometry[finalIntersect,column1]-photometry[finalIntersect,column2],photometry[finalIntersect,column1],'k.',markersize=2)
    plt.xlabel(color1,fontsize=18)
    plt.ylabel(filter1,fontsize=18)
    plt.axis(axis1)
    plt.gca().invert_yaxis()
    #plt.savefig(color1+'WellMeasured.pdf')
    #plt.show()

    #@@@PLOT@@@
    plt.figure(figsize=(8,16))
    plt.plot(photometry[finalIntersect,column2]-photometry[finalIntersect,column3],photometry[finalIntersect,column2],'k.',markersize=2)
    plt.xlabel(color2,fontsize=18)
    plt.ylabel(filter2,fontsize=18)
    plt.axis(axis2)
    plt.gca().invert_yaxis()
    #plt.savefig(color2+'WellMeasured.pdf')
    #plt.show()

    #@@@PLOT@@@
    plt.figure(figsize=(8,16))
    fig_size = plt.rcParams["figure.figsize"]
    #fig_size[0] = 6
    #fig_size[1] = 12
    plt.plot(photometry[finalIntersect,column1]-photometry[finalIntersect,column3],photometry[finalIntersect,column1],'k.',markersize=2)
    plt.xlabel(color3,fontsize=18)
    plt.ylabel(filter1,fontsize=18)
    plt.axis(axis3)
    plt.gca().invert_yaxis()
    #plt.savefig(color3+'WellMeasured.pdf')
    #plt.show()


    print('Fitting yay!')
    Color336606,medianOffset336606, MS336 = MSridgelineFitting(filter1,finalIntersect, 'FALSE')
    Color606814,medianOffset606814, MS606 = MSridgelineFitting(filter2,finalIntersect, 'FALSE')
    Color336814,medianOffset336814, MS814 = MSridgelineFitting(filter3,finalIntersect, 'FALSE')


    #fitsize is the size of your magnitude bins for your rotated histogram median ridgeline fit
    #the smaller, the more accurate but it gets much less smooth
    #binsize is the size of the magnitude bin for my histograms to see if there are multiple populations
    #it's probably good to run each cluster with a variety of values for each of these parameters

    fitsize=0.3
    binsize=0.5

    distanceFit336606,selectedIndex336=rotatedHistograms(filter1,finalIntersect,Color336606,distanceError814336606, 'FALSE')
    distanceFit606814,selectedIndex606=rotatedHistograms(filter2,finalIntersect,Color606814,distanceError336606814, 'FALSE')
    distanceFit336814,redIndex,blueIndex,extremelyRedIndex,selectedIndex814=rotatedHistograms(filter3,finalIntersect,Color336814,distanceError606336814, 'FALSE')





    finalIntersectPhotometry=photometry[finalIntersect,:]
    VCluster = finalIntersectPhotometry[:,6]


    blackIndex=np.arange(VCluster.shape[0])
    redBlueIndex=np.concatenate((redIndex,blueIndex))
    blackIndex=np.delete(blackIndex,redBlueIndex)

    #do the histograms again without the red and blue stars included

    noRedBlue=np.delete(finalIntersectPhotometry, (redBlueIndex),axis=0)

    plottingQuestion = 'FALSE'
    print('More plots')
    if plottingQuestion == 'True':

        plt.figure()
        plt.scatter(distanceFit336606,VCluster,c='k',marker='.')
        plt.scatter(distanceFit336606[redIndex],VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit336606[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-606)')
        plt.ylabel('606')
        plt.show()

        plt.figure()
        plt.scatter(distanceFit336606,VCluster,c='k',marker='.')
        plt.scatter(distanceFit336606[redIndex],VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit336606[extremelyRedIndex],VCluster[extremelyRedIndex],c='g',marker='.',edgecolor='face')
        plt.scatter(distanceFit336606[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-606)')
        plt.ylabel('606')
        plt.show()

        plt.figure()
        plt.scatter(distanceFit606814,VCluster,c='k',marker='.')
        plt.scatter(distanceFit606814[redIndex]-0.02,VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit606814[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(606-814)')
        plt.ylabel('606')
        plt.show()

        plt.figure()
        plt.scatter(distanceFit606814,VCluster,c='k',marker='.')
        plt.scatter(distanceFit606814[redIndex]-0.02,VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit606814[extremelyRedIndex]-0.02,VCluster[extremelyRedIndex],c='g',marker='.',edgecolor='face')
        plt.scatter(distanceFit606814[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(606-814)')
        plt.ylabel('606')
        plt.show()

        plt.figure()
        plt.scatter(distanceFit336814,VCluster,c='k',marker='.')
        plt.scatter(distanceFit336814[redIndex],VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit336814[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-814)')
        plt.ylabel('606')
        plt.show()

        plt.figure()
        plt.scatter(distanceFit336814,VCluster,c='k',marker='.')
        plt.scatter(distanceFit336814[redIndex],VCluster[redIndex],c='r',marker='.',edgecolor='face')
        plt.scatter(distanceFit336814[extremelyRedIndex],VCluster[extremelyRedIndex],c='g',marker='.',edgecolor='face')
        plt.scatter(distanceFit336814[blueIndex],VCluster[blueIndex],c='b',marker='.',edgecolor='face')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-814)')
        plt.ylabel('606')
        plt.show()


        plt.figure()
        plt.scatter(finalIntersectPhotometry[redIndex,2],finalIntersectPhotometry[redIndex,3],c='r',alpha=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[blueIndex,2],finalIntersectPhotometry[blueIndex,3],c='b',alpha=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[blackIndex,2],finalIntersectPhotometry[blackIndex,3],c='k',alpha=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[extremelyRedIndex,2],finalIntersectPhotometry[extremelyRedIndex,3],c='g',alpha=0.1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[MS336,xcoord], finalIntersectPhotometry[MS336,ycoord],c=medianOffset336606)
        plt.title('F336W-F606W')
        plt.xlabel('x')
        plt.ylabel('y')
        cbar = plt.colorbar()
        cbar.set_label('Color Correction',rotation=270,labelpad=10)
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[MS606,xcoord], finalIntersectPhotometry[MS606,ycoord],c=medianOffset606814)
        plt.title('F606W-F814W')
        plt.xlabel('x')
        plt.ylabel('y')
        cbar = plt.colorbar()
        cbar.set_label('Color Correction',rotation=270,labelpad=10)
        plt.show()

        plt.figure()
        plt.scatter(finalIntersectPhotometry[MS814,xcoord], finalIntersectPhotometry[MS814,ycoord],c=medianOffset336814)
        plt.title('F336W-F814W')
        plt.xlabel('x')
        plt.ylabel('y')
        cbar = plt.colorbar()
        cbar.set_label('Color Correction',rotation=270,labelpad=10)
        plt.show()

        plt.figure()
        plt.scatter(distanceFit336606,finalIntersectPhotometry[:,8],c='k',marker='.')
        #plt.scatter(distanceFit336606[redIndex],finalIntersectPhotometry[redIndex,8],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit336606[blueIndex],finalIntersectPhotometry[blueIndex,8],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError814336606/2, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-606)')
        plt.ylabel('814')
        plt.show()


        plt.figure()
        plt.scatter(distanceFit336814,finalIntersectPhotometry[:,6],c='k',marker='.')
        #plt.scatter(distanceFit336814[redIndex],finalIntersectPhotometry[redIndex,6],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit336814[blueIndex],finalIntersectPhotometry[blueIndex,6],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError606336814/2, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-814)')
        plt.ylabel('606')
        plt.show()


        plt.figure()
        plt.scatter(distanceFit606814,finalIntersectPhotometry[:,4],c='k',marker='.')
        #plt.scatter(distanceFit606814[redIndex],finalIntersectPhotometry[redIndex,4],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit606814[blueIndex],finalIntersectPhotometry[blueIndex,4],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError336606814/2, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(606-814)')
        plt.ylabel('336')
        plt.show()


        plt.figure()
        plt.scatter(distanceFit336606,finalIntersectPhotometry[:,6],c='k',marker='.')
        #plt.scatter(distanceFit336606[redIndex],finalIntersectPhotometry[redIndex,6],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit336606[blueIndex],finalIntersectPhotometry[blueIndex,6],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError606606814/1.7, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-606)')
        plt.ylabel('606')
        plt.show()


        plt.figure()
        plt.scatter(distanceFit336814,finalIntersectPhotometry[:,6],c='k',marker='.')
        #plt.scatter(distanceFit336814[redIndex],finalIntersectPhotometry[redIndex,6],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit336814[blueIndex],finalIntersectPhotometry[blueIndex,6],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError606606814/1.8, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(336-814)')
        plt.ylabel('606')
        plt.show()


        plt.figure()
        plt.scatter(distanceFit606814,finalIntersectPhotometry[:,6],c='k',marker='.')
        #plt.scatter(distanceFit606814[redIndex],finalIntersectPhotometry[redIndex,6],c='r',marker='.',edgecolor='face')
        #plt.scatter(distanceFit606814[blueIndex],finalIntersectPhotometry[blueIndex,6],c='b',marker='.',edgecolor='face')
        plt.errorbar(yzeros, magvalue, xerr=distanceError606606814/2, elinewidth=2, ecolor='cyan', color='cyan')
        plt.axis((-0.5,0.5,24.5,28))
        plt.gca().invert_yaxis()
        plt.xlabel('f(606-814)')
        plt.ylabel('606')
        plt.show()












