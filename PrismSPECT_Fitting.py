## PrismSPECT_Fitting.py will fit the lines that you want from the PrismSPECT 
## text files.  Part of HADHES version 0.2.0.
## T.E.Steinberger July 2016

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from scipy.special import wofz
import Tkinter as tk
import tkMessageBox as mb
import tkFileDialog as fd


## The initial function creates crucial arrays and values that will be used
## in the rest of the suite of codes.  Temperature and Density matrix are
## just arrays of the temperatures and densities that you used. Number_of_
## " : are self explanatory parameters. Color_Max gives a parameter for the
## colors array to run till.  This is what changes the colors in for loops
## when we have a bunch of files.  area, and width are empty arrays made 
## that will be filled in later for efficiency.
def Initial():
    Temperature_Matrix = np.array([30, 32.5, 35,37.5, 40,42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70, 72.5, 75, 77.5, 80, 82.5, 85, 87.5, 90, 92.5, 95])
    Density_Matrix = np.array([1e20, 1.25e20, 1.5e20, 1.75e20, 2e20, 2.25e20, 2.5e20, 2.75e20, 3e20])
    Number_of_Temperatures = len(Temperature_Matrix)
    Number_of_Densities = len(Density_Matrix)
    Number_of_Files = Number_of_Temperatures*Number_of_Densities
    Color_Max = max(Number_of_Densities, Number_of_Temperatures)
    cmap = plt.get_cmap('gnuplot')
    colors = [cmap(i) for i in np.linspace(0, 1, Color_Max)]
    width1 = np.zeros(Number_of_Files)
    width2 = np.zeros(Number_of_Files)
    area1 = np.zeros(Number_of_Files)
    area2 = np.zeros(Number_of_Files)
    Fits1 = np.zeros((Number_of_Files, 8))
    Fits2 = np.zeros((Number_of_Files, 8))
    Real_Wavelength = np.genfromtxt(Directory + '/Averaged_Data.txt', usecols = (0))
    Real_Transmission = np.genfromtxt(Directory + '/Averaged_Data.txt', usecols = (1))
    Window_Length = np.abs(Real_Wavelength[0] - Real_Wavelength[len(Real_Wavelength)-1])
    return Number_of_Temperatures, Temperature_Matrix, Number_of_Densities, Density_Matrix, Number_of_Files, colors, area1, area2, width1, width2, Real_Wavelength, Real_Transmission, Window_Length, Fits1, Fits2

## Used for opening old files.
global File_Location

## This function is what allows you to grab points from a plot.  Returns 
## only the x values (xcoords) but one could modify this to return the 
## y values as well.  This is used to restrict the Transmission arrays
def Get_Points(x, y, x1, y1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.plot(x1, y1)
    xcoords = []
    ycoords = []
    def onclick(event):
        ix, iy = (float(event.xdata), float(event.ydata))
        xcoords.append(ix)
        ycoords.append(iy)
        return xcoords, ycoords
    cid = fig.canvas.mpl_connect('key_press_event', onclick)
    plt.title('Select ranges for lines of interest.')
    plt.show()
    return xcoords[0], xcoords[1], xcoords[2], xcoords[3]

## This function returns the parameters for the deconvolution of three
## peaks. Will return an array that will be used for the scipy fitting.
## [intensity1, position1, width1, ...] is the format (goes to three). 
def Get_Params_Three(x, y, x1, y1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.plot(x1, y1)
    xcoords = []
    ycoords = []
    def onclick(event):
        ix, iy = (float(event.xdata), float(event.ydata))
        xcoords.append(ix)
        ycoords.append(iy)
        return xcoords, ycoords
    cid = fig.canvas.mpl_connect('key_press_event', onclick)
    plt.title('Select the parameter points')
    plt.show()
    Width_Array1, Width_Array2, Width_Array3 = (xcoords[2] - xcoords[0]), np.arange(xcoords[3], xcoords[4], np.abs(xcoords[3]-xcoords[4])/10), (xcoords[10] - xcoords[8])
    Intensity_Array1, Intensity_Array2, Intensity_Array3 = (1 - ycoords[1]), np.arange((1-ycoords[5]), (1-ycoords[6]),  np.abs((1-ycoords[5])-(ycoords[6]))/5), (1 - ycoords[9])
    Postion1, Position2, Position3 = xcoords[1], xcoords[7], xcoords[9]
    return Width_Array1, Width_Array2, Width_Array3, Intensity_Array1, Intensity_Array2, Intensity_Array3, Postion1, Position2, Position3

## This is the same as the one above.  The only difference is that this 
## function will return the parameters for two lines instead of three.
def Get_Params_Two(x, y, x1, y1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.plot(x1, y1)
    xcoords = []
    ycoords = []
    def onclick(event):
        ix, iy = (float(event.xdata), float(event.ydata))
        xcoords.append(ix)
        ycoords.append(iy)
        return xcoords, ycoords
    cid = fig.canvas.mpl_connect('key_press_event', onclick)
    plt.title('Select the parameter points')
    plt.show()
    Width_Array1, Width_Array2 = np.arange(xcoords[0], xcoords[1], np.abs(xcoords[0]-xcoords[1])/10), np.arange(xcoords[5], xcoords[6], np.abs(xcoords[5]-xcoords[6])/10)
    Intensity_Array1, Intensity_Array2 = np.arange((1-ycoords[2]), (1-ycoords[3]), np.abs(((1-ycoords[3])-(1-ycoords[2])))/5), np.arange((1-ycoords[7]),(1- ycoords[8]), np.abs(((1-ycoords[7])-(1-ycoords[8])))/5)
    Position1, Position2 = xcoords[4], xcoords[9]
    return Width_Array1, Width_Array2, Intensity_Array1, Intensity_Array2, Position1, Position2

## This is the same for Get_Params_Three and Two.  Only difference is that## it returns parameters for one line instead of two or three.
def Get_Params_One(x, y, x1, y1):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y)
    ax.plot(x1, y1)
    xcoords = []
    ycoords = []
    def onclick(event):
        ix, iy = (float(event.xdata), float(event.ydata))
        xcoords.append(ix)
        ycoords.append(iy)
        return xcoords, ycoords
    cid = fig.canvas.mpl_connect('key_press_event', onclick)
    plt.title('Select the parameter points')
    plt.show()
    Width_Array = np.arange(xcoords[0], xcoords[1], np.abs(xcoords[0]-xcoords[1])/10)
    Intensity_Array = np.arange(1-(ycoords[2]), (1-ycoords[3]), np.abs((1-(ycoords[3])-(1-ycoords[2])))/5)
    Position = xcoords[4]
    return Width_Array, Intensity_Array, Position

## This is the function that defines what a voigt profile is.  This one
## will be used for the two or three peaks.
def voigt(x, pars):
    a = pars[0] #Intensity
    b = pars[1] #Mean Center
    c = pars[2] #Width of Lorentzian
    z = ((x-b) + 1j*complex(c))/(np.sqrt(2)*c)
    v = -np.real(wofz(z))/(c*np.sqrt(2*np.pi)) #Voigt Profile
    V = (v/((abs(min(v)))/a))+1 #Normailized to be between 0 and 1
    return V #Lorentzian

## This is the function that will return the fitted voigt profile for
## this is used for fitting just one line not two or three.
def voigt1(x, *pars):
    a = pars[0] #Intensity
    b = pars[1] #Mean Center
    c = pars[2] #Width of Lorentzian
    z = ((x-b) + 1j*complex(c))/(np.sqrt(2)*c)
    v = -np.real(wofz(z))/(c*np.sqrt(2*np.pi)) #Voigt Profile
    V = (v/((abs(min(v)))/a))+1 #Normailized to be between 0 and 1
    return V #Lorentzian

## Takes the parameters that are grabbed from the plot and the voigt 
## function and fits three peaks of which the center peak will be 
## returned.
def three_peaks(x, *pars):
    a0 = pars[0] #Intensity of first peak
    b0 = pars[1] #Mean center of first peak
    c0 = pars[2] #Width of first peak
    a1 = pars[3] #Intensity of second peak
    b1 = pars[4] #Mean center of second peak
    c1 = pars[5] #Width of second peak
    a2 = pars[6] #Intensity of third peak
    b2 = pars[7] #Mean center of third peak
    c2 = pars[8] #Width of third peak
    peak1 = voigt(x, [a0, b0, c0]) #First peak
    peak2 = voigt(x, [a1, b1, c1]) #Second peak
    peak3 = voigt(x, [a2, b2, c2]) #third peak
    return (peak1 + peak2 + peak3)-2 #Convolved peaks

## Same as the above except for two peaks instead of three peaks.
def two_peaks(x, *pars):
    a0 = pars[0] #Intensity of first peak
    b0 = pars[1] #Mean center of first peak
    c0 = pars[2] #Width of first peak
    a1 = pars[3] #Intensity of second peak
    b1 = pars[4] #Mean center of second peak
    c1 = pars[5] #Width of second peak
    peak1 = voigt(x, [a0, b0, c0]) #First peak
    peak2 = voigt(x, [a1, b1, c1]) #Second peak
    return (peak1 + peak2)-1 #Convolved peaks

## This will take the fit from the three peaks and calculate the area to
## the window that the real data has.
def Area_Three(x, y, parameters, Window_Length, Real_Wavelength):
    a = parameters[0:3] # Parameters for first line
    b = parameters[3:6] # Parameters for second lines
    c = parameters[6:9] # Parameters for third line
    peak1 = voigt(x, a) # First fit
    peak2 = voigt(x, b) # Second fit
    peak3 = voigt(x, c) # Third fit
    area1 = (Window_Length) - abs(np.trapz(peak1[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])], x[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])]))
    return area1

## Same as Area_Three except for two peaks (returns both areas based on
## the ones that you pick.
def Area_Two(x, y, parameters, Window_Length, Real_Wavelength):
    a = parameters[0:3]
    b = parameters[3:6]
    peak1 = voigt(x, a)
    peak2 = voigt(x, b)
    area1 = (Window_Length) - abs(np.trapz(peak1[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])], x[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])]))
    area2 = (Window_Length) - abs(np.trapz(peak2[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])], x[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])]))
    return [area1, area2]

## Same as Area_Three and Two but will return area for on line instead of
## two or three.
def Area_One(x, y, Window_Length, Real_Wavelength, *parameters):
    peak1 = voigt1(x, *parameters)
    area1 = (Window_Length) - abs(np.trapz(peak1[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])], x[(x<=Real_Wavelength[len(Real_Wavelength)-1]) & (x>=Real_Wavelength[0])]))
    return area1

## Plots the most extreme PrismSPECT file and the least extreme PrismSPECT 
## file so that you can accurately choose line parameters from two extremes
def Extreme_Plot(n):
    Data_least = np.loadtxt(File_Location + '/PRISMSPECT0.ppd')
    Data_most = np.loadtxt(File_Location + '/PRISMSPECT' + str(n-1) + '.ppd')
    xl = Data_least[:, 0]
    yl = Data_least[:, 1]
    xm = Data_most[:, 0]
    ym = Data_most[:, 1]
    return xl, yl, xm, ym

## Ask where to draw from.
root = tk.Tk()  # Has to do with making a window NOT pop up
root.withdraw()
Directory = fd.askdirectory(title = 'Directory for saved averaged data.') # ask where to save
root.destroy()

## This is the start of fitting the PrismSPECT data.
Number_of_Temperatures, Temperature_Matrix, Number_of_Densities, Density_Matrix, Number_of_Files, colors, area1, area2, width1, width2, Real_Wavelength, Real_Transmission, Window_Length, Fits1, Fits2 = Initial()
u = 0 # Setting a counter variable.
Temperatures = np.zeros(Number_of_Files).reshape(Number_of_Files, 1) # for text files
Densities = np.zeros(Number_of_Files).reshape(Number_of_Files, 1) # for text files. 

## This while loop creats temperature and density arrays that will be
## used in the text file that is saved.
while u/Number_of_Densities < Number_of_Temperatures:
    for i in range(0 + u, Number_of_Densities+u):
        Temperatures[i] = Temperature_Matrix[u/Number_of_Densities]
        Densities[i] = Density_Matrix[i-u]
    u = u + Number_of_Densities

root = tk.Tk()  # Has to do with making a window NOT pop up
root.withdraw()
File_Location = fd.askdirectory(title = 'Directory containing PRISMSpect files') # ask where to save
root.destroy()

## Get the number of lines to be deconvoluted and the ranges.
Peaks1, Peaks2 = input('How many peaks are to be deconvoluted for the first and second lines (separate by comma)? 1, 2, or 3?\n')
left1, right1, left2, right2 = Get_Points(*Extreme_Plot(Number_of_Files))

## Start of fitting for the first peak.
if Peaks1 == 1:
    Area_Params = ['# Area Parameters for first line']
    Width_Array, Intensity_Array, Position = Get_Params_One(*Extreme_Plot(Number_of_Files))
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left1)|(Wavelength>right1)] = 1
        tick = 0
        Fit_Parameters = np.zeros((((len(Width_Array)/2)*len(Intensity_Array)), 3))
        Area_1 = np.zeros(((len(Width_Array)/2)*len(Intensity_Array)))
        for intensity in range(0, len(Intensity_Array)):
            for width in range(0, len(Width_Array)/2):
                popt, pcov = so.curve_fit(voigt1, Wavelength, Transmission, [Intensity_Array[intensity], Position, (np.abs(Width_Array[width] - Width_Array[(len(Width_Array)-1)-width]))], maxfev = 100000)
                if any(popt < 0):
                    popt = [Intensity_Array[intensity], Position, (np.abs(Width_Array[width] - Width_Array[(len(Width_Array)-1)-width]))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick] = Area_One(Wavelength, Transmission, Window_Length, Real_Wavelength, *popt)
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array)/2)*len(Intensity_Array))))*100) + "% for complete for file " + str(i+1)
        print "File " + str(i+1) + " of " + str(Number_of_Files) + " complete."
        Average_Intensity, StDev_Intensity = np.average(Fit_Parameters[:, 0]), np.std(Fit_Parameters[:, 0])
        Average_Position, StDev_Position = np.average(Fit_Parameters[:, 1]), np.std(Fit_Parameters[:, 1])
        Average_Width, StDev_Width = np.average(Fit_Parameters[:, 2]), np.std(Fit_Parameters[:, 2])
        Area1_Average, StDev_Final_Area1 = np.average(Area_1), np.std(Area_1)
        width1[i] = 2*Average_Width*(1 + np.sqrt(1.386294361))
        area1[i] = Area1_Average
        Fits1[i, :] = [Average_Intensity, StDev_Intensity, Average_Position, StDev_Position, Average_Width, StDev_Width,  Area1_Average, StDev_Final_Area1]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_Intensity, Average_Position, Average_Width]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()

## This begins the section that will deconvolve two lines next to eachother
## and return parameters for either the left or the right line.
if Peaks1 == 2:
    Width_Array1, Width_Array2, Intensity_Array1, Intensity_Array2, Position1, Position2 = Get_Params_Two(*Extreme_Plot(Number_of_Files))
    Fit_Parameters = np.zeros((((len(Width_Array1)/2)*len(Intensity_Array1)), 6))
    Area_1 = np.zeros(((len(Width_Array1)/2)*len(Intensity_Array1), 2))
    which_peak = raw_input('Which peak do you want to look at, left or right?\n')
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left1)|(Wavelength>right1)] = 1
        tick = 0
        for intensity in range(0, len(Intensity_Array1)):
            for width in range(0, len(Width_Array1)/2):
                popt, pcov = so.curve_fit(two_peaks, Wavelength, Transmission, [Intensity_Array1[intensity], Position1, (np.abs(Width_Array1[width] - Width_Array1[(len(Width_Array1)-1)-width])), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width]))], maxfev = 100000)
                if any(popt < 0):
                    popt = [Intensity_Array1[intensity], Position1, (np.abs(Width_Array1[width] - Width_Array1[(len(Width_Array1)-1)-width])), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width]))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick, :] = Area_Two(Wavelength, Transmission, popt, Window_Length, Real_Wavelength)
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array1)/2)*len(Intensity_Array1))))*100) + "% for complete for file " + str(i+1)
        AreaL, StDevL, AreaR, StDevR = np.average(Area_1[:, 0]), np.std(Area_1[:, 0]), np.average(Area_1[:, 1]), np.std(Area_1[:, 1])
        Average_WidthL, StDevWidthL, Average_WidthR, StDevWidthR = np.average(Fit_Parameters[:, 2]), np.std(Fit_Parameters[:, 2]), np.average(Fit_Parameters[:, 5]), np.std(Fit_Parameters[:, 5])
        Average_IntensityL, StDevIntensityL, Average_IntensityR, StDevIntensityR = np.average(Fit_Parameters[:, 0]), np.std(Fit_Parameters[:, 0]), np.average(Fit_Parameters[:, 3]), np.std(Fit_Parameters[:, 3])
        Average_PositionL, StDevPositionL, Average_PositionR, StDevPositionR = np.average(Fit_Parameters[:, 1]), np.std(Fit_Parameters[:, 1]), np.average(Fit_Parameters[:, 4]), np.std(Fit_Parameters[:, 4])
        print "File " + str(i+1) + " complete"
        if which_peak == 'left':
            width1[i] = 2*Average_WidthL*(1 + np.sqrt(1.386294361))
            area1[i] = AreaL
            Fits1[i, :] = [Average_IntensityL, StDevIntensityL, Average_PositionL, StDevPositionL, Average_WidthL, StDevWidthL, AreaL, StDevL]
        elif which_peak == 'right':
            width1[i] = 2*Average_WidthR*(1 + np.sqrt(1.386294361))
            area1[i] = AreaR
            Fits1[i, :] = [Average_IntensityR, StDevIntensityR, Average_PositionR, StDevPositionR, Average_WidthR, StDevWidthR, AreaR, StDevR]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_IntensityL, Average_PositionL, Average_WidthL]))
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_IntensityR, Average_PositionR, Average_WidthR]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()

## This part of the script will deconvolute three lines next to each other
## and return the middle line parameters.
elif Peaks1 == 3:
    Width_Array1, Width_Array2, Width_Array3, Intensity_Array1, Intensity_Array2, Intensity_Array3, Position1, Position2, Position3 = Get_Params_Three(*Extreme_Plot(Number_of_Files))
    Fit_Parameters = np.zeros((((len(Width_Array2)/2)*len(Intensity_Array2)), 9))
    Area_1 = np.zeros(((len(Width_Array2)/2)*len(Intensity_Array2)))
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left1)|(Wavelength>right1)] = 1
        tick = 0
        for intensity in range(0, len(Intensity_Array2)):
            for width in range(0, len(Width_Array2)/2):
                popt, pcov = so.curve_fit(three_peaks, Wavelength, Transmission, [Intensity_Array1, Position1, (np.abs(Width_Array1)), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width])), Intensity_Array3, Position3, (np.abs(Width_Array3))], maxfev = 100000)
                if any(popt < 0):
                    popt  = [Intensity_Array1, Position1, (np.abs(Width_Array1)), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width])), Intensity_Array3, Position3, (np.abs(Width_Array3))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick] = Area_Three(Wavelength, Transmission, popt, Window_Length, Real_Wavelength)
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array2)/2)*len(Intensity_Array2))))*100) + "% for complete for file " + str(i+1)
        Average_Width, StDevWidth, Average_Intensity, StDevIntensity, Average_Position, StDevPosition = np.average(Fit_Parameters[:, 5]), np.std(Fit_Parameters[:, 5]), np.average(Fit_Parameters[:, 3]), np.std(Fit_Parameters[:, 3]), np.average(Fit_Parameters[:, 4]), np.std(Fit_Parameters[:, 4])
        Average_Area, StDevArea = np.average(Area_1), np.std(Area_1)
        print "File " + str(i+1) + " complete"
        width1[i] = 2*Average_Intensity*(1 + np.sqrt(1.386294361))
        area1[i] = Average_Area
        Fits1[i, :] = [Average_Intensity, StDevIntensity, Average_Position, StDevPosition, Average_Width, StDevWidth,  Average_Area, StDevArea]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_Intensity, Average_Position, Average_Width]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()


## Start of fitting for the second peak.
if Peaks2 == 1:
    Width_Array, Intensity_Array, Position = Get_Params_One(*Extreme_Plot(Number_of_Files))
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left2)|(Wavelength>right2)] = 1
        tick = 0
        Fit_Parameters = np.zeros((((len(Width_Array)/2)*len(Intensity_Array)), 3))
        Area_1 = np.zeros(((len(Width_Array)/2)*len(Intensity_Array)))
        for intensity in range(0, len(Intensity_Array)):
            for width in range(0, len(Width_Array)/2):
                popt, pcov = so.curve_fit(voigt1, Wavelength, Transmission, [Intensity_Array[intensity], Position, (np.abs(Width_Array[width] - Width_Array[(len(Width_Array)-1)-width]))], maxfev = 100000)
                if any(popt < 0):
                    popt = [Intensity_Array[intensity], Position, (np.abs(Width_Array[width] - Width_Array[(len(Width_Array)-1)-width]))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick] = Area_One(Wavelength, Transmission, Window_Length, Real_Wavelength, *popt)    
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array)/2)*len(Intensity_Array))))*100) + "% for complete for file " + str(i+1)
        print "File " + str(i+1) + " of " + str(Number_of_Files) + " complete."
        Average_Intensity, StDev_Intensity = np.average(Fit_Parameters[:, 0]), np.std(Fit_Parameters[:, 0])
        Average_Position, StDev_Position = np.average(Fit_Parameters[:, 1]), np.std(Fit_Parameters[:, 1])
        Average_Width, StDev_Width = np.average(Fit_Parameters[:, 2]), np.std(Fit_Parameters[:, 2])
        Area1_Average, StDev_Final_Area1 = np.average(Area_1), np.std(Area_1)
        width2[i] = 2*Average_Width*(1 + np.sqrt(1.386294361))
        area2[i] = Area1_Average
        Fits2[i, :] = [Average_Intensity, StDev_Intensity, Average_Position, StDev_Position, Average_Width, StDev_Width, Area1_Average, StDev_Final_Area1]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_Intensity, Average_Position, Average_Width]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()

## This begins the section that will deconvolve two lines next to eachother
## and return parameters for either the left or the right line.
if Peaks2 == 2:
    Width_Array1, Width_Array2, Intensity_Array1, Intensity_Array2, Position1, Position2 = Get_Params_Two(*Extreme_Plot(Number_of_Files))
    Fit_Parameters = np.zeros((((len(Width_Array1)/2)*len(Intensity_Array1)), 6))
    Area_1 = np.zeros(((len(Width_Array1)/2)*len(Intensity_Array1), 2))
    which_peak = raw_input('Which peak do you want to look at, left or right?\n')
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left2)|(Wavelength>right2)] = 1
        tick = 0
        for intensity in range(0, len(Intensity_Array1)):
            for width in range(0, len(Width_Array1)/2):
                popt, pcov = so.curve_fit(two_peaks, Wavelength, Transmission, [Intensity_Array1[intensity], Position1, (np.abs(Width_Array1[width] - Width_Array1[(len(Width_Array1)-1)-width])), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width]))], maxfev = 100000)
                if any(popt < 0):
                    popt = [Intensity_Array1[intensity], Position1, (np.abs(Width_Array1[width] - Width_Array1[(len(Width_Array1)-1)-width])), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width]))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick, :] = Area_Two(Wavelength, Transmission, popt, Window_Length, Real_Wavelength)
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array1)/2)*len(Intensity_Array1))))*100) + "% for complete for file " + str(i+1)
        AreaL, StDevL, AreaR, StDevR = np.average(Area_1[:, 0]), np.std(Area_1[:, 0]), np.average(Area_1[:, 1]), np.std(Area_1[:, 1])
        Average_WidthL, StDevWidthL, Average_WidthR, StDevWidthR = np.average(Fit_Parameters[:, 2]), np.std(Fit_Parameters[:, 2]), np.average(Fit_Parameters[:, 5]), np.std(Fit_Parameters[:, 5])
        Average_IntensityL, StDevIntensityL, Average_IntensityR, StDevIntensityR = np.average(Fit_Parameters[:, 0]), np.std(Fit_Parameters[:, 0]), np.average(Fit_Parameters[:, 3]), np.std(Fit_Parameters[:, 3])
        Average_PositionL, StDevPositionL, Average_PositionR, StDevPositionR = np.average(Fit_Parameters[:, 1]), np.std(Fit_Parameters[:, 1]), np.average(Fit_Parameters[:, 4]), np.std(Fit_Parameters[:, 4])
        print "File " + str(i+1) + " complete"
        if which_peak == 'left':
            width2[i] = 2*Average_WidthL*(1 + np.sqrt(1.386294361))
            area2[i] = AreaL
            Fits2[i, :] = [Average_IntensityL, StDevIntensityL, Average_PositionL, StDevPositionL, Average_WidthL, StDevWidthL, AreaL, StDevL]
        elif which_peak == 'right':
            width2[i] = 2*Average_WidthR*(1 + np.sqrt(1.386294361))
            area2[i] = AreaR
            Fits2[i, :] = [Average_IntensityR, StDevIntensityR, Average_PositionR, StDevPositionR, Average_WidthR, StDevWidthR, AreaR, StDevR]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_IntensityL, Average_PositionL, Average_WidthL]))
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_IntensityR, Average_PositionR, Average_WidthR]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()

## This part of the script will deconvolute three lines next to each other
## and return the middle line parameters.
elif Peaks2 == 3:
    Width_Array1, Width_Array2, Width_Array3, Intensity_Array1, Intensity_Array2, Intensity_Array3, Position1, Position2, Position3 = Get_Params_Three(*Extreme_Plot(Number_of_Files))
    Fit_Parameters = np.zeros((((len(Width_Array2)/2)*len(Intensity_Array2)), 9))
    Area_1 = np.zeros(((len(Width_Array2)/2)*len(Intensity_Array2)))
    for i in range(Number_of_Files):
        Wavelength = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0))
        Transmission = np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1))
        Transmission[(Wavelength<left2)|(Wavelength>right2)] = 1
        tick = 0
        for intensity in range(0, len(Intensity_Array2)):
            for width in range(0, len(Width_Array2)/2):
                popt, pcov = so.curve_fit(three_peaks, Wavelength, Transmission, [Intensity_Array1, Position1, (np.abs(Width_Array1)), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width])), Intensity_Array3, Position3, (np.abs(Width_Array3))], maxfev = 100000)
                if any(popt < 0):
                    popt = [Intensity_Array1, Position1, (np.abs(Width_Array1)), Intensity_Array2[intensity], Position2, (np.abs(Width_Array2[width] - Width_Array2[(len(Width_Array2)-1)-width])), Intensity_Array3, Position3, (np.abs(Width_Array3))]
                Fit_Parameters[tick, :] = popt
                Area_1[tick] = Area_Three(Wavelength, Transmission, popt, Window_Length, Real_Wavelength)
                tick = tick + 1
                print str((float(tick)/float(((len(Width_Array2)/2)*len(Intensity_Array2))))*100) + "% for complete for file " + str(i+1)
        Average_Width, StDevWidth, Average_Intensity, StDevIntensity, Average_Position, StDevPosition = np.average(Fit_Parameters[:, 5]), np.std(Fit_Parameters[:, 5]), np.average(Fit_Parameters[:, 3]), np.std(Fit_Parameters[:, 3]), np.average(Fit_Parameters[:, 4]), np.std(Fit_Parameters[:, 4])
        Average_Area, StDevArea = np.average(Area_1), np.std(Area_1)
        print "File " + str(i+1) + " complete"
        width2[i] = 2*Average_Intensity*(1 + np.sqrt(1.386294361))
        area2[i] = Average_Area
        Fits2[i, :] = [Average_Intensity, StDevIntensity, Average_Position, StDevPosition, Average_Width, StDevWidth, Average_Area, StDevArea]
    plt.plot(Wavelength, voigt1(Wavelength, *[Average_Intensity, Average_Position, Average_Width]))
    plt.plot(np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (0)), np.genfromtxt(File_Location + '/PRISMSPECT' + str(i) + '.ppd', usecols = (1)))
    plt.show()

## Make the arrays used for plotting--namely the width arrays
## and ratio array.  The array must be reshaped so that it is
## a 2-dimensional array and that the rows and column correspond
## to constant dens or temp data.
Ratio_Matrix = (np.array(area1)/np.array(area2)).reshape((Number_of_Temperatures, Number_of_Densities))
Width_Matrix1 = np.array(width1).reshape((Number_of_Temperatures, Number_of_Densities))
Width_Matrix2 = np.array(width2).reshape((Number_of_Temperatures, Number_of_Densities))
Intensities_Matrix = np.array(Fits1[:, 0]/Fits2[:, 0]).reshape((Number_of_Temperatures, Number_of_Densities))

## Plot the stuff.  Will plot ratio vs temp and tensity as well as
## width of the two lines vs temp and density.
root = tk.Tk()  # Has to do with making a window NOT pop up
root.withdraw()
Directory = fd.askdirectory(title = 'Data save drectory') # ask where to save
root.destroy()

## Find out if we are fitting for temperature ratios or density widths.
Purpose = raw_input('Enter fitting purpose--(temp or dens):\n')

if Purpose == 'temp':
    ## Ratio vs. density with lines of constant temp
    plt.figure(figsize = (16,16))
    for i in range(0, Number_of_Temperatures):
        Ratio_Constant_Temp = Ratio_Matrix[i, :]
        plt.plot(Density_Matrix, Ratio_Constant_Temp, label = str(float(int(Temperature_Matrix[i]))) + ' eV', color = colors[i], marker = 's', ms = 15, linewidth = 3)
    plt.xlabel('Density (ions/cc)', fontsize = 20)
    plt.ylabel('Ratio', fontsize = 20)
    plt.tick_params(which = 'both', labelsize = 18)
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, prop = {'size':25})
    Plot_Info = np.concatenate((Density_Matrix.reshape(Number_of_Densities, 1), Ratio_Constant_Temp.reshape(Number_of_Densities, 1)
    ), axis = 1)
    np.savetxt(Directory + "/RATIO_V_DENSITY_CONSTANT_TEMP" + str(float(int(Temperature_Matrix[i]))) + "eV.txt", Plot_Info)
    plt.savefig(Directory + '/RATIO_V_DENSITY_CONSTANT_TEMP', bbox_inches = 'tight')

    ## Ratio vs. temperature with lines of constant dens
    plt.figure(figsize = (16,16))
    for i in range(0, Number_of_Densities):
        Ratio_Constant_Dens = Ratio_Matrix[:, i]
        plt.plot(Temperature_Matrix, Ratio_Constant_Dens, label = str(float(int(Density_Matrix[i]))) + ' ions/cc', color = colors[i], marker = 's', ms = 15, linewidth = 3)
    plt.xlabel('Temperature (eV)', fontsize = 20)
    plt.ylabel('Ratio', fontsize = 20)
    plt.tick_params(which = 'both', labelsize = 18)
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, prop = {'size':25})
    Plot_Info = np.concatenate((Temperature_Matrix.reshape(Number_of_Temperatures, 1), Ratio_Constant_Dens.reshape(Number_of_Temperatures, 1)), axis = 1)
    np.savetxt(Directory + "/RATIO_V_TEMPERATURE_CONSTANT_DENS" + str(float(int(Density_Matrix[i]))) + "ionscc.txt", Plot_Info)
    plt.savefig(Directory + '/RATIO_V_TEMPERATURE_CONSTANT_DENS', bbox_inches = 'tight')

    plt.figure(figsize = (16,16))
    for i in range(0, Number_of_Densities):
        Intensity_Constant_Dens = Intensities_Matrix[:, i]
        plt.plot(Temperature_Matrix, Intensity_Constant_Dens, label = str(float(int(Density_Matrix[i]))) + ' ions/cc', color = colors[i], marker = 's', ms = 15, linewidth = 3)
    plt.xlabel('Temperature (eV)', fontsize = 20)
    plt.ylabel('Intensity Ratio', fontsize = 20)
    plt.tick_params(which = 'both', labelsize = 18)
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, prop = {'size':25})
    Plot_Info = np.concatenate((Temperature_Matrix.reshape(Number_of_Temperatures, 1), Intensity_Constant_Dens.reshape(Number_of_Temperatures, 1)), axis = 1)
    np.savetxt(Directory + "/INTENSITY_V_TEMPERATURE_CONSTANT_DENS" + str(float(int(Density_Matrix[i]))) + "ionscc.txt", Plot_Info)
    plt.savefig(Directory + '/INTENSITY_V_TEMPERATURE_CONSTANT_DENS', bbox_inches = 'tight')
elif Purpose == 'dens':
    ## Width vs density with lines of constant temp for first line.
    plt.figure(figsize = (16, 16))
    for i in range(0, Number_of_Temperatures):
        Width = Width_Matrix1[i, :]
        plt.plot(Density_Matrix, Width, label = str(float(int(Temperature_Matrix[i]))), color = colors[i], marker = 's', ms = 15, linewidth = 3)
    plt.xlabel('Density (ions/cc)')
    plt.ylabel('Width (Angstroms)')
    plt.title('Width vs Density for Line 1')
    plt.tick_params(which = 'both', labelsize = 18)
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, prop = {'size':25})
    Plot_Info = np.concatenate((Density_Matrix.reshape(Number_of_Densities, 1), Width.reshape(Number_of_Densities, 1)), axis = 1)
    np.savetxt(Directory + "/Width1_vs_Density_Constant_Temp" + str(float(int(Temperature_Matrix[i]))) + "eV.txt", Plot_Info)
    plt.savefig(Directory + "/Width1_v_Density_Constant_Temp", bbox_inches = 'tight')

    ## Width vs density with lines of constant temp for second line.
    plt.figure(figsize = (16, 16))
    for i in range(0, Number_of_Temperatures):
        Width = Width_Matrix2[i, :]
        plt.plot(Density_Matrix, Width, label = str(float(int(Temperature_Matrix[i]))), color = colors[i], marker = 's', ms = 15, linewidth = 3)
    plt.xlabel('Density (ions/cc)')
    plt.ylabel('Width (Angstroms)')
    plt.title('Width vs Density for Line 2')
    plt.tick_params(which = 'both', labelsize = 18)
    plt.legend(bbox_to_anchor = (1.05, 1), loc = 2, prop = {'size':25})
    Plot_Info = np.concatenate((Density_Matrix.reshape(Number_of_Densities, 1), Width.reshape(Number_of_Densities, 1)), axis = 1)
    np.savetxt(Directory + "/Width2_vs_Density_Constant_Temp" + str(float(int(Temperature_Matrix[i]))) + "eV.txt", Plot_Info)
    plt.savefig(Directory + "/Width2_v_Density_Constant_Temp", bbox_inches = 'tight')

## Finally to save text files that will probable be used in the next part of HADHES.
Line_Parameters = np.concatenate((np.array(Temperatures).reshape(len(Temperatures), 1), np.array(Densities).reshape(len(Temperatures), 1), Fits1, Fits2, np.array(width1).reshape(len(Temperatures), 1), np.array(width2).reshape(len(Temperatures), 1), (area1/area2).reshape(len(Temperatures), 1)), axis = 1)

## Ask the name of the file to be outputted.
filename = raw_input('Enter the name of the output file:\n')

np.savetxt(Directory + '/' + filename, Line_Parameters, fmt = '%1.2e', header = 'Output file from PrismSPECT_Fitting.py\nTemperatures Densities Intensity 1,  St Dev Intensity 1, Position 1, St Dev Position 1, Width 1, St Dev Width 1, Area 1, St Dev Area 1, Intensity 2, St Dev Intensity 2, Position 2, St Dev Position 2, Width 2, St Dev Width 2, Area 2, St Dev Area 2, FWHM1, FWHM2\n')

Continue = raw_input('PrismSPECT_Fitting.py has finished running--press any letter key to continue.\n')
del Continue
