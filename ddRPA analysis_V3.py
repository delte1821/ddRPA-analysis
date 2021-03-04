# Ji Wook Choi - BNTL Sogang University - 2021

###########################################################################################

# Importing required packages

import numpy as np
import cv2
import plotly.graph_objs as go
import plotly as py
import math
from operator import truediv
import sqlite3

# ==============================================================================================================
# Defining Variables can be changed by user
Nimg    = 8             # Number of images
ImgMag = 10             # Magnification of images (4X = 4, 10X = 10)

# Convert um to px
if ImgMag == 4:
    UP = 0.453          # UP: Constant of the ratio of um:px
elif ImgMag ==10:
    UP = 1.12
else:
    UP = 0.453

# Variables for circle detection
dp1         = 1     # Ratio of accumulator, one means same resolution as input image, bigger numbers mean image is reduced 
param11     = 10     # Threshold passed to Canny edge detector
param21     = 10     # Accumulatoir threshold for circles, the lower the more false circles are recognised
minDiameter1  = 80    # Minimum radius of found circles in um
maxDiameter1  = 140     # Maximum radius of found circles in um
minDist1    = minDiameter1    # Minimum distance between two circle centers in um   

# Defining appendixes required for opening images
name1 = "R"         # Appendix for reference dye images
name2 = ".jpg"      # Appendix for image filetype
tag   = "-marked"   # Appendix for images in which found circles are marked

# ===============================================================================================================
# Defining constants that should not be changed by the user
Npos = 0.0 # Counter for positive wells, float so probability can be calculated
Nneg = 0.0 # Counter for negative wells, float so probability can be calculated

# ===============================================================================================================

# Function for Callculating flourescence intensity inside of cicle

def CircleIntensity(centerx,centery,radius,image, color):
    # Callculates the intensity indside of a circle
    #
    # Parameters:
    # centerx = x-coordinte of circle
    # centery = y-coordinate of circle
    # radius = Radius of circle
    # image = image for callculating intensity (brightness saved as 8 bit, i.e. from 0 to 255
    #
    # Returns:
    # Intensity = average intensity value of circle in the tested image

    if color == "B":
        coval = 0
    elif color == "G":
        coval = 1
    elif color == "R":
        coval = 1
    else:
        print("Color not RGB, assuming color is green")
        coval = 1
    
    # Definging required parameters
    npixels = 0.0     # Count for pixels to find average
    brightness = 0.0  # Count for brightness of pixel

    # Creating square around circle
    #if (centerx >=2 and centery >=2 and centerx <= image.shape[0]-1
    for x in range(centerx-radius-2,centerx+radius+2):        # Varying through x
            for y in range(centery-radius-2,centery+radius+2):       # Varying through y
                 if (x <= image.shape[1]-1 and y <= image.shape[0]-1 and x>=0 and y>=0):  # Making sure coordinate in image
                     pixeldistance = math.sqrt((centerx - x)**2 + (centery - y)**2)     # Pythagoras to find radius from iterated pixcel to center of circle
                     if pixeldistance <  radius:                                        # If Pixel is in circle add to intensity callculation
                         pixel = image[y,x]
                         brightness = brightness + float(pixel[coval])/255                  # Updating total brightness
                         npixels = npixels + 1                                          # Updating total pixcel count
    if npixels == 0:
        npixels = 1 # Preventing error, division by zero
        
    Intensity = brightness / npixels                                                    # Callculating average intesnity of circle

    if Intensity == 0:      # Preventing division by zero
        Intensity = 0.00000001 

    return Intensity

# ==============================================================================================================

# Function for plotting histograms
def Histogram(Data):

    # Making histogram
    data = [
        go.Histogram(
            x=Data, xbins=dict(start=0, size=0.005, end=1.01) 
        )
    ]

    layout = go.Layout(
        title='Histogram of flourescence intensity of wells',
        xaxis=dict(title='Flourescence intensity',
                   range=[0,1.01]),
        yaxis=dict(title='Well count'), 
        )

    fig = go.Figure(data=data, layout=layout)

    py.offline.plot(fig)

# ===============================================================================================================

# Function for esimating Concentration

def ConcCallculation(pHat,Npart,Vol):
    # Function callculating the concentration of outcome of dPCR
    #
    # Parameters:
    # pHat  =   estimated propability of positive partition
    # Npart =   total number of partitions
    # Vol   =   Volume of partitions in uL
    #
    # Returns:
    # C_est = callculated concentration in #particles/uL
    # C_low = lower confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    # C_upp = upper confidence intervall of calculated concentration (95% z-distribution) in # particles/uL
    #
    #######################################################

    # Defingin constants
    zc = 1.96   # 95% confidence intervall z-distribution

    # Callculation of confidence interval on pHat
    pHat_Dev = zc * math.sqrt((pHat * (1-pHat))/Npart)  # Deviation on expected result
    p_hat_low = pHat - pHat_Dev  # Lower bound of p_hat
    p_hat_upp = pHat + pHat_Dev  # Upper bound of p_hat

    # Callculating mean number of molecules per patition including 95%
    # confidence intervall
    lambda1 = -math.log(1-pHat)     # average number of molecules per division as per Poission distribution
    lambda_low = -math.log(1-p_hat_low)  # lower bound of average number of molecules per division
    lambda_upp = -math.log(1-p_hat_upp)  # upper bound of average number of molecules per division

    # Callculating concentrations in mol/uL from lambda values including
    # confidence intervalls
    C_est = lambda1 / Vol       # Esitmated concentration
    C_low = lambda_low / Vol    # Estimated lower bound of concentration
    C_upp = lambda_upp / Vol    # Estimated higher bound of concentration

    return C_est, C_low, C_upp

# ==============================================================================================================

# Main code file calling the previous defined functions

# Creating database for Saving results
conn = sqlite3.connect('OutputData.db')
c = conn.cursor()
c.execute('DROP TABLE IF EXISTS Wells')
c.execute('DROP TABLE IF EXISTS Concentration')
c.execute('CREATE TABLE IF NOT EXISTS Wells (ImmageNumber, Xcoordinate, Ycoordinate, Radius, Intensity)')
c.execute('CREATE TABLE IF NOT EXISTS Concentration (Threshold, PosWells, NegWells, Propability, LowBound, UppBound, Conc)')
conn.commit()

#################################################################################################################
# Loop going through all the images taken, identifying circles, drawing cicles, and measuring intensities

for i in range(0,Nimg):
    
    # Print current image for debugging purposes
    

    # Creating filenames for images to be opened
    nameFlu = name1 + str(i+1) + name2
    nameFlutag = str(i+1) + name1 + tag + name2
    print("Currently on ", nameFlu)

    # Converting images for detection of circles and measurement of flourescence

    # Fluorescence dye image
    imgFlu  = cv2.imread(nameFlu, 0)       # Opening image
    imgFlu  = cv2.medianBlur(imgFlu,5)  # Smothening image data
    imgFlu2 = cv2.imread(nameFlu)       # Image for marking circles in
    imgFlu3 = cv2.imread(nameFlu)       # Measuring reference dye intesity in

    # Fitting circles using Hough transform from package cv2
    minRadius1 = int(UP * minDiameter1 / 2)
    maxRadius1 = int(UP * maxDiameter1 / 2)
    minDist1 = int(UP * minDist1)
    # function calls as follows: cv2.HoughCircles(image, method, dp, minDist, circles, param1, param2, minRadius, maxRadius)
    circles = cv2.HoughCircles(imgFlu,cv2.HOUGH_GRADIENT,dp1,minDist1,param1=param11,param2=param21,minRadius=minRadius1,maxRadius=maxRadius1)

    # Drawing fitted circles to reference dye image
    circles = np.uint16(np.around(circles)) # Preparing data for plotting
    for j in circles[0,:]:
        # draw the outer circle
        cv2.circle(imgFlu2,(j[0],j[1]),j[2],(0,255,0),2)
        # draw the center of the circle
        cv2.circle(imgFlu2,(j[0],j[1]),2,(0,0,255),3)
    # Saving image with marked circles
    cv2.imwrite(nameFlutag,imgFlu2)

    # Measuring intensity in marked cicles
    intensity = []  # Empty list to save circle intensities into
    for j in circles[0,:]:
        intensity.append(CircleIntensity(j[0],j[1],j[2],imgFlu3,name1))
    
 
    # Saving data to output database
    k=0     # Running index to sync intensity list with circle list
    for j in circles[0,:]:
        c.execute('''INSERT INTO Wells VALUES (?, ?, ?, ?, ?)''', (str(i+1) ,str(j[0]) ,str(j[1]) ,str(j[2]) ,str(intensity[k]))) 
        k = k+1

conn.commit()

#################################################################################################################
# Loop going through all the images taken, identifying circles, drawing cicles, and measuring intensities

# Importing data from database
data = c.execute('''Select * FROM Wells''')
data = [float(item[4]) for item in data.fetchall()]

# Mapping data
plotdata = list(map(float, data))

# Plotting histogram
Histogram(plotdata)
while True:
    
    # Defining manual threshold
    Threshold1 = eval(input('Please enter minimum threshold: '))
    if Threshold1 == 4:
        break
    Threshold2 = eval(input('Please enter maximum threshold: '))

    # Counting positive and negative partitions
    for x in plotdata:
        if (Threshold1 < x) & (x <= Threshold2):        # Wells contianing now particle/cell have a higher intensity
            Nneg = Nneg + 1        
        elif (x < 1.0) & (x > Threshold2):         # If well is positive average brightness is low
            Npos = Npos + 1


    # Calculating volume of well
    dropletradius = UP * (minRadius1 + maxRadius1) / 2
    VolWell = (4/3) * math.pi * dropletradius ** 3 # Volume of a well in cubic micrometer
    VolWell = VolWell * math.pow(10,-9)

    # Callculation of concentrations according to Poission distribution
    pHat = Npos / (Npos + Nneg) # Estimated propability of positive well
    Npart = Npos + Nneg # Total number of detected wells
    C_est, C_low, C_upp = ConcCallculation(pHat,Npart,VolWell) # Callculation of concentrations according to predifined function
    Stdev = C_est - C_low
    CV = 100 * Stdev / C_est

    print(('Radius of droplet:                       {0:.2f}' .format(dropletradius)))
    print(('Number of Poitive calls:                 {0:.2f}' .format(Npos)))
    print(('Number of Negative calls:                {0:.2f}' .format(Nneg)))
    print(('Propability of positive partition:       {0:.2f} \n' .format(pHat)))
    # Outputting results
    print(('Estimated concentration:                 {0:.3f} copies/uL' .format(C_est)))
    print(('Standard deviation:                      {0:.3f} copies/uL' .format(Stdev)))
    print(('Coefficient of variation (CV):           {0:.3f} %' .format(CV)))
    print(('Lower bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_low)))
    print(('Upper bound of 95% Confidence intervall: {0:.3f} copies/uL' .format(C_upp)))

    c.execute(''' INSERT INTO Concentration VALUES (?, ?, ?, ?, ?, ?, ?)''', (Threshold2 , Npos, Nneg, pHat, C_low, C_upp, C_est))

    Npos = 0.0  # Initializing parameter
    Nneg = 0.0  # Initializing parameter


conn.commit()
conn.close()



    
