import sys
import csv
import math
import cmath
from copy import deepcopy
import numpy as np
from scipy.signal import savgol_filter

#****************** Helper functions**********************
def smooth_col_2_3(datatable, window = 5, order = 1, re_normalize = True):
    real_list = []
    imag_list = []
    for i in range(0,len(datatable)):
        mag = datatable[i][1]
        phi = datatable[i][2] / 180 * math.pi
        real_list.append(mag * math.cos(phi))
        imag_list.append(mag * math.sin(phi))
    real_list = savgol_filter(real_list, window, order)
    imag_list = savgol_filter(imag_list, window, order)
    for i in range(0,len(datatable)):
        x = complex(real_list[i], imag_list[i])
        datatable[i][1] = abs(x)
        datatable[i][2] = cmath.phase(x) * 180 / math.pi
    if (re_normalize):
        normalize(datatable)

def normalize(table):
    BasePhase = table[-1][2]
    BaseMag = table[-1][1]
    for i in range(0,len(table)):
        table[i][1] /= BaseMag
        table[i][2] -= BasePhase

def parseRawCSV(filename):
    dataTable = []
    try:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
        
            for row in csv_reader:
                row_float = []
                for i in range(0, len(row)):
                    try:
                        row_float.append(float(row[i]))
                    except ValueError:
                        continue
                if len(row_float) > 0:
                    dataTable.append(row_float)
    except:
        return dataTable
    return dataTable

def writeCSV(filename, data):
    with open(filename, mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', lineterminator = '\n')
        csv_writer.writerow(['Frequency', 'Magnitude', 'Phase(degrees)'])
        for i in range(0, len(data)):
            csv_writer.writerow(data[i])
    csv_file.close()

def PolDivide(numerator, denominator):
    for i in range(0, len(numerator)):
        freq = numerator[i][0]
        numerator[i][1] /= interpolate(freq, denominator, 1)
        numerator[i][2] -= interpolate(freq, denominator, 2)

def PolInv(baseTable):
    ret = deepcopy(baseTable)
    for i in range(0, len(baseTable)):
        ret[i][1] = 1 / baseTable[i][1]
        ret[i][2] = -baseTable[i][2]
    return ret

def InvertPhase(baseTable):
    for i in range(0, len(baseTable)):
        baseTable[i][2] *= -1

def PolMult(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        baseTable[i][1] *= interpolate(freq, auxTable, 1)
        baseTable[i][2] += interpolate(freq, auxTable, 2)

def PolAdd(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = interpolate(freq, auxTable, 2) / 180 * math.pi
        aux_mag = interpolate(freq, auxTable, 1)
        aux = complex(aux_mag * math.cos(base_phi), aux_mag * math.sin(base_phi))
        sum = base + aux
        baseTable[i][1] = abs(sum)
        baseTable[i][2] = cmath.phase(sum) * 180 / math.pi

def PolSubtract(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = interpolate(freq, auxTable, 2) / 180 * math.pi
        aux_mag = interpolate(freq, auxTable, 1)
        aux = complex(aux_mag * math.cos(base_phi), aux_mag * math.sin(base_phi))
        diff = base - aux
        baseTable[i][1] = abs(diff)
        baseTable[i][2] = cmath.phase(diff) * 180 / math.pi

def interpolate(x, lookupTable, ycol = 1, semilog = True):
    #make sure table is sorted in descending order (by first column)
    if(lookupTable[0][0] < lookupTable[-1][0]):
        lookupTable.sort(key = sortByFirstElement)

    #if x is outside the bounds of the table, return one of the bounding values
    if x >= lookupTable[0][0]:
        return lookupTable[0][ycol]
    elif x <= lookupTable[-1][0]:
        return lookupTable[-1][ycol]
    #else interpolate
    else:
        indexR = 1
        for i in range(1, len(lookupTable)):
            if x > lookupTable[indexR][0]:
                break
            else:
                indexR += 1
        indexL = indexR - 1
        xL = 0.0; xR = 0.0; xi = x
        if (semilog):
            xL = math.log10(lookupTable[indexL][0])
            xR = math.log10(lookupTable[indexR][0])
            xi = math.log10(x)
        else:
            xL = lookupTable[indexL][0]
            xR = lookupTable[indexR][0]
        yL = lookupTable[indexL][ycol]
        yR = lookupTable[indexR][ycol]
        yi = yL + (xi - xL) / (xR - xL) * (yR - yL)
        return yi

def sortByFirstElement(val):
    return val[0]


def ScalarMult(baseTable, scalar):
    for i in range(0, len(baseTable)):
        baseTable[i][1] *= scalar

def popSweep(masterTable):
    ret = []
    while True:
        mag = 0
        phase = 0
        freq = 0
        IsigMag = 0
        N = 0
        while True:
            row = masterTable.pop(0)
            freq = row[0]
            mag += row[1]
            phase += row[2]
            if len(row) > 3:
                IsigMag += row[3]
            N += 1
            if (len(masterTable) == 0) or (masterTable[0][0] != freq):
                break
        #if len(row) > 3:
        #    ret.append([freq, mag / N, phase / N, IsigMag / N])
        #else:
        #    ret.append([freq, mag / N, phase / N])
        ret.append([freq, mag / N, phase / N])
        if (len(masterTable) == 0) or (masterTable[0][0] > freq):
            break
    return ret


#*************** Main script ****************************
basefolder = sys.argv[1]

# import raw data from Squidstat experiment file
rawDataTier1 = parseRawCSV(sys.argv[2])
InvertPhase(rawDataTier1)

rawDataTier2 = parseRawCSV(sys.argv[3])
InvertPhase(rawDataTier2)

rawDataTier3 = parseRawCSV(sys.argv[4])
InvertPhase(rawDataTier3)

rawData10R = parseRawCSV(sys.argv[5])
InvertPhase(rawData10R)

rawData100R = parseRawCSV(sys.argv[6])
InvertPhase(rawData100R)

rawData1k = parseRawCSV(sys.argv[7])
InvertPhase(rawData1k)

rawData100k = parseRawCSV(sys.argv[8])
InvertPhase(rawData100k)

rawData10M = parseRawCSV(sys.argv[9])
InvertPhase(rawData10M)

# parse and normalize sweeps from each gain setting
range2 = popSweep(rawDataTier1)
H_BaselineV = deepcopy(range2)
H_BaselineI = PolInv(H_BaselineV)
V2 = popSweep(rawDataTier1)
V5 = popSweep(rawDataTier1)
V10 = popSweep(rawDataTier1)
I2 = PolInv(popSweep(rawDataTier1))
I5 = PolInv(popSweep(rawDataTier1))
I10 = PolInv(popSweep(rawDataTier1))
V20 = popSweep(rawDataTier2)
V50 = popSweep(rawDataTier2)
V100 = popSweep(rawDataTier2)
I20 = PolInv(popSweep(rawDataTier2))
I50 = PolInv(popSweep(rawDataTier2))
I100 = PolInv(popSweep(rawDataTier2))
V200 = popSweep(rawDataTier3)
V500 = popSweep(rawDataTier3)
V1000 = popSweep(rawDataTier3)
I200 = PolInv(popSweep(rawDataTier3))
I500 = PolInv(popSweep(rawDataTier3))
I1000 = PolInv(popSweep(rawDataTier3))
normalize(H_BaselineV)
normalize(H_BaselineI)
normalize(V2)
normalize(V5)
normalize(V10)
normalize(V20)
normalize(V50)
normalize(V100)
normalize(V200)
normalize(V500)
normalize(V1000)
normalize(I2)
normalize(I5)
normalize(I10)
normalize(I20)
normalize(I50)
normalize(I100)
normalize(I200)
normalize(I500)
normalize(I1000)

#calculate tier 1 gains (x2 thru x10)
PolDivide(V2, H_BaselineV)
PolDivide(V5, H_BaselineV)
PolDivide(V10, H_BaselineV)
PolDivide(I2, H_BaselineI)
PolDivide(I5, H_BaselineI)
PolDivide(I10, H_BaselineI)

#calculate tier2 gains (x20 thru x100)
PolDivide(V20, H_BaselineV)
PolMult(V20, I10)
PolDivide(V50, H_BaselineV)
PolMult(V50, I10)
PolDivide(V100, H_BaselineV)
PolMult(V100, I10)
PolDivide(I20, H_BaselineI)
PolMult(I20, V10)
PolDivide(I50, H_BaselineI)
PolMult(I50, V10)
PolDivide(I100, H_BaselineI)
PolMult(I100, V10)

#calculate tier2 gains (x200 thru x1000)
PolDivide(V200, H_BaselineV)
PolMult(V200, I100)
PolDivide(V500, H_BaselineV)
PolMult(V500, I100)
PolDivide(V1000, H_BaselineV)
PolMult(V1000, I100)
PolDivide(I200, H_BaselineI)
PolMult(I200, V100)
PolDivide(I500, H_BaselineI)
PolMult(I500, V100)
PolDivide(I1000, H_BaselineI)
PolMult(I1000, V100)

#smooth data
smooth_col_2_3(I2)
smooth_col_2_3(I5)
smooth_col_2_3(I10)
smooth_col_2_3(I20)
smooth_col_2_3(I50)
smooth_col_2_3(I100)
smooth_col_2_3(I200)
smooth_col_2_3(I500)
smooth_col_2_3(I1000)
smooth_col_2_3(V2)
smooth_col_2_3(V5)
smooth_col_2_3(V10)
smooth_col_2_3(V20)
smooth_col_2_3(V50)
smooth_col_2_3(V100)
smooth_col_2_3(V200)
smooth_col_2_3(V500)
smooth_col_2_3(V1000)

# parse sweeps from each range setting
R2_V10_I1_10R = popSweep(rawData10R)
R1_V10_I1_10R = popSweep(rawData10R)
R1_V10_I10_10R = popSweep(rawData10R)
R0_V10_I10_10R = popSweep(rawData10R)
R2_V10_I1_100R = popSweep(rawData100R)
R3_V10_I1_100R = popSweep(rawData100R)
R3_V10_I1_1k = popSweep(rawData1k)
R4_V10_I1_1k = popSweep(rawData1k)
R5_V10_I1_1k = popSweep(rawData1k)
R5_V10_I1_100k = popSweep(rawData100k)
R6_V10_I1_100k = popSweep(rawData100k)
R6_V10_I100_10M = popSweep(rawData10M)
R8_V10_I1_10M = popSweep(rawData10M)


#range2
normalize(range2)
range2 = PolInv(range2)

#range1
range1 = deepcopy(R1_V10_I1_10R)
PolDivide(range1, R2_V10_I1_10R)
normalize(range1)
range1 = PolInv(range1)
PolMult(range1, range2)

#range0
range0 = deepcopy(R0_V10_I10_10R)
PolDivide(range0, R1_V10_I10_10R)
normalize(range0)
range0 = PolInv(range0)
PolMult(range0, range1)

#range3
range3 = deepcopy(R3_V10_I1_100R)
PolDivide(range3, R2_V10_I1_100R)
normalize(range3)
range3 = PolInv(range3)
PolMult(range3, range2)

#range4
range4 = deepcopy(R4_V10_I1_1k)
PolDivide(range4, R3_V10_I1_1k)
normalize(range4)
range4 = PolInv(range4)
PolMult(range4, range3)

#range5
range5 = deepcopy(R5_V10_I1_1k)
PolDivide(range5, R3_V10_I1_1k)
normalize(range5)
range5 = PolInv(range5)
PolMult(range5, range3)

#range6/7
range6 = deepcopy(R6_V10_I1_100k)
normalize(R5_V10_I1_100k)
smooth_col_2_3(R5_V10_I1_100k)
PolDivide(range6, R5_V10_I1_100k)
normalize(range6)
range6 = PolInv(range6)
PolMult(range6, range5)

#range8/9
range8 = deepcopy(R8_V10_I1_10M)
normalize(R6_V10_I100_10M)
smooth_col_2_3(R6_V10_I100_10M)
PolDivide(range8, R6_V10_I100_10M)
normalize(range8)
range8 = PolInv(range8)
PolMult(range8, range6)

# write the resulting data to CSV files
writeCSV(basefolder + 'WEgain1_Igain2.csv', I2)
writeCSV(basefolder + 'WEgain1_Igain5.csv', I5)
writeCSV(basefolder + 'WEgain1_Igain10.csv', I10)
writeCSV(basefolder + 'WEgain1_Igain20.csv', I20)
writeCSV(basefolder + 'WEgain1_Igain50.csv', I50)
writeCSV(basefolder + 'WEgain1_Igain100.csv', I100)
writeCSV(basefolder + 'WEgain1_Igain200.csv', I200)
writeCSV(basefolder + 'WEgain1_Igain500.csv', I500)
writeCSV(basefolder + 'WEgain1_Igain1000.csv', I1000)
writeCSV(basefolder + 'WEgain2_Igain1.csv', V2)
writeCSV(basefolder + 'WEgain5_Igain1.csv', V5)
writeCSV(basefolder + 'WEgain10_Igain1.csv', V10)
writeCSV(basefolder + 'WEgain20_Igain1.csv', V20)
writeCSV(basefolder + 'WEgain50_Igain1.csv', V50)
writeCSV(basefolder + 'WEgain100_Igain1.csv', V100)
writeCSV(basefolder + 'WEgain200_Igain1.csv', V200)
writeCSV(basefolder + 'WEgain500_Igain1.csv', V500)
writeCSV(basefolder + 'WEgain1000_Igain1.csv', V1000)
writeCSV(basefolder + 'Range0.csv', range0)
writeCSV(basefolder + 'Range1.csv', range1)
writeCSV(basefolder + 'Range2.csv', range2)
writeCSV(basefolder + 'Range3.csv', range3)
writeCSV(basefolder + 'Range4.csv', range4)
writeCSV(basefolder + 'Range5.csv', range5)
writeCSV(basefolder + 'Range6.csv', range6)
writeCSV(basefolder + 'Range7.csv', range6)
writeCSV(basefolder + 'Range8.csv', range8)
writeCSV(basefolder + 'Range9.csv', range8)

# write "dummy" files
dummyData = [[1000000, 1, 0], [1, 1, 0]]
writeCSV(basefolder + 'V_inputBuffer.csv', dummyData)
writeCSV(basefolder + 'WEgain1_Igain1.csv', dummyData)

print("Script complete, press enter to close this window...")
input()

#V1_I1 = deepcopy(H_BaselineV)
#PolDivideCol_2_3(V1_I1, H_BaselineV)


# data testing section*************************************************************************
#rawSquidstatData = parseRawCSV(basefolder + '/squidstatRawData2.csv') # re-import data
#InvertPhase(rawSquidstatData)
#results = []

#V1 = V1_I1
#V2 = V2_I1
#V5 = V5_I1
#V10 = V10_I1
#V20 = V20_I10
#V50 = V50_I10
#V100 = V100_I10
#I1 = V1_I1
#I2 = I2_V1
#I5 = I5_V1
#I10 = I10_V1
#I20 = I20_V10
#I50 = I50_V10
#I100 = I100_V10

##group 1
#VgainList = [V1, V2, V5, V10, V1, V1, V1]
#IgainList = [I1, I1, I1, I1, I2, I5, I10]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

##group 2
#VgainList = [V10, V20, V50, V100, V10, V10, V10]
#IgainList = [I10, I10, I10, I10, I20, I50, I100]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

##group 3
#VgainList = [V10, V20, V50, V100, V50, V50, V50]
#IgainList = [I50, I50, I50, I50, I20, I50, I100]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

##group 4
#VgainList = [V10, V20, V50, V100, V20, V20, V20]
#IgainList = [I20, I20, I20, I20, I20, I50, I100]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

##group 5
#VgainList = [V10, V20, V50, V100, V100, V100, V100]
#IgainList = [I100, I100, I100, I100, I20, I50, I100]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

##group 6
#VgainList = [V10, V20, V50, V100, V1, V1, V1]
#IgainList = [I1, I1, I1, I1, I20, I50, I100]
#for i in range(0, len(VgainList)):
#    data = popSweep(rawSquidstatData)
#    PolDivideCol_2_3(data, H_BaselineV)
#    PolDivideCol_2_3(data, VgainList[i])
#    PolMultCol_2_3(data, IgainList[i])
#    results += data

#writeCSV(basefolder + '/testresults.csv', results)