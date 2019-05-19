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

def normalize(table, lowestMagFrequency = 0):
    BasePhase = table[-1][2]
    BaseMagIndex = findIndex(table, lowestMagFrequency)
    BaseMag = table[BaseMagIndex][1]
    for i in range(0,len(table)):
        if i <= BaseMagIndex:
            table[i][1] /= BaseMag
        else:
            table[i][1] = 1
        table[i][2] -= BasePhase

def findIndex(table, minFreq):
    index = 0
    for i in range(0,len(table)):
        if table[i][0] > minFreq:
            index = i
            continue
        else:
            break
    return index


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

def PolDivideCol_2_3(numerator, denominator):
    for i in range(0, len(numerator)):
        if len(denominator) > i:
            numerator[i][1] /= denominator[i][1]
            numerator[i][2] -= denominator[i][2]

def PolInvCol_2_3(baseTable):
    ret = deepcopy(baseTable)
    for i in range(0, len(baseTable)):
        ret[i][1] = 1 / baseTable[i][1]
        ret[i][2] = -baseTable[i][2]
    return ret

def InvertPhase(baseTable):
    for i in range(0, len(baseTable)):
        baseTable[i][2] *= -1

def PolMultCol_2_3(baseTable, auxTable):
    length = min(len(baseTable), len(auxTable))
    for i in range(0, length):
        baseTable[i][1] *= auxTable[i][1]
        baseTable[i][2] += auxTable[i][2]

def PolAddCol_2_3(baseTable, auxTable):
    length = min(len(baseTable), len(auxTable))
    for i in range(0, length):
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = auxTable[i][2] / 180 * math.pi
        aux = complex(auxTable[i][1] * math.cos(base_phi), auxTable[i][1] * math.sin(base_phi))
        sum = base + aux
        baseTable[i][1] = abs(sum)
        baseTable[i][2] = cmath.phase(sum) * 180 / math.pi

def PolSubtractCol_2_3(baseTable, auxTable):
    length = min(len(baseTable), len(auxTable))
    for i in range(0, length):
        base_phi = baseTable[i][2] / 180 * math.pi
        base = complex(baseTable[i][1] * math.cos(base_phi), baseTable[i][1] * math.sin(base_phi))
        aux_phi = auxTable[i][2] / 180 * math.pi
        aux = complex(auxTable[i][1] * math.cos(base_phi), auxTable[i][1] * math.sin(base_phi))
        diff = base - aux
        baseTable[i][1] = abs(diff)
        baseTable[i][2] = cmath.phase(diff) * 180 / math.pi

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
        if len(row) > 3:
            ret.append([freq, mag / N, phase / N, IsigMag / N])
        else:
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

rawDataRange1 = parseRawCSV(sys.argv[5])
InvertPhase(rawDataRange1)

rawDataRange3 = parseRawCSV(sys.argv[6])
InvertPhase(rawDataRange3)

rawDataRange4 = parseRawCSV(sys.argv[7])
InvertPhase(rawDataRange4)

rawDataRange5 = parseRawCSV(sys.argv[8])
InvertPhase(rawDataRange5)

Range2_Rerr = parseRawCSV(sys.argv[9])
Range3_Rerr = parseRawCSV(sys.argv[10])
Range4_Rerr = parseRawCSV(sys.argv[11])
Range1_Rerr = parseRawCSV(sys.argv[12])
Range5_Rerr = parseRawCSV(sys.argv[13])

#debugging only
#Range2_Rerr = [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]
#Range3_Rerr = Range2_Rerr
#Range4_Rerr = Range2_Rerr
#Range1_Rerr = Range2_Rerr
#Range5_Rerr = Range2_Rerr

# parse and normalize sweeps from each gain setting
range2 = popSweep(rawDataTier1)
H_BaselineV = deepcopy(range2)
H_BaselineI = PolInvCol_2_3(H_BaselineV)
V2 = popSweep(rawDataTier1)
V5 = popSweep(rawDataTier1)
V10 = popSweep(rawDataTier1)
I2 = PolInvCol_2_3(popSweep(rawDataTier1))
I5 = PolInvCol_2_3(popSweep(rawDataTier1))
I10 = PolInvCol_2_3(popSweep(rawDataTier1))
V20 = popSweep(rawDataTier2)
V50 = popSweep(rawDataTier2)
V100 = popSweep(rawDataTier2)
I20 = PolInvCol_2_3(popSweep(rawDataTier2))
I50 = PolInvCol_2_3(popSweep(rawDataTier2))
I100 = PolInvCol_2_3(popSweep(rawDataTier2))
V200 = popSweep(rawDataTier3)
V500 = popSweep(rawDataTier3)
V1000 = popSweep(rawDataTier3)
I200 = PolInvCol_2_3(popSweep(rawDataTier3))
I500 = PolInvCol_2_3(popSweep(rawDataTier3))
I1000 = PolInvCol_2_3(popSweep(rawDataTier3))
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
PolDivideCol_2_3(V2, H_BaselineV)
PolDivideCol_2_3(V5, H_BaselineV)
PolDivideCol_2_3(V10, H_BaselineV)
PolDivideCol_2_3(I2, H_BaselineI)
PolDivideCol_2_3(I5, H_BaselineI)
PolDivideCol_2_3(I10, H_BaselineI)

#calculate tier2 gains (x20 thru x100)
PolDivideCol_2_3(V20, H_BaselineV)
PolMultCol_2_3(V20, I10)
PolDivideCol_2_3(V50, H_BaselineV)
PolMultCol_2_3(V50, I10)
PolDivideCol_2_3(V100, H_BaselineV)
PolMultCol_2_3(V100, I10)
PolDivideCol_2_3(I20, H_BaselineI)
PolMultCol_2_3(I20, V10)
PolDivideCol_2_3(I50, H_BaselineI)
PolMultCol_2_3(I50, V10)
PolDivideCol_2_3(I100, H_BaselineI)
PolMultCol_2_3(I100, V10)

#calculate tier2 gains (x200 thru x1000)
PolDivideCol_2_3(V200, H_BaselineV)
PolMultCol_2_3(V200, I100)
PolDivideCol_2_3(V500, H_BaselineV)
PolMultCol_2_3(V500, I100)
PolDivideCol_2_3(V1000, H_BaselineV)
PolMultCol_2_3(V1000, I100)
PolDivideCol_2_3(I200, H_BaselineI)
PolMultCol_2_3(I200, V100)
PolDivideCol_2_3(I500, H_BaselineI)
PolMultCol_2_3(I500, V100)
PolDivideCol_2_3(I1000, H_BaselineI)
PolMultCol_2_3(I1000, V100)

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
range1 = popSweep(rawDataRange1)
range3 = popSweep(rawDataRange3)
range4ctrl = popSweep(rawDataRange4)
range4 = popSweep(rawDataRange4)
range5ctrl = popSweep(rawDataRange5)
range5 = popSweep(rawDataRange5)
ScalarMult(Range2_Rerr, 100 / 10)
ScalarMult(Range3_Rerr, 1000 / 10)
ScalarMult(Range4_Rerr, 10000 / 10)
ScalarMult(Range5_Rerr, 100000 / 10)
#todo: is it best to smooth before or after corrections?
smooth_col_2_3(Range1_Rerr, re_normalize = False)
smooth_col_2_3(Range2_Rerr, re_normalize = False)
smooth_col_2_3(Range3_Rerr, re_normalize = False)
smooth_col_2_3(Range4_Rerr, re_normalize = False)
smooth_col_2_3(Range5_Rerr, re_normalize = False)

#debugging: recalculate Range5_Rerr based on Range4_Rerr, not on VNA measurement
r5Norm = deepcopy(range5); PolDivideCol_2_3(r5Norm, range5ctrl); normalize(r5Norm)
Range5_Rerr = deepcopy(Range4_Rerr); PolMultCol_2_3(Range5_Rerr, r5Norm); smooth_col_2_3(Range5_Rerr, re_normalize = False)

#range0/1
PolDivideCol_2_3(range1, V10); PolMultCol_2_3(range1, I100)     #divide out Igain = 100 and Vgain = 10
#todo: PolAddCol_2_3(range1, Range1_Rerr)                              #account for R_err
normalize(range1)
range1 = PolInvCol_2_3(range1)
range0 = deepcopy(range1)

#range2
PolSubtractCol_2_3(range2, Range2_Rerr)     #PolAddCol_2_3(range2, Range2_Rerr)                              #account for R_err
normalize(range2)
range2 = PolInvCol_2_3(range2)

#range3
PolDivideCol_2_3(range3, V10);                                                                                #divide out Igain = 10 and Vgain = 10
PolSubtractCol_2_3(range3, Range3_Rerr)      #PolAddCol_2_3(range3, Range3_Rerr)                              #account for R_err
normalize(range3)
range3 = PolInvCol_2_3(range3)

##range4
#PolDivideCol_2_3(range4ctrl, V10); PolMultCol_2_3(range4ctrl, I10); PolAddCol_2_3(range4ctrl, Range2_Rerr)  #divide out Igain = 10 and Vgain = 10; account for R_err
#PolDivideCol_2_3(range4, V10); PolAddCol_2_3(range4, Range4_Rerr)                                           #divide out Vgain = 10; account for R_err
#PolDivideCol_2_3(range4, range4ctrl); PolMultCol_2_3(range4, PolInvCol_2_3(range3))                         #normalize to a load of 100Ohm
#normalize(range4)
#range4 = PolInvCol_2_3(range4)

#alternate range4
PolDivideCol_2_3(range4ctrl, V10); PolMultCol_2_3(range4ctrl, I10); PolSubtractCol_2_3(range4ctrl, Range3_Rerr)  #divide out Igain = 10 and Vgain = 10; account for R_err
PolDivideCol_2_3(range4, V10); PolSubtractCol_2_3(range4, Range4_Rerr)                                           #divide out Vgain = 10; account for R_err
load_adj_1k_to_100R = deepcopy(range4ctrl); PolDivideCol_2_3(load_adj_1k_to_100R, PolInvCol_2_3(range3))    #calculate adjustment from 1kOhm to 100Ohm
smooth_col_2_3(load_adj_1k_to_100R)                                                                         #remove noise from control runs
PolDivideCol_2_3(range4, load_adj_1k_to_100R)
normalize(range4)
range4 = PolInvCol_2_3(range4)

##range5
#PolDivideCol_2_3(range5ctrl, V10); PolAddCol_2_3(range5ctrl, Range4_Rerr)  #divide out Igain = 10 and Vgain = 10; account for R_err
#PolDivideCol_2_3(range5, V10); PolAddCol_2_3(range5, Range5_Rerr)                                           #divide out Vgain = 10; account for R_err
#PolDivideCol_2_3(range5, range5ctrl); PolMultCol_2_3(range5, PolInvCol_2_3(range4))                         #normalize to a load of 100Ohm
#normalize(range5)
#range5 = PolInvCol_2_3(range5)

#alternate range5
#PolDivideCol_2_3(range5ctrl, V10); PolSubtractCol_2_3(range5ctrl, Range4_Rerr)  #divide out Igain = 10 and Vgain = 10; account for R_err
PolDivideCol_2_3(range5, V10);  PolSubtractCol_2_3(range5, Range5_Rerr)                                     #divide out Vgain = 10; account for R_err
#load_adj_10k_to_100R = deepcopy(range5ctrl); PolDivideCol_2_3(load_adj_10k_to_100R, PolInvCol_2_3(range4))  #calculate adjustment from 1kOhm to 100Ohm
#smooth_col_2_3(load_adj_10k_to_100R)                                                                        #remove noise from control runs
PolDivideCol_2_3(range5, load_adj_1k_to_100R) #PolDivideCol_2_3(range5, load_adj_10k_to_100R)
normalize(range5)
range5 = PolInvCol_2_3(range5)

# calculate R_err tables
PolDivideCol_2_3(Range1_Rerr, PolInvCol_2_3(range1))
PolDivideCol_2_3(Range2_Rerr, PolInvCol_2_3(range2))
PolDivideCol_2_3(Range3_Rerr, PolInvCol_2_3(range3))
PolDivideCol_2_3(Range4_Rerr, PolInvCol_2_3(range4))
PolDivideCol_2_3(Range5_Rerr, PolInvCol_2_3(range5))
#todo: is it better to smooth here, after corrections are performed?
#smooth_col_2_3(Range2_Rerr, re_normalize = False)
#smooth_col_2_3(Range4_Rerr, re_normalize = False)

# write the resulting data to CSV files
writeCSV(basefolder + 'Igain2.csv', I2)
writeCSV(basefolder + 'Igain5.csv', I5)
writeCSV(basefolder + 'Igain10.csv', I10)
writeCSV(basefolder + 'Igain20.csv', I20)
writeCSV(basefolder + 'Igain50.csv', I50)
writeCSV(basefolder + 'Igain100.csv', I100)
writeCSV(basefolder + 'Igain200.csv', I200)
writeCSV(basefolder + 'Igain500.csv', I500)
writeCSV(basefolder + 'Igain1000.csv', I1000)
writeCSV(basefolder + 'WEgain2.csv', V2)
writeCSV(basefolder + 'WEgain5.csv', V5)
writeCSV(basefolder + 'WEgain10.csv', V10)
writeCSV(basefolder + 'WEgain20.csv', V20)
writeCSV(basefolder + 'WEgain50.csv', V50)
writeCSV(basefolder + 'WEgain100.csv', V100)
writeCSV(basefolder + 'WEgain200.csv', V200)
writeCSV(basefolder + 'WEgain500.csv', V500)
writeCSV(basefolder + 'WEgain1000.csv', V1000)
writeCSV(basefolder + 'Range0.csv', range0)
writeCSV(basefolder + 'Range1.csv', range1)
writeCSV(basefolder + 'Range2.csv', range2)
writeCSV(basefolder + 'Range3.csv', range3)
writeCSV(basefolder + 'Range4.csv', range4)
writeCSV(basefolder + 'Range5.csv', range5)
writeCSV(basefolder + 'Range1_Rerr.csv', Range1_Rerr)
writeCSV(basefolder + 'Range2_Rerr.csv', Range2_Rerr)
writeCSV(basefolder + 'Range3_Rerr.csv', Range3_Rerr)
writeCSV(basefolder + 'Range4_Rerr.csv', Range4_Rerr)
writeCSV(basefolder + 'Range5_Rerr.csv', Range5_Rerr)

# write "dummy" files
dummyData = [[1000000, 1, 0], [1, 1, 0]]
writeCSV(basefolder + 'V_inputBuffer.csv', dummyData)
writeCSV(basefolder + 'Igain1.csv', dummyData)
writeCSV(basefolder + 'WEgain1.csv', dummyData)

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