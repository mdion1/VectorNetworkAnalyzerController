import sys
import csv
from copy import deepcopy
import numpy as np
#from scipy.signal import savgol_filter

#****************** Helper functions**********************
#def smooth_col_2_3(datatable):
#    y1data = []
#    y2data = []
#    for i in range(0,len(datatable)):
#        y1data.append(datatable[i][1])
#        y2data.append(datatable[i][2])
#    y1data = savgol_filter(y1data, 5, 2) # window size 51, polynomial order 3
#    y2data = savgol_filter(y2data, 5, 2) # window size 51, polynomial order 3
#    for i in range(0,len(datatable)):
#        datatable[i][1] = y1data[i]
#        datatable[i][2] = y2data[i]

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

def calculate(blockMeas, roughMeas, baseline, outputPwrColumn):
    normalize(blockMeas)
    PolDivideCol_2_3(blockMeas, baseline)

def calculateAlt1(blockMeas, baseline):
    normalize(blockMeas, lowestMagFrequency = 9999)
    PolDivideCol_2_3(blockMeas, baseline)


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
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Frequency', 'Magnitude', 'Phase(degrees)', 'Output power (dB)'])
        for i in range(0, len(data)):
            csv_writer.writerow(data[i])
    csv_file.close()

def PolDivideCol_2_3(numerator, denominator):
    for i in range(0, len(denominator)):
        numerator[i][1] /= denominator[i][1]
        numerator[i][2] -= denominator[i][2]

def PolInvCol_2_3(baseTable):
    ret = deepcopy(baseTable)
    for i in range(0, len(baseTable)):
        ret[i][1] = 1 / baseTable[i][1]
        ret[i][2] = -baseTable[i][2]
    return ret

def PolMultCol_2_3(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        baseTable[i][1] *= auxTable[i][1]
        baseTable[i][2] += auxTable[i][2]

def Append4thCol(baseTable, auxTable):
    for i in range(0, len(baseTable)):
        baseTable[i].append(auxTable[i][0])


#*************** Main script ****************************

basefolder = sys.argv[1]
alternate = True #alternate = False

# import "lumped" measurements (all of the signal blocks measured together)
H_Baseline = parseRawCSV(basefolder + '/BaselineRaw.csv')
normalize(H_Baseline)
control = parseRawCSV(basefolder + '/control.csv')
control_10DB_A = parseRawCSV(basefolder + '/control_attenA_10DB.csv')
control_10DB_B = parseRawCSV(basefolder + '/control_attenB_10DB.csv')
normalize(control)
normalize(control_10DB_A)
normalize(control_10DB_B)
PolDivideCol_2_3(H_Baseline, control)
H_BaselineInv = PolInvCol_2_3(H_Baseline)

H_TIA_1_1 = deepcopy(H_BaselineInv)
H_TIA_1_2 = PolInvCol_2_3(parseRawCSV(basefolder + '/TIA_stage1setting2raw.csv'))
#H_TIA_1_3 = parseRawCSV(basefolder + '/TIA_stage1setting3raw.csv')
#H_TIA_1_4 = parseRawCSV(basefolder + '/TIA_stage1setting4raw.csv')
H_TIA_2_1 = deepcopy(H_BaselineInv)
H_TIA_2_2 = PolInvCol_2_3(parseRawCSV(basefolder + '/TIA_stage2setting2raw.csv'))
#H_IBuf = PolInvCol_2_3(parseRawCSV(basefolder + '/I_inputBufferRaw.csv'))
H_Igain1 = deepcopy(H_BaselineInv)
H_Igain2 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain2raw.csv'))
H_Igain5 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain5raw.csv'))
H_Igain10 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain10raw.csv'))
H_Igain20 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain20raw.csv'))
H_Igain50 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain50raw.csv'))
H_Igain100 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain100raw.csv'))
H_Igain200 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain200raw.csv'))
H_Igain500 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain500raw.csv'))
H_Igain1000 = PolInvCol_2_3(parseRawCSV(basefolder + '/Igain1000raw.csv'))
PolMultCol_2_3(H_TIA_1_2, control)
PolMultCol_2_3(H_TIA_2_2, control)
#PolMultCol_2_3(H_IBuf, control)
PolMultCol_2_3(H_Igain2, control)
PolMultCol_2_3(H_Igain5, control)
PolMultCol_2_3(H_Igain10, control)
PolMultCol_2_3(H_Igain20, control)
PolMultCol_2_3(H_Igain50, control)
PolMultCol_2_3(H_Igain100, control)
PolMultCol_2_3(H_Igain200, control)
PolMultCol_2_3(H_Igain500, control)
PolMultCol_2_3(H_Igain1000, control_10DB_B)

H_VBuf = deepcopy(H_Baseline)
H_Vgain1 = deepcopy(H_Baseline)
H_Vgain2 = parseRawCSV(basefolder + '/WEgain2raw.csv')
H_Vgain5 = parseRawCSV(basefolder + '/WEgain5raw.csv')
H_Vgain10 = parseRawCSV(basefolder + '/WEgain10raw.csv')
H_Vgain20 = parseRawCSV(basefolder + '/WEgain20raw.csv')
H_Vgain50 = parseRawCSV(basefolder + '/WEgain50raw.csv')
H_Vgain100 = parseRawCSV(basefolder + '/WEgain100raw.csv')
H_Vgain200 = parseRawCSV(basefolder + '/WEgain200raw.csv')
H_Vgain500 = parseRawCSV(basefolder + '/WEgain500raw.csv')
H_Vgain1000 = parseRawCSV(basefolder + '/WEgain1000raw.csv')
PolDivideCol_2_3(H_Vgain2, control)
PolDivideCol_2_3(H_Vgain5, control)
PolDivideCol_2_3(H_Vgain10, control)
PolDivideCol_2_3(H_Vgain20, control)
PolDivideCol_2_3(H_Vgain50, control)
PolDivideCol_2_3(H_Vgain100, control)
PolDivideCol_2_3(H_Vgain200, control)
PolDivideCol_2_3(H_Vgain500, control)
PolDivideCol_2_3(H_Vgain1000, control_10DB_A)

# calculate the transfer function of each block
calculate(H_TIA_1_1, H_BaselineInv)
calculateAlt1(H_TIA_1_2, H_BaselineInv)
#calculate(H_TIA_1_3, H_BaselineInv)
#calculate(H_TIA_1_4, H_BaselineInv)
calculate(H_TIA_2_1, H_BaselineInv)
calculateAlt1(H_TIA_2_2, H_BaselineInv)
#calculate(H_IBuf, H_BaselineInv)
calculate(H_Igain1, H_BaselineInv)
calculate(H_Igain2, H_BaselineInv)
calculate(H_Igain5, H_BaselineInv)
calculate(H_Igain10, H_BaselineInv)
calculate(H_Igain20, H_BaselineInv)
calculate(H_Igain50, H_BaselineInv)
calculate(H_Igain100, H_BaselineInv)
calculate(H_Igain200, H_BaselineInv)
calculate(H_Igain500, H_BaselineInv)
calculate(H_Igain1000, H_BaselineInv)
calculate(H_VBuf, H_Baseline)
calculate(H_Vgain1, H_Baseline)
calculate(H_Vgain2, H_Baseline)
calculate(H_Vgain5, H_Baseline)
calculate(H_Vgain10, H_Baseline)
calculate(H_Vgain20, H_Baseline)
calculate(H_Vgain50, H_Baseline)
calculate(H_Vgain100, H_Baseline)
calculate(H_Vgain200, H_Baseline)
calculate(H_Vgain500, H_Baseline)
calculate(H_Vgain1000, H_Baseline)

# calculate the "residual" transfer function
H_residual = deepcopy(H_Baseline)  

# write the resulting data to CSV files
writeCSV(basefolder + '/Residual.csv', H_residual)
writeCSV(basefolder + '/TIA_stage1setting1.csv', H_TIA_1_1)
writeCSV(basefolder + '/TIA_stage1setting2.csv', H_TIA_1_2)
#writeCSV(basefolder + '/TIA_stage1setting3.csv', H_TIA_1_3)
#writeCSV(basefolder + '/TIA_stage1setting4.csv', H_TIA_1_4)
writeCSV(basefolder + '/TIA_stage2setting1.csv', H_TIA_2_1)
writeCSV(basefolder + '/TIA_stage2setting2.csv', H_TIA_2_2)
#writeCSV(basefolder + '/I_inputBuffer.csv', H_IBuf)
writeCSV(basefolder + '/Igain1.csv', H_Igain1)
writeCSV(basefolder + '/Igain2.csv', H_Igain2)
writeCSV(basefolder + '/Igain5.csv', H_Igain5)
writeCSV(basefolder + '/Igain10.csv', H_Igain10)
writeCSV(basefolder + '/Igain20.csv', H_Igain20)
writeCSV(basefolder + '/Igain50.csv', H_Igain50)
writeCSV(basefolder + '/Igain100.csv', H_Igain100)
writeCSV(basefolder + '/Igain200.csv', H_Igain200)
writeCSV(basefolder + '/Igain500.csv', H_Igain500)
writeCSV(basefolder + '/Igain1000.csv', H_Igain1000)
writeCSV(basefolder + '/V_inputBuffer.csv', H_VBuf)
writeCSV(basefolder + '/WEgain1.csv', H_Vgain1)
writeCSV(basefolder + '/WEgain2.csv', H_Vgain2)
writeCSV(basefolder + '/WEgain5.csv', H_Vgain5)
writeCSV(basefolder + '/WEgain10.csv', H_Vgain10)
writeCSV(basefolder + '/WEgain20.csv', H_Vgain20)
writeCSV(basefolder + '/WEgain50.csv', H_Vgain50)
writeCSV(basefolder + '/WEgain100.csv', H_Vgain100)
writeCSV(basefolder + '/WEgain200.csv', H_Vgain200)
writeCSV(basefolder + '/WEgain500.csv', H_Vgain500)
writeCSV(basefolder + '/WEgain1000.csv', H_Vgain1000)