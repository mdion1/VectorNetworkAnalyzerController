import sys
import csv
from copy import deepcopy
import numpy as np
from scipy.signal import savgol_filter

#****************** Helper functions**********************
def smooth_col_2_3(datatable):
    y1data = []
    y2data = []
    for i in range(0,len(datatable)):
        y1data.append(datatable[i][1])
        y2data.append(datatable[i][2])
    y1data = savgol_filter(y1data, 5, 2) # window size 51, polynomial order 3
    y2data = savgol_filter(y2data, 5, 2) # window size 51, polynomial order 3
    #for i in range(0,len(datatable)):
    #    datatable[i][1] = y1data[i]
    #    datatable[i][2] = y2data[i]

def calculate(blockMeas, roughMeas, baseline, outputPwrColumn):
    if alternate:
        PolDivideCol_2_3(blockMeas, baseline)
        smooth_col_2_3(blockMeas)
        #PolMultCol_2_3(blockMeas, roughMeas)
        Append4thCol(blockMeas, outputPwrColumn)
    else:
        PolDivideCol_2_3(blockMeas, baseline)
        smooth_col_2_3(blockMeas)
        PolMultCol_2_3(blockMeas, roughMeas)
        Append4thCol(blockMeas, outputPwrColumn)


def parseRawCSV(filename):
    dataTable = []
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
    return dataTable

def writeCSV(filename, data):
    with open(filename, mode = 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')
        csv_writer.writerow(['Frequency', 'Magnitude', 'Phase(degrees)', 'Output power (dB)'])
        for i in range(0, len(data)):
            csv_writer.writerow(data[i])
    csv_file.close()

def PolDivideCol_2_3(numerator, denominator):
    for i in range(0, len(numerator)):
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
control = parseRawCSV(basefolder + '/control.csv')
PolDivideCol_2_3(H_Baseline, control)
H_BaselineInv = PolInvCol_2_3(H_Baseline)

H_TIA_1_1 = deepcopy(H_BaselineInv)
H_TIA_1_2 = PolInvCol_2_3(parseRawCSV(basefolder + '/TIA_stage1setting2raw.csv'))
#H_TIA_1_3 = parseRawCSV(basefolder + '/TIA_stage1setting3raw.csv')
#H_TIA_1_4 = parseRawCSV(basefolder + '/TIA_stage1setting4raw.csv')
H_TIA_2_1 = deepcopy(H_BaselineInv)
H_TIA_2_2 = PolInvCol_2_3(parseRawCSV(basefolder + '/TIA_stage2setting2raw.csv'))
#H_IBuf = PolInvCol_2_3(parseRawCSV(basefolder + '/I_inputBufferRaw.csv'))
H_Igain1_1 = deepcopy(H_BaselineInv)
H_Igain1_2 = PolInvCol_2_3(parseRawCSV(basefolder + '/IgainBlock1_gain2raw.csv'))
H_Igain1_5 = PolInvCol_2_3(parseRawCSV(basefolder + '/IgainBlock1_gain5raw.csv'))
H_Igain1_10 = PolInvCol_2_3(parseRawCSV(basefolder + '/IgainBlock1_gain10raw.csv'))
H_Igain2_1 = deepcopy(H_BaselineInv)
H_Igain2_10 = PolInvCol_2_3(parseRawCSV(basefolder + '/IgainBlock2_gain10_alt1raw.csv'))
H_Igain3_1 = deepcopy(H_BaselineInv)
H_Igain3_10 = PolInvCol_2_3(parseRawCSV(basefolder + '/IgainBlock3_gain10_alt2raw.csv'))
PolMultCol_2_3(H_TIA_1_2, control)
PolMultCol_2_3(H_TIA_2_2, control)
#PolMultCol_2_3(H_IBuf, control)
PolMultCol_2_3(H_Igain1_2, control)
PolMultCol_2_3(H_Igain1_5, control)
PolMultCol_2_3(H_Igain1_10, control)
PolMultCol_2_3(H_Igain2_10, control)
PolMultCol_2_3(H_Igain3_10, control)


H_VBuf = deepcopy(H_Baseline)
H_Vgain1_1 = deepcopy(H_Baseline)
H_Vgain1_2 = parseRawCSV(basefolder + '/WEgainBlock1_gain2raw.csv')
H_Vgain1_5 = parseRawCSV(basefolder + '/WEgainBlock1_gain5raw.csv')
H_Vgain1_10 = parseRawCSV(basefolder + '/WEgainBlock1_gain10raw.csv')
H_Vgain2_1 = deepcopy(H_Baseline)
H_Vgain2_10 = parseRawCSV(basefolder + '/WEgainBlock2_gain10_alt1raw.csv')
H_Vgain3_1 = deepcopy(H_Baseline)
H_Vgain3_10 = parseRawCSV(basefolder + '/WEgainBlock3_gain10_alt2raw.csv')
PolDivideCol_2_3(H_Vgain1_2, control)
PolDivideCol_2_3(H_Vgain1_5, control)
PolDivideCol_2_3(H_Vgain1_10, control)
PolDivideCol_2_3(H_Vgain2_10, control)
PolDivideCol_2_3(H_Vgain3_10, control)

# import individual, "rough" measurements (each of the signal blocks measured one at a time)
h_TIA_1_1 = parseRawCSV(basefolder + '/TIA_stage1setting1rough.csv')
h_TIA_2_1 = parseRawCSV(basefolder + '/TIA_stage2setting1rough.csv')
PolDivideCol_2_3(h_TIA_1_1, control)
PolDivideCol_2_3(h_TIA_2_1, control)
h_TIA_1_2_combine = deepcopy(h_TIA_1_1)
PolMultCol_2_3(h_TIA_1_2_combine, h_TIA_2_1)
h_Igain1_1 = parseRawCSV(basefolder + '/IgainBlock1_gain1rough.csv')
h_Igain2_1 = parseRawCSV(basefolder + '/IgainBlock2_gain1rough.csv')
h_Igain3_1 = parseRawCSV(basefolder + '/IgainBlock3_gain1rough.csv')
PolDivideCol_2_3(h_Igain1_1, control)
PolDivideCol_2_3(h_Igain2_1, control)
PolDivideCol_2_3(h_Igain3_1, control)

h_VBuf = parseRawCSV(basefolder + '/V_inputBufferRough.csv')
h_Vgain1_1 = parseRawCSV(basefolder + '/WEgainBlock1_gain1rough.csv')
h_Vgain2_1 = parseRawCSV(basefolder + '/WEgainBlock2_gain1rough.csv')
h_Vgain3_1 = parseRawCSV(basefolder + '/WEgainBlock3_gain1rough.csv')
PolDivideCol_2_3(h_VBuf, control)
PolDivideCol_2_3(h_Vgain1_1, control)
PolDivideCol_2_3(h_Vgain2_1, control)
PolDivideCol_2_3(h_Vgain3_1, control)

#import output power columns
outputPwr_TIA_1_1 = parseRawCSV(basefolder + "/TIA_stage1setting1_baseline.csv")
outputPwr_TIA_1_2 = parseRawCSV(basefolder + "/TIA_stage1setting2_baseline.csv")
#outputPwr_TIA_1_3 = parseRawCSV(basefolder + "/TIA_stage1setting3_baseline.csv")
#outputPwr_TIA_1_4 = parseRawCSV(basefolder + "/TIA_stage1setting4_baseline.csv")
outputPwr_TIA_2_1 = parseRawCSV(basefolder + "/TIA_stage2setting1_baseline.csv")
outputPwr_TIA_2_2 = parseRawCSV(basefolder + "/TIA_stage2setting2_baseline.csv")
outputPwr_Igain1_1 = parseRawCSV(basefolder + "/IgainBlock1_gain1_baseline.csv")
outputPwr_Igain1_2 = parseRawCSV(basefolder + "/IgainBlock1_gain2_baseline.csv")
outputPwr_Igain1_5 = parseRawCSV(basefolder + "/IgainBlock1_gain5_baseline.csv")
outputPwr_Igain1_10 = parseRawCSV(basefolder + "/IgainBlock1_gain10_baseline.csv")
outputPwr_Igain2_1 = parseRawCSV(basefolder + "/IgainBlock2_gain1_baseline.csv")
outputPwr_Igain2_10 = parseRawCSV(basefolder + "/IgainBlock2_gain10_baseline.csv")
outputPwr_Igain3_1 = parseRawCSV(basefolder + "/IgainBlock3_gain1_baseline.csv")
outputPwr_Igain3_10 = parseRawCSV(basefolder + "/IgainBlock3_gain10_baseline.csv")
#outputPwr_IBuf = parseRawCSV(basefolder + "/I_inputBuffer_baseline.csv")
outputPwr_VBuf = parseRawCSV(basefolder + "/V_inputBuffer_baseline.csv")
outputPwr_Vgain1_1 = parseRawCSV(basefolder + "/WEgainBlock1_gain1_baseline.csv")
outputPwr_Vgain1_2 = parseRawCSV(basefolder + "/WEgainBlock1_gain2_baseline.csv")
outputPwr_Vgain1_5 = parseRawCSV(basefolder + "/WEgainBlock1_gain5_baseline.csv")
outputPwr_Vgain1_10 = parseRawCSV(basefolder + "/WEgainBlock1_gain10_baseline.csv")
outputPwr_Vgain2_1 = parseRawCSV(basefolder + "/WEgainBlock2_gain1_baseline.csv")
outputPwr_Vgain2_10 = parseRawCSV(basefolder + "/WEgainBlock2_gain10_baseline.csv")
outputPwr_Vgain3_1 = parseRawCSV(basefolder + "/WEgainBlock3_gain1_baseline.csv")
outputPwr_Vgain3_10 = parseRawCSV(basefolder + "/WEgainBlock3_gain10_baseline.csv")

# calculate the transfer function of each block
calculate(H_TIA_1_1, h_TIA_1_1, H_BaselineInv, outputPwr_TIA_1_1)
calculate(H_TIA_1_2, h_TIA_1_1, H_BaselineInv, outputPwr_TIA_1_2)
#calculate(H_TIA_1_3, h_TIA_1_1, H_BaselineInv, outputPwr_TIA_1_3)
#calculate(H_TIA_1_4, h_TIA_1_1, H_BaselineInv, outputPwr_TIA_1_4)
calculate(H_TIA_2_1, h_TIA_2_1, H_BaselineInv, outputPwr_TIA_2_1)
calculate(H_TIA_2_2, h_TIA_2_1, H_BaselineInv, outputPwr_TIA_2_2)
#calculate(H_IBuf, h_TIA_1_2_combine, H_BaselineInv, outputPwr_IBuf)
calculate(H_Igain1_1, h_Igain1_1, H_BaselineInv, outputPwr_Igain1_1)
calculate(H_Igain1_2, h_Igain1_1, H_BaselineInv, outputPwr_Igain1_2)
calculate(H_Igain1_5, h_Igain1_1, H_BaselineInv, outputPwr_Igain1_5)
calculate(H_Igain1_10, h_Igain1_1, H_BaselineInv, outputPwr_Igain1_10)
calculate(H_Igain2_1, h_Igain2_1, H_BaselineInv, outputPwr_Igain2_1)
calculate(H_Igain2_10, h_Igain2_1, H_BaselineInv, outputPwr_Igain2_10)
calculate(H_Igain3_1, h_Igain3_1, H_BaselineInv, outputPwr_Igain3_1)
calculate(H_Igain3_10, h_Igain3_1, H_BaselineInv, outputPwr_Igain3_10)
calculate(H_VBuf, h_VBuf, H_Baseline, outputPwr_VBuf)
calculate(H_Vgain1_1, h_Vgain1_1, H_Baseline, outputPwr_Vgain1_1)
calculate(H_Vgain1_2, h_Vgain1_1, H_Baseline, outputPwr_Vgain1_2)
calculate(H_Vgain1_5, h_Vgain1_1, H_Baseline, outputPwr_Vgain1_5)
calculate(H_Vgain1_10, h_Vgain1_1, H_Baseline, outputPwr_Vgain1_10)
calculate(H_Vgain2_1, h_Vgain2_1, H_Baseline, outputPwr_Vgain2_1)
calculate(H_Vgain2_10, h_Vgain2_1, H_Baseline, outputPwr_Vgain2_10)
calculate(H_Vgain3_1, h_Vgain3_1, H_Baseline, outputPwr_Vgain3_1)
calculate(H_Vgain3_10, h_Vgain3_1, H_Baseline, outputPwr_Vgain3_10)

# calculate the "residual" transfer function
H_residual = deepcopy(H_Baseline)
if not alternate:
    PolDivideCol_2_3(H_residual, h_VBuf)
    PolDivideCol_2_3(H_residual, h_Igain1_1)
    PolDivideCol_2_3(H_residual, h_Igain2_1)
    PolDivideCol_2_3(H_residual, h_Igain3_1)
    PolMultCol_2_3(H_residual, h_TIA_1_1)
    PolMultCol_2_3(H_residual, h_TIA_2_1)
    PolMultCol_2_3(H_residual, h_Igain1_1)
    PolMultCol_2_3(H_residual, h_Igain2_1)
    PolMultCol_2_3(H_residual, h_Igain3_1)
    

# write the resulting data to CSV files
writeCSV(basefolder + '/Residual.csv', H_residual)
writeCSV(basefolder + '/TIA_stage1setting1.csv', H_TIA_1_1)
writeCSV(basefolder + '/TIA_stage1setting2.csv', H_TIA_1_2)
#writeCSV(basefolder + '/TIA_stage1setting3.csv', H_TIA_1_3)
#writeCSV(basefolder + '/TIA_stage1setting4.csv', H_TIA_1_4)
writeCSV(basefolder + '/TIA_stage2setting1.csv', H_TIA_2_1)
writeCSV(basefolder + '/TIA_stage2setting2.csv', H_TIA_2_2)
#writeCSV(basefolder + '/I_inputBuffer.csv', H_IBuf)
writeCSV(basefolder + '/IgainBlock1_gain1.csv', H_Igain1_1)
writeCSV(basefolder + '/IgainBlock1_gain2.csv', H_Igain1_2)
writeCSV(basefolder + '/IgainBlock1_gain5.csv', H_Igain1_5)
writeCSV(basefolder + '/IgainBlock1_gain10.csv', H_Igain1_10)
writeCSV(basefolder + '/IgainBlock2_gain1.csv', H_Igain2_1)
writeCSV(basefolder + '/IgainBlock2_gain10_alt1.csv', H_Igain2_10)
writeCSV(basefolder + '/IgainBlock3_gain1.csv', H_Igain3_1)
writeCSV(basefolder + '/IgainBlock3_gain10_alt2.csv', H_Igain3_10)
writeCSV(basefolder + '/V_inputBuffer.csv', H_VBuf)
writeCSV(basefolder + '/WEgainBlock1_gain1.csv', H_Vgain1_1)
writeCSV(basefolder + '/WEgainBlock1_gain2.csv', H_Vgain1_2)
writeCSV(basefolder + '/WEgainBlock1_gain5.csv', H_Vgain1_5)
writeCSV(basefolder + '/WEgainBlock1_gain10.csv', H_Vgain1_10)
writeCSV(basefolder + '/WEgainBlock2_gain1.csv', H_Vgain2_1)
writeCSV(basefolder + '/WEgainBlock2_gain10_alt1.csv', H_Vgain2_10)
writeCSV(basefolder + '/WEgainBlock3_gain1.csv', H_Vgain3_1)
writeCSV(basefolder + '/WEgainBlock3_gain10_alt2.csv', H_Vgain3_10)
