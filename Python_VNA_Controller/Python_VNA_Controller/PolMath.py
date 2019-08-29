from copy import deepcopy
import math
import cmath

def PolDivide(numerator, denominator):
    denom_copy = deepcopy(denominator)
    removeOffsetAt10kHz(denom_copy)
    for i in range(0, len(numerator)):
        freq = numerator[i][0]
        numerator[i][1] /= interpolate(freq, denom_copy, 1)
        numerator[i][2] -= interpolate(freq, denom_copy, 2)

def PolInv(baseTable):
    ret = deepcopy(baseTable)
    for i in range(0, len(baseTable)):
        ret[i][1] = 1 / baseTable[i][1]
        ret[i][2] = -baseTable[i][2]
    return ret

def PolMult(baseTable, auxTable):
    aux_copy = deepcopy(auxTable)
    removeOffsetAt10kHz(aux_copy)
    for i in range(0, len(baseTable)):
        freq = baseTable[i][0]
        baseTable[i][1] *= interpolate(freq, aux_copy, 1)
        baseTable[i][2] += interpolate(freq, aux_copy, 2)

def find10kHzIndex(sweep):
    for i in range(0, len(sweep)):
        if sweep[i][0] > 10000:
            continue
        else:
            return i - 1

def removeOffsetAt10kHz(sweep):
    index = find10kHzIndex(sweep)
    if index == -1:
        return

    initialOffset = sweep[index][1] / sweep[index + 1][1]
    for i in range(0, index + 1):
        sweep[i][1] /= initialOffset

def getV1_I1_offset(range2sweep):
    index = find10kHzIndex(range2sweep)
    ratio = range2sweep[index][1] / range2sweep[index + 1][1]
    freq1 = range2sweep[index][0]
    freq2 = range2sweep[index + 1][0]
    ret = [[freq1, ratio, 0], [freq2, 1, 0]]
    return ret

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