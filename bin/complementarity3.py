#!/usr/bin/env python3


from OligoCalc import Input
from OligoCalc import Start
from OligoCalc import MakeComplement
from sys import argv


debugHairpin = 0
debugDimers = 0
doTiming = 0
theAlignedArray = 0
theHairpinArray = 0
minAlignLen = 0
minHairpinLen = 0
maxMismatchNum = 1
# hairpins must have this many bases between self-annealed sequences
bubbleSize = 3


def add(a, b):
    a_is_none = a is None
    b_is_none = b is None
    if a_is_none or b_is_none:
        return None
    return a + b


def op_monad(op):
    def apply(a, b):
        a_is_none = a is None
        b_is_none = b is None
        if a_is_none or b_is_none:
            return None
        return op(a, b)
    return apply


def makeMatrix(matLength):
    return [[None for _ in range(matLength)] for _ in range(matLength)]
    # var theMatrix = new Array(matLength)
    # for (var i=0
    #      i < matLength
    #      i++) {
    #     // increment column
    #     theMatrix[i] = new Array(matLength)
    # }
    # return theMatrix


def fillMatchMatrix(cols, rows, mat):
    d = len(cols)
    # print(rows)
    if (d < 4):
        return
    global broadMatch

    if (broadMatch):
        # // Do the degenerate thing!
        for i in range(d):
            # // increment column
            for j in range(d):
                # // increment row
                if (isBaseEqual(cols[i], rows[j])):
                    mat[i][j] = 1
                    if (i > 0 and j > 0):
                        mat[i][j] = add(mat[i][j], mat[i - 1][j - 1])
                        # mat[i][j] += mat[i - 1][j - 1]
                        # // (increment diagonal values)
                else:
                    mat[i][j] = 0
                    if (i > 1 and j > 1):
                        if (mat[i - 1][j - 1] > mat[i - 2][j - 2] and mat[i - 1][j - 1] > 1 and i < d - 1 and j < d - 1):
                            # // allow one base mismatch only if there are at least 2 matched base on 5' and at least 1 matched base on 3'
                            mat[i][j] = mat[i - 1][j - 1]
                        elif (i < d - 1 and j < d - 1):
                            mat[i - 1][j - 1] = 0
    else:
        for i in range(2):
            for j in range(2):
                # // increment column
                # // increment row
                if (cols[i] == rows[j]):
                    mat[i][j] = 1
                    if (i and j):
                        # mat[i][j] += mat[i - 1][j - 1]
                        mat[i][j] = add(mat[i][j], mat[i - 1][j - 1])
                    # // (increment diagonal values)
                else:
                    mat[i][j] = 0
        for i in range(2, d - 1):
            # // increment column
            for j in range(2, d - 1):
                # // increment row
                if (cols[i] == rows[j]):
                    # mat[i][j] = mat[i - 1][j - 1] + 1
                    mat[i][j] = add(mat[i - 1][j - 1], 1)
                    # // (increment diagonal values)
                else:
                    mat[i][j] = 0
                    if (op_monad(lambda a, b: a > b)(mat[i - 1][j - 1], 1) and cols[i + 1] == rows[j + 1]):
                        # // allow one base mismatch only if there are at least 2 matched base on 5' and at least 1 matched base on 3'
                        mat[i][j] = mat[i - 1][j - 1]
        i = d - 1
        j = i
        # // increment column
        # // increment row
        if (cols[i] == rows[j]):
            mat[i][j] = 1
            # mat[i][j] += mat[i - 1][j - 1]
            mat[i][j] = add(mat[i][j], mat[i - 1][j - 1])
            # // (increment diagonal values)
        else:
            mat[i][j] = 0


def makeAlignedArray(mat, minLen, maxMisMatch):
    # // assumes an orthogonal matrix
    # /* theAlignedArray is a bit strange in the second dimension. Assume it is a length 5 array called 'theResults'
    # theResults[0] == start index
    # theResults[1] == start matching index in reverse complement seq
    # theResults[2] == end index of aligned bases (inclusive)
    # theResults[3] == end matching index in reverse complement Seq
    # theResults[4] == number of mismatches
    # */

    matLength = len(mat)

    # count = 0
    theResults = []
    i = None
    j = None
    k = None
    mismatches = None

    # s = ""
    for i in range(matLength):
        for j in range(matLength):
            if (mat[i][j] == 1):  # //potential start of an alignment
                # print(f'{i} {j}')
                mismatches = 0
                hasMatch = 1
                lastMatch = 1
                maxInc = matLength - (j if i <= j else i)
                # s = f'{s} {maxInc}'

                k = 1
                while k < maxInc:
                    hasMatch = mat[i + k][j + k]

                    # s = f'{s} {hasMatch}'
                    if (not hasMatch):
                        # print(hasMatch)
                        break
                    # s = f'{s} {hasMatch}'
                    if (hasMatch <= lastMatch):
                        if (mismatches >= maxMisMatch):
                            break
                        mismatches += 1
                    lastMatch = hasMatch
                    k += 1

                # k = maxInc
                # print(maxMisMatch)
                # s = f'{s} {k}'
                # print(mismatches)

                if (k - mismatches >= minLen):
                    # if count == len(theResults):
                    #     theResults.append(None)
                    # theResults[-1] = [0, 0, 0, 0, 0]
                    # print(i, j, k, maxInc)
                    theResults.append([0, 0, 0, 0, 0])
                    theResults[-1][0] = i  # ;    //start index
                    # ;    //start matching index in reverse complement seq
                    theResults[-1][1] = j
                    # ; //end index of aligned bases (inclusive)
                    theResults[-1][2] = i + k - 1
                    # ; //end matching index in reverse complement Seq
                    theResults[-1][3] = j + k - 1
                    theResults[-1][4] = mismatches  # ;  //mismatch counts
                    # count += 1
    # print(s)
    # for line in theResults:
        # print(line)
    return theResults


def sortAlignedArray(alignedArray):
    # // assumes an orthogonal matrix
    # / * theAlignedArray is a bit strange in the second dimension. Assume it is a length 5 array called 'theResults'
    # theResults[0] == start index
    # theResults[1] == start matching index in reverse complement seq
    # theResults[2] == end index of aligned bases(inclusive)
    # theResults[3] == end matching index in reverse complement Seq
    # theResults[4] == number of mismatches
    # * /
    if (len(alignedArray) > 2):
        if (1 == 2):
            # print("Che")
            tempArray = [0, 0, 0, 0, 0]
            # swapped = 0
            # run_once = True
            # // bubble sort
            while True:
                swapped = 0
                for n in range(len(alignedArray) - 2):
                    if (alignedArray[n][2] - alignedArray[n][0] < alignedArray[n + 1][2] - alignedArray[n + 1][0]):
                        for i in range(5):
                            tempArray[i] = alignedArray[n][i]
                            alignedArray[n][i] = alignedArray[n + 1][i]
                            alignedArray[n + 1][i] = tempArray[i]
                        swapped = 1
                if swapped != 1:
                    break
        else:
            # print("Guevara")
            alignedArray.sort(key=arrayOrder())
            # alignedArray = sorted(alignedArray, key=arrayOrder())  # MYMOD
    # for line in alignedArray:
        # print(line)
    return alignedArray


def comparator(mycmp):
    class K:
        def __init__(self, obj, *args):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0

    return K


def arrayOrder():
    def mycmp(a, b):
        # //size plus position
        # print(a, b)
        # print('----------------------')
        v = (1 if ((a[2] - a[0]) < (b[2] - b[0])) else (-1 if ((a[2] -
                                                                a[0]) > (b[2] - b[0])) else (a[0] - a[1]) - (b[0] - b[1])))
        return v

    return comparator(mycmp)


def calculateMatrices(theOligo, theComplement):
    # var theStart = datetime.now()
    if (len(theOligo.Sequence) != len(theComplement.Sequence)):
        raise ValueError(
            "Error! Primer and its complement are different lengths!")

    # setup d*d matrix
    matrix = makeMatrix(len(theOligo.Sequence))

    fillMatchMatrix(theOligo.seqArray, theComplement.seqArray, matrix)

    global theAlignedArray
    global maxMismatchNum
    global minAlignLen
    theAlignedArray = makeAlignedArray(
        matrix, minAlignLen, maxMismatchNum)

    theAlignedArray = sortAlignedArray(theAlignedArray)


def display3EndDimer(theOligo, theAlignedArray):
    d = len(theOligo.Sequence)
    returnString = ""
    N = 0

    # 3' complementarity
    returnString += "3' Complementarity: "
    # print()
    for n in range(len(theAlignedArray) - 1):
        # print(theAlignedArray[n])
        if(theAlignedArray[n][2] == d-1):  # end position of match in original seq
            N += 1
    # print()
    return N


def reverseString(string):
    return string[::-1]


def stringToArray(string):
    return list(string)


def calcPrimer(form):
    # theStart = None
    # theEnd = None
    # calcStart = datetime.now()

    if (len(form.oligoBox) < 8):
        raise ValueError(
            "Please enter at least 8 bases before checking for self-complementarity!")

    theOligo, theComplement = Start(form)
    # print(theOligo.Sequence)
    # print(theComplement.Sequence)

    if (theOligo.seqArray):
        del theOligo.seqArray
    if (theOligo.revSeqArray):
        del theOligo.revSeqArray
    if (theComplement.seqArray):
        del theComplement.seqArray
    if (theComplement.revSeqArray):
        del theComplement.revSeqArray

    theOligo.revSequence = reverseString(theOligo.Sequence)
    theComplement.revSequence = reverseString(theComplement.Sequence)

    theOligo.seqArray = stringToArray(theOligo.Sequence)
    theOligo.revSeqArray = stringToArray(theOligo.revSequence)
    theComplement.seqArray = stringToArray(theComplement.Sequence)
    theComplement.revSeqArray = stringToArray(theComplement.revSequence)

    # // change if removing the hairpin selection options

    global minAlignLen
    global minHairpinLen

    minAlignLen = int(form.selfComp)  # // for 3' complementarity
    minHairpinLen = int(form.hairpin)  # // for hairpin

    global broadMatch
    broadMatch = False

    calculateMatrices(theOligo, theComplement)

    global theAlignedArray
    r = display3EndDimer(theOligo, theAlignedArray)
    # print(r)

    return r


if __name__ == "__main__":
    sequence = argv[1]
    form = Input(oligoBox=sequence)
    if len(argv) == 4:
        form.selfComp = int(argv[2])
        form.hairpin = int(argv[3])
    else:
        form.selfComp = 5
        form.hairpin = 4
    print(calcPrimer(form))
