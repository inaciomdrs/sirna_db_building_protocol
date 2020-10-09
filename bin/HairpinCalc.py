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


def reverseString(string):
    return string[::-1]


def stringToArray(string):
    return list(string)


def isBaseEqual(c1, c2):
    if (c1 == c2):
        return True
    global broadMatch
    if (broadMatch):
        if (c1 == 'N' or c2 == 'N'):
            return True

        equA = "AMRWVHD"
        # // lack of 'M' caught by Paul Wayper. Thanks Paul!
        equT = "TWYKHDB"
        equG = "GRSKVDB"
        equC = "CMSYVHB"
        # // lack of 'M' caught by Paul Wayper. Thanks Paul!

        if (c1 == 'A'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                return False
        if (c1 == 'T'):
            try:
                equT.index(c2)
                return True
            except ValueError:
                return False
        if (c1 == 'G'):
            try:
                equG.index(c2)
                return True
            except ValueError:
                return False
        if (c1 == 'C'):
            try:
                equC.index(c2)
                return True
            except ValueError:
                return False
        if (c1 == 'M'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equC.index(c2)
                    return True
                except ValueError:
                    return False
        if (c1 == 'R'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equG.index(c2)
                    return True
                except ValueError:
                    return False
        if (c1 == 'W'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equT.index(c2)
                    return True
                except ValueError:
                    return False
        if (c1 == 'S'):
            try:
                equG.index(c2)
                return True
            except ValueError:
                try:
                    equC.index(c2)
                    return True
                except ValueError:
                    return False
        if (c1 == 'Y'):
            try:
                equT.index(c2)
                return True
            except ValueError:
                try:
                    equC.index(c2)
                    return True
                except ValueError:
                    return False
        if (c1 == 'K'):
            try:
                equT.index(c2)
                return True
            except ValueError:
                try:
                    equG.index(c2)
                    return True
                except ValueError:
                    return False

        if (c1 == 'V'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equG.index(c2)
                    return True
                except ValueError:
                    try:
                        equC.index(c2)
                        return True
                    except ValueError:
                        return False
        if (c1 == 'H'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equT.index(c2)
                    return True
                except ValueError:
                    try:
                        equC.index(c2)
                        return True
                    except ValueError:
                        return False
        if (c1 == 'D'):
            try:
                equA.index(c2)
                return True
            except ValueError:
                try:
                    equT.index(c2)
                    return True
                except ValueError:
                    try:
                        equG.index(c2)
                        return True
                    except ValueError:
                        return False
        if (c1 == 'B'):
            try:
                equT.index(c2)
                return True
            except ValueError:
                try:
                    equG.index(c2)
                    return True
                except ValueError:
                    try:
                        equC.index(c2)
                        return True
                    except ValueError:
                        return False
    return False


def getIndexOf(seq, subSeq, startIndex, minMatch):
    # // look for subSeq in seq
    # / * returns an array where
    # theResult[0] is the index of the first match of subseq that is of at least length minMatch in seq
    # theResult[1] is the length of the match
    # * /
    theResult = [-1, -1]
    # theResult[0] = -1
    # theResult[1] = -1
    global broadMatch
    if (not broadMatch):
        for k in range(minMatch, len(subSeq)+1):
            # // can replace this with seq.search for GREP capabilities
            try:
                theMatch = seq.index(subSeq[0:k], startIndex)
            except ValueError:
                break
            theResult[0] = theMatch
            theResult[1] = k
            # if (debugHairpin) primerWin.document.write("(" + theMatch + "," + k + ") ")
    else:
        for i in range(startIndex, len(seq)):
            if (isBaseEqual(seq[i], subSeq[0])):
                for j in range(len(subSeq)):
                    if (not isBaseEqual(seq[i + j], subSeq[j])):
                        break
                    elif (j >= minMatch - 1):
                        theResult[0] = theMatch
                        theResult[1] = k
                if (j == subSeq.length):
                    theResult[0] = theMatch
                    theResult[1] = k
    # if (debugHairpin) primerWin.document.write("TheResult[0]=" + theResult[0] + " (first match); TheResult[1]=" + theResult[1] + ";<br>")
    return theResult


def DoHairpinArrayInsert(a, b, c, d, results):
    arrayCount = len(results)
    if (a >= c or a >= b or c >= d or b >= c):
        # if (debugHairpin) primerWin.document.write("DoHairpinArrayInsert: ERROR IN VALUES PASSED! [0]=" + a + "; [1]=" + b + "[2]=" + c + "; [3]=" + d + ";<br>\n")
        # print('DoHairpinArrayInsert - Return on 1st if')
        return results

    for i in range(arrayCount):
        if (results[i][0] <= a and results[i][1] >= b and results[i][2] <= c and results[i][3] >= d):
            # print('DoHairpinArrayInsert - Return on 2nd if')
            return results
        if (results[i][0] >= a and results[i][1] <= b and results[i][2] >= c and results[i][3] <= d):
            # print('DoHairpinArrayInsert - Return on 3rd if')
            results[i][0] = a
            results[i][1] = b
            results[i][2] = c
            results[i][3] = d
            # if (debugHairpin) primerWin.document.write("DoHairpinArrayInsert: position " + i + " in results replaced with [0]=" + a + "; [1]=" + b + "[2]=" + c + "; [3]=" + d + ";<br>")
            return results

    # results[arrayCount] = [0, 0, 0, 0]
    results.append([0, 0, 0, 0])  # MYMOD
    results[arrayCount][0] = a
    results[arrayCount][1] = b
    results[arrayCount][2] = c
    results[arrayCount][3] = d
    # print('DoHairpinArrayInsert - Return on final if')
    # if (debugHairpin) primerWin.document.write("DoHairpinArrayInsert: arrayCount=" + arrayCount + "; [0]=" + a + "; [1]=" + b + "[2]=" + c + "; [3]=" + d + ";<br>")
    return results


def calcHairpin(theFullSequence, minHairpinLength):
    # /* compare theCompSeq with theFullSeq starting at theFullSeq[startPos]. Successful matches must be at least minMatch long * /
    # /* The resulting array is an array of arrays. each result should be an array of 4 integers
    # result[0]: position of start of match in sequence
    # result[1]: position of end of match
    # result[2]: position of the start of the complement(really the end since it would be 3'-5')
    # result[3]: position of the end of the complement(really the start since it would be 3'-5')
    # */

    theFullComplement = MakeComplement(theFullSequence, True)
    theResults = []

    # if (debugHairpin) primerWin.document.write("<PRE>")
    # if (debugHairpin) primerWin.document.write("calcHairpin: theFullSequence  =" + theFullSequence + "; theFullSequence.length=" + theFullSequence.length + "; minHairpinLen" + minHairpinLen + ";<br>")
    # if (debugHairpin) primerWin.document.write("calcHairpin: theFullComplement=" + theFullComplement + "; theFullComplement.length=" + theFullComplement.length + ";<br>")

    theResult = None
    count = None
    compPos = None
    seqPos = None

    global bubbleSize
    # // makes sure that we do not anneal the full length of the primer - that should come out in the dimerization report
    maxSeqLength = abs(len(theFullSequence) / 2) - bubbleSize
    # print(f'maxSeqLength = {maxSeqLength}')
    maxMatch = 0

    # console.log(maxSeqLength)
    compPos = 0
    seqPos = None
    while compPos < len(theFullComplement) - 2 * minHairpinLength:
        # print(f' <><> seqPos = {seqPos} | compPos = {compPos}')
        # for compPos in range(len(theFullComplement) - 2 * minHairpinLength):
        maxMatch = 0
        seqPos = 0
        while seqPos < len(theFullSequence) - maxSeqLength:
            # for seqPos in range(len(theFullSequence) - maxSeqLength):
            # print(f'seqPos = {seqPos} | compPos = {compPos}')
            # if (debugHairpin) primerWin.document.write("calcHairpin: compPos=" + compPos + "; seqPos=" + seqPos + ";<br>")
            theResult = getIndexOf(
                theFullSequence[0:seqPos + int(maxSeqLength)],
                theFullComplement[compPos:len(theFullComplement)],
                seqPos, minHairpinLength)
            # for line in theResult:
                # print(line, end=' ')
            # print()
            if (theResult[0] > -1):
                # // theResult[0] is the index of the first match of theFullComplement that is of at least length minHairpinLength in theFullSequence
                # // theResult[1] is the length of the match

                theResults = DoHairpinArrayInsert(
                    theResult[0],
                    theResult[0] + theResult[1] - 1,
                    len(theFullSequence) - compPos - theResult[1],
                    len(theFullSequence) - compPos - 1,
                    theResults
                )
                # print(f'theResult = {theResult}')
                # print(f'maxMatch = {maxMatch}')
                if (theResult[1] > maxMatch):
                    maxMatch = theResult[1]
                # ; // move forward to guarantee nothing else is found that is a reasonable match
                seqPos = theResult[0] + theResult[1] - minHairpinLength

                # print(f">{seqPos}, {minHairpinLength}, {maxSeqLength}<")
                if (seqPos + minHairpinLength >= maxSeqLength):
                    # ; // move compPos forward to stop identical checks if long match was found!
                    compPos += maxMatch - minHairpinLength
                    # print(f' <> seqPos = {seqPos} | compPos = {compPos}')
                    break  # ; // we have moved far enough on the primer to guarentee we have everything - further would give us the reverse match
            else:
                if (maxMatch > minHairpinLength):
                    # ; // move compPos forward to stop identical checks if long match was found!
                    compPos += maxMatch - minHairpinLength
                break  # ; // not found in the rest of the sequence!
            seqPos += 1
        compPos += 1
    # if (debugHairpin) primerWin.document.write("<\/PRE>")
    return theResults


def calcPrimer(form):
    if (len(form.oligoBox) < 8):
        raise ValueError(
            "Please enter at least 8 bases before checking for self-complementarity!")

    theOligo, theComplement = Start(form)

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

    theHairpinArray = calcHairpin(theOligo.Sequence, minHairpinLen)
    r = len(theHairpinArray)
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
