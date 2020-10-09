#!/usr/bin/env python3
from OligoCalc import Input
from OligoCalc import Start
from OligoCalc import MakeComplement
from sys import argv
# from datetime import datetime


# GLOBALS

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


def calcDegeneratePrimers(theOligo, theComplement):
    if (not theOligo.hasIUpacBase):
        return

    global broadMatch
    broadMatch = True  # True: do all degenerate comparisons
    calculateMatrices(theOligo, theComplement)
    # anchorString = "<font COLOR='green'>-----------------------------<BR>-----------------------------<BR><\/FONT>"
    # doc.write(anchorString.anchor("allMatches"))

    # doc.write("Your oligo contains degenerated bases.<BR>")
    # doc.write("This section displays <font COLOR='RED'>all potential matches <\/FONT>in the case of degenerated bases.<BR>")
    # doc.write("For Example:<BR> <PRE>   'N' matches 'A','T','G','C', or 'N';<BR>   'R' matches 'T','C','S','N';<BR>   'W' matches 'A','T','W','N'; etc.<\/PRE><P>")
    # hrefString = "view strict matches only"
    # doc.writeln(
    #     "Scroll up to view <FONT COLOR='green'> <a href='#strictMatches'>strict Matches<\/a> <\/FONT>")

    # if (!isCompatible) {
    #     doc.write("<p><b>Sorry, the hairpin loop calculation is only available if you are using IE or Netscape 4.x or higher!!\n<\/B><br>")
    # } else {
    #     doc.writeln(displayHairpin(theHairpinArray, theOligo.Sequence))
    # }
    print(display3EndDimer(theOligo, theAlignedArray))
    print(displayAllDimers(theAlignedArray, theOligo.Sequence, theOligo.revSequence))


def calculateMatrices(theOligo, theComplement):
    # var theStart = datetime.now()
    if (len(theOligo.Sequence) != len(theComplement.Sequence)):
        raise ValueError(
            "Error! Primer and its complement are different lengths!")

    # setup d*d matrix
    matrix = makeMatrix(len(theOligo.Sequence))

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - makeMatrix took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()

    # //populates the matrix

    fillMatchMatrix(theOligo.seqArray, theComplement.seqArray, matrix)

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - fillMatchMatrix took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()

    # if (isIE)
    # if (theAlignedArray) delete theAlignedArray

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - delete theAlignedArray took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    # //exam the matrix for 3 prime complementary

    global theAlignedArray
    global maxMismatchNum
    global minAlignLen
    theAlignedArray = makeAlignedArray(matrix, minAlignLen, maxMismatchNum)

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - makeAlignedArray took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    # delete matrix
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - delete matrix took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    # theAlignedArray = sortAlignedArray(theAlignedArray)
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - sortAlignedArray took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    # //exam the sequence for potential hairpins
    # if (isCompatible) {
    #     if (isIE)
    #         if (theHairpinArray) delete theHairpinArray
    #     theHairpinArray = calcHairpin(theOligo.Sequence, minHairpinLen)
    # }
    global theHairpinArray
    global minHairpinLen
    theHairpinArray = calcHairpin(theOligo.Sequence, minHairpinLen)

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calculateMatrices - calcHairpin took " + (theEnd - theStart) + " ms<br>")


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
    # count = None
    compPos = None
    seqPos = None

    # // makes sure that we do not anneal the full length of the primer - that should come out in the dimerization report
    maxSeqLength = abs(len(theFullSequence) // 2) - bubbleSize
    maxMatch = 0

    # console.log(maxSeqLength)

    for compPos in range(len(theFullComplement) - 2 * minHairpinLength):
        maxMatch = 0
        for seqPos in range(len(theFullSequence) - maxSeqLength):
            # if (debugHairpin) primerWin.document.write("calcHairpin: compPos=" + compPos + "; seqPos=" + seqPos + ";<br>")
            theResult = getIndexOf(
                theFullSequence[0:seqPos + maxSeqLength],
                theFullComplement[compPos:len(theFullComplement)],
                seqPos, minHairpinLength)
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
                if (theResult[1] > maxMatch):
                    maxMatch = theResult[1]
                # ; // move forward to guarantee nothing else is found that is a reasonable match
                seqPos = theResult[0] + theResult[1] - minHairpinLength
                if (seqPos + minHairpinLength >= maxSeqLength):
                    # ; // move compPos forward to stop identical checks if long match was found!
                    compPos += maxMatch - minHairpinLength
                    break  # ; // we have moved far enough on the primer to guarentee we have everything - further would give us the reverse match
            else:
                if (maxMatch > minHairpinLength):
                    # ; // move compPos forward to stop identical checks if long match was found!
                    compPos += maxMatch - minHairpinLength
                break  # ; // not found in the rest of the sequence!
    # if (debugHairpin) primerWin.document.write("<\/PRE>")
    return theResults


def DoHairpinArrayInsert(a, b, c, d, results):
    arrayCount = len(results)
    if (a >= c or a >= b or c >= d or b >= c):
        # if (debugHairpin) primerWin.document.write("DoHairpinArrayInsert: ERROR IN VALUES PASSED! [0]=" + a + "; [1]=" + b + "[2]=" + c + "; [3]=" + d + ";<br>\n")
        return results

    for i in range(arrayCount):
        if (results[i][0] <= a and results[i][1] >= b and results[i][2] <= c and results[i][3] >= d):
            return results
        if (results[i][0] >= a and results[i][1] <= b and results[i][2] >= c and results[i][3] <= d):
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
    # if (debugHairpin) primerWin.document.write("DoHairpinArrayInsert: arrayCount=" + arrayCount + "; [0]=" + a + "; [1]=" + b + "[2]=" + c + "; [3]=" + d + ";<br>")
    return results


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


def makeMatrix(matLength):
    return [[0 for _ in range(matLength)] for _ in range(matLength)]
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
                        mat[i][j] += mat[i - 1][j - 1]
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
                        mat[i][j] += mat[i - 1][j - 1]
                    # // (increment diagonal values)
                else:
                    mat[i][j] = 0
        for i in range(2, d - 1):
            # // increment column
            for j in range(2, d - 1):
                # // increment row
                if (cols[i] == rows[j]):
                    mat[i][j] = mat[i - 1][j - 1] + 1
                    # // (increment diagonal values)
                else:
                    mat[i][j] = 0
                    if (mat[i - 1][j - 1] > 1 and cols[i + 1] == rows[j + 1]):
                        # // allow one base mismatch only if there are at least 2 matched base on 5' and at least 1 matched base on 3'
                        mat[i][j] = mat[i - 1][j - 1]
        i = d - 1
        j = i
        # // increment column
        # // increment row
        if (cols[i] == rows[j]):
            mat[i][j] = 1
            mat[i][j] += mat[i - 1][j - 1]
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
    count = 0
    theResults = []
    i = None
    j = None
    k = None
    mismatches = None
    for i in range(matLength):
        for j in range(matLength):
            if (mat[i][j] == 1):  # //potential start of an alignment
                mismatches = 0
                hasMatch = 1
                lastMatch = 1
                maxInc = matLength - (j if i <= j else i)
                for k in range(1, maxInc):
                    hasMatch = mat[i + k][j + k]
                    if (not hasMatch):
                        break
                    if (hasMatch <= lastMatch):
                        if (mismatches >= maxMisMatch):
                            break
                        mismatches += 1
                    lastMatch = hasMatch

                if (k - mismatches >= minLen):
                    if count == len(theResults):
                        theResults.append(None)
                    theResults[count] = [0, 0, 0, 0, 0]
                    theResults[count][0] = i  # ;    //start index
                    # ;    //start matching index in reverse complement seq
                    theResults[count][1] = j
                    # ; //end index of aligned bases (inclusive)
                    theResults[count][2] = i + k - 1
                    # ; //end matching index in reverse complement Seq
                    theResults[count][3] = j + k - 1
                    theResults[count][4] = mismatches  # ;  //mismatch counts
                    count += 1

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
            tempArray = [0, 0, 0, 0, 0]
            swapped = 0
            run_once = True
            # // bubble sort
            while (swapped == 1) or run_once:
                run_once = False
                swapped = 0
                for n in range(len(alignedArray) - 2):
                    if (alignedArray[n][2] - alignedArray[n][0] < alignedArray[n + 1][2] - alignedArray[n + 1][0]):
                        for i in range(5):
                            tempArray[i] = alignedArray[n][i]
                            alignedArray[n][i] = alignedArray[n + 1][i]
                            alignedArray[n + 1][i] = tempArray[i]
                        swapped = 1
        else:
            alignedArray = sorted(alignedArray, key=arrayOrder())  # MYMOD

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
        return (1 if ((a[2] - a[0]) < (b[2] - b[0])) else (-1 if ((a[2] - a[0]) > (b[2] - b[0])) else (a[0] - a[1]) - (b[0] - b[1])))

    return comparator(mycmp)


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


# // return string containing all hairpins
def displayHairpin(theHairpinArray, theSequence):
    returnString = ""
    # s1 = ""
    # d = len(theSequence)
    # i, j = 0, 0
    theHairpinArrayLength = len(theHairpinArray)
    # Potential hairpin formation: "
    returnString = returnString + "Number of harpins: "
    if (theHairpinArrayLength > 0):
        if (theHairpinArrayLength > 1):
            if (theHairpinArray[theHairpinArrayLength - 1][1] == theHairpinArray[theHairpinArrayLength - 2][1] and (theHairpinArray[theHairpinArrayLength - 1][2] == theHairpinArray[theHairpinArrayLength - 2][2])):
                theHairpinArrayLength = theHairpinArrayLength - 1
                # // get rid of the last one
        # for i in range(theHairpinArrayLength):
        #     # // add a bar between 2 legs of the hairpin if bases in the 2nd leg is contiguous to the 1st leg
        #     # // substring wants a value from the start location to 1+the end location
        #     s1 = theSequence.substring(0, theHairpinArray[i][0]) +
        #     theSequence.substring(theHairpinArray[i][0], theHairpinArray[i][1] + 1).fontcolor("red") +
        #     ((theHairpinArray[i][1] + 1 >= theHairpinArray[i][2]) ? "-": "") +
        #     theSequence.substring(theHairpinArray[i][1] + 1, theHairpinArray[i][2]) +
        #     theSequence.substring(theHairpinArray[i][2], theHairpinArray[i][3] + 1).fontcolor("red") +
        #     theSequence.substring(theHairpinArray[i][3] + 1, d)
        #     returnString = returnString + "5' " + s1 + " 3'<BR>"

        returnString += str(theHairpinArrayLength)
        # returnString = returnString + "</PRE>"
    else:
        returnString += " 0"

    return returnString


def display3EndDimer(theOligo, theAlignedArray):
    d = len(theOligo.Sequence)
    returnString = ""
    # s1 = ""
    # s2 = ""
    # // 3' complementarity
    returnString += "3' Complementarity: "
    N = 0

    for n in range(len(theAlignedArray)):
        # // end position of match in original seq
        if (theAlignedArray[n][2] == d - 1):
            N += 1

    returnString += f"{N}"
    return returnString


# // all possible dimerization sites
def displayAllDimers(theAlignedArray, theSequence, reversedSeq):
    # / * theAlignedArray is a bit strange in the second dimension. Assume it is a length 5 array called 'theResults'
    # theResults[0] == start index
    # theResults[1] == start matching index in reverse complement seq
    # theResults[2] == end index of aligned bases(inclusive)
    # theResults[3] == end matching index in reverse complement Seq
    # theResults[4] == number of mismatches
    # * /
    # d = len(theSequence)
    returnString = ""
    # s1 = ""
    # s2 = ""
    # maxoffset = 0
    # offset, j, n = None, None, None
    # offsetStr, maxoffsetStr = None, None
    # // all other possible alignment sites
    returnString += "All potential self-annealing sites are marked in red (allowing 1 mis-match): "
    if (len(theAlignedArray) > 1):
        returnString += str(len(theAlignedArray))
    else:
        returnString += "0"

    # returnString += "\n"
    return returnString


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

    # // minAlignLen = parseInt(form.selfComp.value)
    # // minHairpinLen = parseInt(form.hairpin.value)
    global broadMatch
    broadMatch = False
    # // True: do all degenerate comparisons
    # // now create window
    # primerWin = window.open(
    # "", 'primer', 'width=700,toolbar=0,location=0,directories=0,status=1,menuBar=1,scrollBars=1,resizable=1')
    # // primerWin.document.open("text/html")
    # print(
    # "<HTML><HEAD><TITLE>Oligo Self Complementarity Check<\/TITLE><\/HEAD>")
    # primerWin.document.write('<STYLE type="text\/css">')
    # primerWin.document.write('<!--')
    # if (isMac and browserVersion < 5.0) {
    # primerWin.document.write(MacStyleSheet)
    # } else {
    # primerWin.document.write(PCStyleSheet)
    # }
    # primerWin.document.write('-->')
    # primerWin.document.write('<\/STYLE>')
    # print("<BODY BGCOLOR=white>")
    # print(
    # "Minimum base pairs required for single primer self-dimerization: " + minAlignLen + ".<BR>")
    # print(
    # "Minimum base pairs required for a hairpin: " + minHairpinLen + ".<BR>")
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calcPrimer - initialize window took " + (theEnd - calcStart) + " ms<br>")
    calculateMatrices(theOligo, theComplement)
    # // do this after the window is in place so we can print diagnostics if necessary
    # anchorString = ""
    # print(anchorString.anchor("strictMatches"))
    # if (theOligo.hasIUpacBase) {
    #     print("Your oligo contains degenerate bases.<BR>")
    #     print(
    #         "The strictMatch section displays only <font COLOR='RED'>perfect matches <\/FONT>in the case of degenerate bases.<BR>")
    #     print(
    #         "For Example:<BR> <PRE>   'N' matches only with 'N';<BR>   'R' matches only with 'S';<BR>   'W' matches only with 'W'; etc.<\/PRE><P>")
    #     hrefString = "view all matches."
    #     print(
    #         "Scroll down to view <font COLOR='GREEN'><a href='#allMatches'>all Matches<\/a><\/FONT>")
    # }
    print(displayHairpin(theHairpinArray, theOligo.Sequence), end='\t')

    # if (!isCompatible) {
    # print(
    # "<p><b>Sorry, the hairpin loop calculation is only available if you are using Netscape or IE 4.0 or higher!\n<\/B><br>")
    # } else {
    # theStart = datetime.now()
    # print(displayHairpin(
    # theHairpinArray, theOligo.Sequence))
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calcPrimer - displayHairpin took " + (theEnd - theStart) + " ms<br>")
    # }
    # theStart = datetime.now()

    print(display3EndDimer(theOligo, theAlignedArray), end='\t')

    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calcPrimer - display3EndDimer took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    print(displayAllDimers(theAlignedArray, theOligo.Sequence, theOligo.revSequence))
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calcPrimer - displayAllDimers took " + (theEnd - theStart) + " ms<br>")
    # theStart = datetime.now()
    if (theOligo.hasIUpacBase):
        calcDegeneratePrimers(theOligo, theComplement)
    # theEnd = datetime.now()
    # if (doTiming) primerWin.document.write("calcPrimer - all calls to close took " + (theEnd - calcStart) + " ms<br>")
    # print("<\/BODY><\/HTML>")
    # primerWin.document.close()
    # primerWin.focus()
    # return False


if __name__ == "__main__":
    sequence = argv[1]
    form = Input(oligoBox=sequence)
    if len(argv) == 3:
        form.selfComp = int(argv[2])
        form.hairpin = int(argv[3])
    else:
        form.selfComp = 3
        form.hairpin = 3
    calcPrimer(form)
