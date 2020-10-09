#!/usr/bin/env python3

# os 3 deltas, os 3 TMs e o  RlnK .
# faz um  pythom para isso
# e depois tem uma função lá no check self-complementarity

from sys import argv
from math import log


class Oligo:
    def __init__(self):
        self.aCount = 0
        self.cCount = 0
        self.gCount = 0
        self.tCount = 0
        self.uCount = 0
        self.iCount = 0
        self.mCount = 0
        self.rCount = 0
        self.wCount = 0
        self.sCount = 0
        self.yCount = 0
        self.kCount = 0
        self.vCount = 0
        self.hCount = 0
        self.dCount = 0
        self.bCount = 0
        self.nCount = 0

        self.gcValmin = 0
        self.gcValmax = 0
        self.mwValmin = 0
        self.mwValmax = 0

        self.concValmin = 0
        self.concValmax = 0
        self.microgramValmin = 0
        self.microgramValmax = 0

        self.deltaHValmin = 0
        self.deltaGValmin = 0
        self.deltaSValmin = 0
        self.deltaHValmax = 0
        self.deltaGValmax = 0
        self.deltaSValmax = 0

        self.RlogK = 0

        self.basicTmValmin = 0
        self.adjustedTmValmin = 0
        self.nearestNeighborTmValmin = 0
        self.basicTmValmax = 0
        self.adjustedTmValmax = 0
        self.nearestNeighborTmValmax = 0

        self.eA_A260min = 0
        self.eC_A260min = 0
        self.eG_A260min = 0
        self.eT_A260min = 0
        self.eU_A260min = 0

        self.eA_A260max = 0
        self.eC_A260max = 0
        self.eG_A260max = 0
        self.eT_A260max = 0
        self.eU_A260max = 0

        # effective counts for MW()
        self.eA_MWmin = 0
        self.eC_MWmin = 0
        self.eG_MWmin = 0
        self.eT_MWmin = 0
        self.eU_MWmin = 0
        self.eI_MWmin = 0

        self.eA_MWmax = 0
        self.eC_MWmax = 0
        self.eG_MWmax = 0
        self.eT_MWmax = 0
        self.eU_MWmax = 0
        self.eI_MWmax = 0
        self.eGC_min = 0
        self.eGC_max = 0

        self.aaCount = 0
        self.atCount = 0
        self.taCount = 0
        self.caCount = 0
        self.gtCount = 0
        self.ctCount = 0
        self.gaCount = 0
        self.cgCount = 0
        self.gcCount = 0
        self.ggCount = 0
        self.IUpairVals_min = [0, 0, 0]
        self.IUpairVals_max = [0, 0, 0]
        self.hasIUpacBase = 0
        self.isDeoxy = True
        self.isSingleStranded = True

        # arrays
        self.seqArray = []
        self.revSeqArray = []

        # publically set values
        self.Sequence = ""
        self.revSequence = ""
        self.ODs = 1
        self.FivePrimeModification = ""
        self.FivePrimeMW = 0
        self.ThreePrimeModification = ""
        self.ThreePrimeMW = 0

        self.famCount = 0
        self.tetCount = 0
        self.hexCount = 0
        self.tamraCount = 0
        self.saltConcentration = 50
        self.primerConcentration = 50


class Input:
    def __init__(
        self, oligoBox="", ODs=1,
        fivePrimeName="", threePrimeName="",
        fivePrimeValue=0, threePrimeValue=0,
        saltConcBox=50, primerConcBox=50,
        deoxy='ssRNA'
    ):
        self.oligoBox = oligoBox
        self.fivePrimeName = fivePrimeName
        self.threePrimeName = threePrimeName
        self.fivePrimeValue = fivePrimeValue
        self.threePrimeValue = threePrimeValue
        self.ODs = ODs
        self.saltConcBox = saltConcBox
        self.primerConcBox = primerConcBox
        self.deoxy = deoxy


def RemoveNonPrintingChars(theString):
    return theString.replace(' ', '')


def IsIUpacBase(theBase):
    return theBase in "MRWSYKVHDBN"


def IsBase(theBase):
    return theBase in 'ATCGU'


def CheckBase(theString):
    returnString = ""
    cnt = 0
    # rcnt = 0
    cha = ""
    theString = theString.upper()
    theString = RemoveNonPrintingChars(theString)

    for cha in theString:
        if IsIUpacBase(cha) or IsBase(cha):
            returnString += cha
            cnt += 1
        elif cha != " " and cha != "\n":
            raise ValueError(f" base {cnt+1}: {cha} is not a valid base!")

    return returnString


def MakeComplement(theSequence, isDNA):
    returnString = ""

    for temp in theSequence:
        if temp == "A":
            temp = "T" if isDNA else "U"
        elif temp == "T":
            temp = "A"
        elif temp == "U":
            temp = "A"
        elif temp == "G":
            temp = "C"
        elif temp == "C":
            temp = "G"
        elif temp == "M":
            temp = "K"
        elif temp == "K":
            temp = "M"
        elif temp == "R":
            temp = "Y"
        elif temp == "Y":
            temp = "R"
        elif temp == "W":
            temp = "W"
        elif temp == "S":
            temp = "S"
        elif temp == "V":
            temp = "B"
        elif temp == "B":
            temp = "V"
        elif temp == "H":
            temp = "D"
        elif temp == "D":
            temp = "H"

        returnString = returnString+temp

    return returnString[::-1]


def AreThereIUpacBases(theSequence):
    return any(IsIUpacBase(base) for base in theSequence)


def CountNeighbors(term, sequence):
    counter = 0
    for i in range(len(sequence)-1):
        if sequence[i:i+2] == term:
            counter += 1
    return counter


def CalcIUpair(base0, base, i, theSequence, choice):
    IUpacBase = ""
    pair1 = ""
    pair2 = ""
    temp1 = [0, 0, 0]
    temp2 = [0, 0, 0]
    reValue = [0, 0, 0]
    # print(theSequence, i)
    base2 = theSequence[i+1]

    if IsIUpacBase(base0):  # if previous base is IUpacBase, do nothing
        return reValue

    if IsIUpacBase(base):
        if base == "M":
            IUpacBase = "AC"
        elif base == "R":
            IUpacBase = "AG"
        elif base == "W":
            IUpacBase = "AT"
        elif base == "S":
            IUpacBase = "CG"
        elif base == "Y":
            IUpacBase = "CT"
        elif base == "K":
            IUpacBase = "GT"
        elif base == "V":
            IUpacBase = "ACG"
        elif base == "H":
            IUpacBase = "ACT"
        elif base == "D":
            IUpacBase = "AGT"
        elif base == "B":
            IUpacBase = "CGT"
        elif base == "N":
            IUpacBase = "ACGT"

        j = 0
        while IUpacBase[j] != "":
            base = IUpacBase[j]
            pair1 = base0+base

            if pair1 == "AA":
                temp1[0] = 1.2
                temp1[1] = 8.0
                temp1[2] = 21.9

            elif pair1 == "AT":
                temp1[0] = 0.9
                temp1[1] = 5.6
                temp1[2] = 15.2
            elif pair1 == "TA":
                temp1[0] = 0.9
                temp1[1] = 6.6
                temp1[2] = 18.4
            elif pair1 == "CA":
                temp1[0] = 1.7
                temp1[1] = 8.2
                temp1[2] = 21.0
            elif pair1 == "GT":
                temp1[0] = 1.5
                temp1[1] = 9.4
                temp1[2] = 25.5
            elif pair1 == "CT":
                temp1[0] = 1.5
                temp1[1] = 6.6
                temp1[2] = 16.4
            elif pair1 == "GA":
                temp1[0] = 1.5
                temp1[1] = 8.8
                temp1[2] = 23.5
            elif pair1 == "CG":
                temp1[0] = 2.8
                temp1[1] = 11.8
                temp1[2] = 29.0
            elif pair1 == "GC":
                temp1[0] = 2.3
                temp1[1] = 10.5
                temp1[2] = 26.4
            elif pair1 == "GG":
                temp1[0] = 2.1
                temp1[1] = 10.9
                temp1[2] = 28.4

            if base2 == "":
                for k in range(2):
                    temp2[k] = 0.0
            elif not IsIUpacBase(base2):
                pair2 = base+base2
                if pair2 == "AA":
                    temp2[0] = 1.2
                    temp2[1] = 8.0
                    temp2[2] = 21.9
                elif pair2 == "AT":
                    temp2[0] = 0.9
                    temp2[1] = 5.6
                    temp2[2] = 15.2
                elif pair2 == "TA":
                    temp2[0] = 0.9
                    temp2[1] = 6.6
                    temp2[2] = 18.4
                elif pair2 == "CA":
                    temp2[0] = 1.7
                    temp2[1] = 8.2
                    temp2[2] = 21.0
                elif pair2 == "GT":
                    temp2[0] = 1.5
                    temp2[1] = 9.4
                    temp2[2] = 25.5
                elif pair2 == "CT":
                    temp2[0] = 1.5
                    temp2[1] = 6.6
                    temp2[2] = 16.4
                elif pair2 == "GA":
                    temp2[0] = 1.5
                    temp2[1] = 8.8
                    temp2[2] = 23.5
                elif pair2 == "CG":
                    temp2[0] = 2.8
                    temp2[1] = 11.8
                    temp2[2] = 29.0
                elif pair2 == "GC":
                    temp2[0] = 2.3
                    temp2[1] = 10.5
                    temp2[2] = 26.4
                elif pair2 == "GG":
                    temp2[0] = 2.1
                    temp2[1] = 10.9
                    temp2[2] = 28.4
            elif IsIUpacBase(base2):
                base0 = base
                base = base2
                i = i + 1
                temp2 = CalcIUpair(base0, base, i, theSequence, choice)
                i = i - 1

            for k in range(3):
                if j == 0:
                    reValue[k] = temp1[k]+temp2[k]
                else:
                    if (choice == "max") and (reValue[k] < (temp1[k]+temp2[k])):
                        reValue[k] = temp1[k]+temp2[k]
                    elif(choice == "min") and (reValue[k] > (temp1[k]+temp2[k])):
                        reValue[k] = temp1[k]+temp2[k]

            j += 1

    return reValue


def OligoCount(oligoObj):
    oligoObj.aCount = oligoObj.Sequence.count("A")
    oligoObj.cCount = oligoObj.Sequence.count("C")
    oligoObj.gCount = oligoObj.Sequence.count("G")
    oligoObj.tCount = oligoObj.Sequence.count("T")
    oligoObj.uCount = oligoObj.Sequence.count("U")
    oligoObj.iCount = oligoObj.Sequence.count("I")
    oligoObj.mCount = oligoObj.Sequence.count("M")
    oligoObj.rCount = oligoObj.Sequence.count("R")
    oligoObj.wCount = oligoObj.Sequence.count("W")
    oligoObj.sCount = oligoObj.Sequence.count("S")
    oligoObj.yCount = oligoObj.Sequence.count("Y")
    oligoObj.kCount = oligoObj.Sequence.count("K")
    oligoObj.vCount = oligoObj.Sequence.count("V")
    oligoObj.hCount = oligoObj.Sequence.count("H")
    oligoObj.dCount = oligoObj.Sequence.count("D")
    oligoObj.bCount = oligoObj.Sequence.count("B")
    oligoObj.nCount = oligoObj.Sequence.count("N")

    # Effective a, c, g, t count for different calculations
    # effective counts for A260min()

    oligoObj.eA_A260min = oligoObj.aCount
    oligoObj.eU_A260min = oligoObj.uCount
    oligoObj.eI_A260min = oligoObj.iCount
    oligoObj.eC_A260min = sum([
        oligoObj.cCount, oligoObj.mCount, oligoObj.sCount,
        oligoObj.yCount, oligoObj.vCount, oligoObj.hCount,
        oligoObj.bCount, oligoObj.nCount
    ])
    oligoObj.eG_A260min = oligoObj.gCount+oligoObj.rCount
    oligoObj.eT_A260min = oligoObj.tCount + \
        oligoObj.wCount+oligoObj.kCount+oligoObj.dCount
    oligoObj.eA_A260max = oligoObj.aCount+oligoObj.mCount+oligoObj.rCount + \
        oligoObj.wCount+oligoObj.vCount+oligoObj.hCount+oligoObj.bCount+oligoObj.nCount
    oligoObj.eC_A260max = oligoObj.cCount
    oligoObj.eG_A260max = oligoObj.gCount + \
        oligoObj.sCount+oligoObj.kCount+oligoObj.dCount
    oligoObj.eT_A260max = oligoObj.tCount+oligoObj.yCount

    # effective counts for MW()
    oligoObj.eA_MWmin = oligoObj.aCount+oligoObj.rCount
    oligoObj.eC_MWmin = oligoObj.eC_A260min
    oligoObj.eU_MWmin = oligoObj.eU_A260min
    oligoObj.eI_MWmin = oligoObj.eI_A260min
    oligoObj.eG_MWmin = oligoObj.gCount
    oligoObj.eT_MWmin = oligoObj.eT_A260min
    oligoObj.eA_MWmax = oligoObj.aCount + \
        oligoObj.mCount+oligoObj.wCount+oligoObj.hCount
    oligoObj.eC_MWmax = oligoObj.cCount
    oligoObj.eG_MWmax = oligoObj.gCount+oligoObj.rCount+oligoObj.sCount + \
        oligoObj.kCount+oligoObj.vCount + oligoObj.dCount+oligoObj.bCount+oligoObj.nCount
    oligoObj.eT_MWmax = oligoObj.tCount+oligoObj.yCount
    oligoObj.eGC_min = oligoObj.gCount + oligoObj.cCount + oligoObj.sCount
    oligoObj.eGC_max = (len(oligoObj.Sequence))-oligoObj.aCount - \
        oligoObj.tCount-oligoObj.wCount-oligoObj.uCount

    # count Nearest Neighbors
    oligoObj.aaCount = CountNeighbors("AA", oligoObj.Sequence)+CountNeighbors(
        "TT", oligoObj.Sequence)+CountNeighbors("UU", oligoObj.Sequence)
    oligoObj.atCount = CountNeighbors(
        "AT", oligoObj.Sequence)+CountNeighbors("AU", oligoObj.Sequence)
    oligoObj.taCount = CountNeighbors(
        "TA", oligoObj.Sequence)+CountNeighbors("UA", oligoObj.Sequence)
    oligoObj.caCount = CountNeighbors("CA", oligoObj.Sequence)+CountNeighbors(
        "TG", oligoObj.Sequence)+CountNeighbors("UG", oligoObj.Sequence)
    oligoObj.gtCount = CountNeighbors("GT", oligoObj.Sequence)+CountNeighbors(
        "AC", oligoObj.Sequence)+CountNeighbors("GU", oligoObj.Sequence)
    oligoObj.ctCount = CountNeighbors("CT", oligoObj.Sequence)+CountNeighbors(
        "AG", oligoObj.Sequence)+CountNeighbors("CU", oligoObj.Sequence)
    oligoObj.gaCount = CountNeighbors("GA", oligoObj.Sequence)+CountNeighbors(
        "TC", oligoObj.Sequence)+CountNeighbors("UC", oligoObj.Sequence)
    oligoObj.cgCount = CountNeighbors("CG", oligoObj.Sequence)
    oligoObj.gcCount = CountNeighbors("GC", oligoObj.Sequence)
    oligoObj.ggCount = CountNeighbors(
        "GG", oligoObj.Sequence)+CountNeighbors("CC", oligoObj.Sequence)

    # Calculate IUpac pairs
    # 08/02/99 fix one bug for nearest Neighbor Calc,

    for j in range(3):
        oligoObj.IUpairVals_min[j] = 0
        oligoObj.IUpairVals_max[j] = 0

        for i in range(1, len(oligoObj.Sequence)-1):
            # first base can not be IUpacbase
            base0 = oligoObj.Sequence[i-1]
            base = oligoObj.Sequence[i]
            temp = []
            temp = CalcIUpair(base0, base, i, oligoObj.Sequence, "min")

            for k in range(3):
                oligoObj.IUpairVals_min[k] += temp[k]

            temp = CalcIUpair(base0, base, i, oligoObj.Sequence, "max")

            for k in range(3):
                oligoObj.IUpairVals_max[j] += temp[j]

    return oligoObj


def GC(oligoObj, choice):
    len_seq = len(oligoObj.Sequence)
    if len_seq > 0:
        if choice == "min":
            return round(100 * oligoObj.eGC_min / len_seq)
        else:
            return round(100*oligoObj.eGC_max / len_seq)
    raise ValueError("You have passed an empty sequence!")


def micrograms(MolWt, Conc):
    # /* MolWt is gms/mol Conc is micromoles/L assume volume is 1 milliliter */
    if MolWt > 0 and Conc > 0:
        return round((MolWt*Conc/100)/10)
    raise ZeroDivisionError("Either MolWt either Conc is Zero!")


def MW(oligoObj, choice):
    mw = None
    if (len(oligoObj.Sequence) > 0):
        if oligoObj.isDeoxy:
            if choice == "min":
                mw = 313.21 * oligoObj.eA_MWmin + 329.21 * oligoObj.eG_MWmin + 289.18 * oligoObj.eC_MWmin + \
                    304.2 * oligoObj.eT_MWmin + 290.169 * \
                    oligoObj.eU_MWmin + 314. * oligoObj.eI_MWmin - 61.96
            else:
                mw = 313.21 * oligoObj.eA_MWmax + 329.21 * oligoObj.eG_MWmax + 289.18 * oligoObj.eC_MWmax + \
                    304.2 * oligoObj.eT_MWmax + 290.169 * \
                    oligoObj.eU_MWmax + 314. * oligoObj.eI_MWmax - 61.96
        else:  # // is riboNucleotide
            if choice == "min":
                mw = 329.21 * oligoObj.eA_MWmin + 345.21 * oligoObj.eG_MWmin + 305.18 * oligoObj.eC_MWmin + \
                    320.2 * oligoObj.eT_MWmin + 306.169 * \
                    oligoObj.eU_MWmin + 330 * oligoObj.eI_MWmin+159
            else:
                mw = 329.21 * oligoObj.eA_MWmax + 345.21 * oligoObj.eG_MWmax + 305.18 * oligoObj.eC_MWmax + \
                    320.2 * oligoObj.eT_MWmax + 306.169 * \
                    oligoObj.eU_MWmax + 330 * oligoObj.eI_MWmax+159
        new_mw = mw+oligoObj.FivePrimeMW
        new_mw += oligoObj.ThreePrimeMW
        new_mw = round(10 * new_mw)/10
        return new_mw
    raise ValueError("You have passed an empty sequence!")


def A260(oligoObj, choice):
    div = None
    if (len(oligoObj.Sequence) > 0):
        if (oligoObj.isSingleStranded):
            if (oligoObj.isDeoxy):
                if (choice == "min"):  # calculates the minimum value, use the e_A260max since it is the divident
                    div = (oligoObj.eA_A260max * 15200 + oligoObj.eG_A260max * 12010 + oligoObj.eC_A260max * 7050 + oligoObj.eT_A260max * 8400 +
                           oligoObj.eU_A260max * 9800 + oligoObj.famCount * 20960 + oligoObj.tetCount * 16255 + oligoObj.hexCount * 31580 + oligoObj.tamraCount * 31980)
                else:  # choice == "max"
                    div = (oligoObj.eA_A260min * 15200 + oligoObj.eG_A260min * 12010 + oligoObj.eC_A260min * 7050 + oligoObj.eT_A260min * 8400 +
                           oligoObj.eU_A260min * 9800 + oligoObj.famCount * 20960 + oligoObj.tetCount * 16255 + oligoObj.hexCount * 31580 + oligoObj.tamraCount * 31980)

            else:
                if (choice == "min"):  # calculates the minimum value, use the e_A260max since it is the divident
                    div = (oligoObj.eA_A260max * 15400 + oligoObj.eG_A260max * 13700 + oligoObj.eC_A260max * 9000 + oligoObj.eT_A260max * 9400 +
                           oligoObj.eU_A260max * 10000 + oligoObj.famCount * 20960 + oligoObj.tetCount * 16255 + oligoObj.hexCount * 31580 + oligoObj.tamraCount * 31980)
                else:  # choice == "max"
                    div = (oligoObj.eA_A260min * 15400 + oligoObj.eG_A260min * 13700 + oligoObj.eC_A260min * 9000 + oligoObj.eT_A260min * 9400 +
                           oligoObj.eU_A260min * 10000 + oligoObj.famCount * 20960 + oligoObj.tetCount * 16255 + oligoObj.hexCount * 31580 + oligoObj.tamraCount * 31980)
        else:
            if oligoObj.isDeoxy:
                div = oligoObj.Sequence.length * 2 * 6400
            else:
                div = oligoObj.Sequence.length * 2 * 8000

        # units are in microMoles/liter
        return (round(oligoObj.ODs * 1000000000 / div) / 1000)

    raise ValueError("You have passed an empty sequence!")


def OligoCalcMWConcAndODs(oligoObj):
    if not oligoObj.hasIUpacBase:
        oligoObj.mwValmin = MW(oligoObj, "min")
        oligoObj.mwValmax = oligoObj.mwValmin
    else:
        oligoObj.mwValmin = MW(oligoObj, "min")
        oligoObj.mwValmax = MW(oligoObj, "max")

        # Now do A260 concentration calculation
        oligoObj.concValmin = A260(oligoObj, "min")
        oligoObj.concValmax = A260(oligoObj, "max")
        oligoObj.microgramValmin = micrograms(
            oligoObj.mwValmin, oligoObj.concValmin
        )  # external call
        oligoObj.microgramValmax = micrograms(
            oligoObj.mwValmax, oligoObj.concValmax
        )


def DeltaH(oligoObj, choice):
    if len(oligoObj.Sequence) > 7:
        val = 0.0
        if oligoObj.isDeoxy:
            val += 8.0*oligoObj.aaCount
            val += 5.6*oligoObj.atCount
            val += 6.6*oligoObj.taCount
            val += 8.2*oligoObj.caCount
            val += 9.4*oligoObj.gtCount
            val += 6.6*oligoObj.ctCount
            val += 8.8*oligoObj.gaCount
            val += 11.8*oligoObj.cgCount
            val += 10.5*oligoObj.gcCount
            val += 10.9*oligoObj.ggCount
        else:
            val += 6.8*oligoObj.aaCount
            val += 9.38*oligoObj.atCount
            val += 7.69*oligoObj.taCount
            val += 10.44*oligoObj.caCount
            val += 11.4*oligoObj.gtCount
            val += 10.48*oligoObj.ctCount
            val += 12.44*oligoObj.gaCount
            val += 10.64*oligoObj.cgCount
            val += 14.88*oligoObj.gcCount
            val += 13.39*oligoObj.ggCount

        if choice == "min":
            val += oligoObj.IUpairVals_min[1]
        else:
            val += oligoObj.IUpairVals_max[1]

        return round((1000*val))/1000

    raise ValueError("Length of sequence is lesser than 8nt")


def DeltaG(oligoObj, choice):
    if len(oligoObj.Sequence) > 7:
        val = -5.0
        # Helix initiation Free Energy of 5 kcal.
        # symmetry function: if symmetrical, subtract another 0.4
        val += 1.2*oligoObj.aaCount
        val += 0.9*oligoObj.atCount
        val += 0.9*oligoObj.taCount
        val += 1.7*oligoObj.caCount
        val += 1.5*oligoObj.gtCount
        val += 1.5*oligoObj.ctCount
        val += 1.5*oligoObj.gaCount
        val += 2.8*oligoObj.cgCount
        val += 2.3*oligoObj.gcCount
        val += 2.1*oligoObj.ggCount

        if choice == "min":
            val += oligoObj.IUpairVals_min[0]
        else:
            val += oligoObj.IUpairVals_max[0]

        return round((1000*val))/1000

    raise ValueError("Length of sequence is lesser than 8nt")


def DeltaS(oligoObj, choice):
    if len(oligoObj.Sequence) > 7:
        val = 0
        if oligoObj.isDeoxy:
            val += 21.9*oligoObj.aaCount
            val += 15.2*oligoObj.atCount
            val += 18.4*oligoObj.taCount
            val += 21.0*oligoObj.caCount
            val += 25.5*oligoObj.gtCount
            val += 16.4*oligoObj.ctCount
            val += 23.5*oligoObj.gaCount
            val += 29.0*oligoObj.cgCount
            val += 26.4*oligoObj.gcCount
            val += 28.4*oligoObj.ggCount
        else:
            val += 19.0*oligoObj.aaCount
            val += 26.7*oligoObj.atCount
            val += 20.5*oligoObj.taCount
            val += 26.9*oligoObj.caCount
            val += 29.5*oligoObj.gtCount
            val += 27.1*oligoObj.ctCount
            val += 32.5*oligoObj.gaCount
            val += 26.7*oligoObj.cgCount
            val += 36.9*oligoObj.gcCount
            val += 32.7*oligoObj.ggCount

        if choice == "min":
            val += oligoObj.IUpairVals_min[2]
        else:
            val += oligoObj.IUpairVals_max[2]
        return round((1000*val))/1000

    raise ValueError("Length of sequence is lesser than 8nt")


def Tm(oligoObj, choice):
    if len(oligoObj.Sequence) > 0:
        if len(oligoObj.Sequence) < 14:
            if choice == "min":
                return round(
                    2 * (len(oligoObj.Sequence)-oligoObj.eGC_min) +
                    4 * (oligoObj.eGC_min)
                )
            else:
                return round(
                    2 * (len(oligoObj.Sequence)-oligoObj.eGC_max) +
                    4 * (oligoObj.eGC_max)
                )
        else:
            if choice == "min":
                return round(
                    64.9 + 41*((oligoObj.eGC_min - 16.4) /
                               len(oligoObj.Sequence))
                )
            else:
                return round(
                    64.9 + 41*((oligoObj.eGC_max - 16.4) /
                               len(oligoObj.Sequence))
                )
    raise ValueError("You have passed an empty sequence!")


def WAKTm(oligoObj, choice):
    if len(oligoObj.Sequence) > 0:
        if oligoObj.isDeoxy:
            if len(oligoObj.Sequence) < 14:
                if choice == "min":
                    return round(2 * (len(oligoObj.Sequence)-oligoObj.eGC_min) + 4 * (oligoObj.eGC_min) + 21.6+(7.21*log(oligoObj.saltConcentration/1000)))
                else:
                    return round(2 * (len(oligoObj.Sequence)-oligoObj.eGC_max) + 4 * (oligoObj.eGC_max) + 21.6+(7.21*log(oligoObj.saltConcentration/1000)))
            else:
                if choice == "min":
                    return round(100.5 + (0.41*oligoObj.gcValmin) - (820 / len(oligoObj.Sequence)) + (7.21*log(oligoObj.saltConcentration/1000)))
                else:
                    return round(100.5 + (0.41*oligoObj.gcValmax) - (820 / len(oligoObj.Sequence)) + (7.21*log(oligoObj.saltConcentration/1000)))
        else:
            if len(oligoObj.Sequence) > 17:
                if choice == "min":
                    return round(79.8 + (0.584*oligoObj.gcValmin) + (11.8*(oligoObj.gcValmin/100)*(oligoObj.gcValmin/100)) - (820 / len(oligoObj.Sequence)) + (8.03*log(oligoObj.saltConcentration/1000)))
                else:
                    return round(79.8 + (0.584*oligoObj.gcValmax) - (820 / len(oligoObj.Sequence)) + (8.03*log(oligoObj.saltConcentration/1000)))
            else:
                raise ValueError("You have passed a sequence too short!")
    raise ValueError("You have passed an empty sequence!")


def NearestNeighborTM(oligoObj, choice):
    theReturn = ""
    if (len(oligoObj.Sequence) > 7):
        # Convert from nanomoles to moles
        K = 1/(oligoObj.primerConcentration*1E-9)
        R = 1.987  # cal/(mole*K)
        RlnK = R*log(K)  # javascript log is the natural log (ln)
        oligoObj.RlogK = round(1000 * RlnK)/1000
        # // Helix initiation Free Energy of 3.4 kcal (Sugimoto et al, 1996)
        # // symmetry function: if symmetrical, subtract another 0.4 - not implemented
        if (choice == "min"):
            theReturn = 1000*((DeltaH(oligoObj, "min")-3.4) /
                              (DeltaS(oligoObj, "min")+RlnK))
            theReturn += -272.9  # Kelvin to Centigrade
            theReturn += 7.21*log(oligoObj.saltConcentration/1000)
            theReturn = round(theReturn)
        else:
            theReturn = 1000*((DeltaH(oligoObj, "max")-3.4) /
                              (DeltaS(oligoObj, "max")+RlnK))
            theReturn += -272.9  # Kelvin to Centigrade
            theReturn += 7.21*log(oligoObj.saltConcentration/1000)
            theReturn = round(theReturn)
    else:
        oligoObj.RlogK = ""

    return theReturn


def OligoCalculate(oligoObj):
    oligoObj = OligoCount(oligoObj)

    # Do GC calculation
    if(oligoObj.hasIUpacBase):
        oligoObj.gcValmin = GC(oligoObj, "min")
        oligoObj.gcValmax = oligoObj.gcValmin
    else:
        oligoObj.gcValmin = GC(oligoObj, "min")
        oligoObj.gcValmax = GC(oligoObj, "max")

    OligoCalcMWConcAndODs(oligoObj)

    # Calculate numbers for Thermodynamic TM calculations
    if not oligoObj.hasIUpacBase:
        oligoObj.deltaHValmin = DeltaH(oligoObj, "min")
        oligoObj.deltaGValmin = DeltaG(oligoObj, "min")
        oligoObj.deltaSValmin = DeltaS(oligoObj, "min")
        oligoObj.deltaHValmax = oligoObj.deltaHValmin
        oligoObj.deltaGValmax = oligoObj.deltaGValmin
        oligoObj.deltaSValmax = oligoObj.deltaSValmin
    else:
        oligoObj.deltaHValmin = DeltaH(oligoObj, "min")
        oligoObj.deltaGValmin = DeltaG(oligoObj, "min")
        oligoObj.deltaSValmin = DeltaS(oligoObj, "min")
        oligoObj.deltaHValmax = DeltaH(oligoObj, "max")
        oligoObj.deltaGValmax = DeltaG(oligoObj, "max")
        oligoObj.deltaSValmax = DeltaS(oligoObj, "max")

    if not oligoObj.hasIUpacBase:
        oligoObj.basicTmValmin = Tm(oligoObj, "min")
        oligoObj.adjustedTmValmin = WAKTm(oligoObj, "min")
        oligoObj.nearestNeighborTmValmin = NearestNeighborTM(oligoObj, "min")
        oligoObj.basicTmValmax = oligoObj.basicTmValmin
        oligoObj.adjustedTmValmax = oligoObj.adjustedTmValmin
        oligoObj.nearestNeighborTmValmax = oligoObj.adjustedTmValmin
    else:
        oligoObj.basicTmValmin = Tm(oligoObj, "min")
        oligoObj.adjustedTmValmin = WAKTm(oligoObj, "min")
        oligoObj.nearestNeighborTmValmin = NearestNeighborTM(oligoObj, "min")
        oligoObj.basicTmValmax = Tm(oligoObj, "max")
        oligoObj.adjustedTmValmax = WAKTm(oligoObj, "max")
        oligoObj.nearestNeighborTmValmax = NearestNeighborTM(oligoObj, "max")

    return oligoObj


def DoOligoCalc(form, theSequence, oligoObj):
    temp = CheckBase(theSequence)  # external call
    # if (temp==-1) {return} // do not continue if the sequence is not valid!

    oligoObj.Sequence = temp
    oligoObj.ODs = form.ODs
    oligoObj.saltConcentration = form.saltConcBox
    oligoObj.primerConcentration = form.primerConcBox

    oligoObj.hasIUpacBase = AreThereIUpacBases(oligoObj.Sequence)

    # if (debug) alert("Sequence length="+len(oligoObj.Sequence)+" Sequence="+oligoObj.Sequence)

    # temp = form.deoxy

    oligoObj.isDeoxy = 'DNA' in form.deoxy
    oligoObj.isSingleStranded = 'ss' in form.deoxy

    oligoObj = OligoCalculate(oligoObj)

    return oligoObj


def GetOligoMods(form, oligoObj):
    oligoObj.FivePrimeModification = form.fivePrimeName
    oligoObj.FivePrimeMW = form.fivePrimeValue
    oligoObj.ThreePrimeModification = form.threePrimeName
    oligoObj.ThreePrimeMW = form.threePrimeValue

    return oligoObj


def Calculate(form, theOligo, theComplement):
    temp = None
    CompString = None

    temp = CheckBase(form.oligoBox)

    if temp and len(form.oligoBox) == 0:
        form.complementBox = ""
        raise ValueError("Empty sequence provided")

    form.oligoBox = temp
    theOligo = DoOligoCalc(form, temp, theOligo)
    theOligo = GetOligoMods(form, theOligo)
    CompString = MakeComplement(form.oligoBox, theOligo.isDeoxy)
    theComplement = DoOligoCalc(form, CompString, theComplement)

    if not theOligo.isSingleStranded:
        # calculate the ds MW now that the complement has been made
        theOligo = DoOligoCalc(form, temp, theOligo)

    return theOligo, theComplement


def Start(form):
    theOligo = Oligo()  # instantiate Oligo object!
    theComplement = Oligo()
    return Calculate(form, theOligo, theComplement)


if __name__ == "__main__":
    sequence = argv[1]
    oligo, complement = Start(Input(oligoBox=sequence))
    fields = [
        'Tm', 'TmSalt', 'TmNearestNeighbor',
        'RlogK', 'deltaG', 'deltaH', 'deltaS'
    ]
    #print('\t'.join(fields))
    print('\t'.join([
        str(oligo.basicTmValmin),
        str(oligo.adjustedTmValmin),
        str(oligo.nearestNeighborTmValmin),
        str(oligo.RlogK),
        str(oligo.deltaGValmin),
        str(oligo.deltaHValmin),
        str(oligo.deltaSValmin),
    ]))
