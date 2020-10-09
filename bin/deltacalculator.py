#!/usr/bin/env python3
from sys import argv


def invertSeq(seq):
    return seq[::-1]


def complementSeqUracil(seq):
    compSeq = ''
    for index in seq:
        if index == 'A':
            compSeq += 'U'
        elif index == 'U':
            compSeq += 'A'
        elif index == 'C':
            compSeq += 'G'
        elif index == 'G':
            compSeq += 'C'
        elif index == '-':
            compSeq += '-'
    return compSeq


def acumGibbsOfEnergy(string):
    GibbsValue = 0.00
    for x in range(0, len(string) - 1):
        if((string[x] == 'A' and string[x + 1] == 'A') or (string[x] == 'A' and string[x + 1] == 'U') or (string[x] == 'U' and string[x + 1] == 'U')):
            GibbsValue += 1.10
        elif((string[x] == 'A' and string[x + 1] == 'C') or (string[x] == 'G' and string[x + 1] == 'U')):
            GibbsValue += 2.40
        elif((string[x] == 'A' and string[x + 1] == 'G') or (string[x] == 'C' and string[x + 1] == 'U')):
            GibbsValue += 1.90
        elif((string[x] == 'C' and string[x + 1] == 'A') or (string[x] == 'C' and string[x + 1] == 'G') or (string[x] == 'U' and string[x + 1] == 'G')):
            GibbsValue += 2.20
        elif((string[x] == 'C' and string[x + 1] == 'C') or (string[x] == 'G' and string[x + 1] == 'G')):
            GibbsValue += 3.30
        elif(string[x] == 'G' and string[x + 1] == 'A'):
            GibbsValue += 2.70
        elif(string[x] == 'G' and string[x + 1] == 'C'):
            GibbsValue += 3.80
        elif(string[x] == 'U' and string[x + 1] == 'A'):
            GibbsValue += 1.40
        elif(string[x] == 'U' and string[x + 1] == 'C'):
            GibbsValue += 2.60
    return GibbsValue


def deltaGibbCalculator(string):
    # Gibbs Calculator
    initialSubString = string[0:5]
    initialSubString = initialSubString
    finalSubString = string[14:19]
    finalSubString = finalSubString
    finalSubString = invertSeq(finalSubString)
    finalSubString = complementSeqUracil(finalSubString)
    initialGibbs = acumGibbsOfEnergy(initialSubString)
    finalGibbs = acumGibbsOfEnergy(finalSubString)
    total = acumGibbsOfEnergy(string)
    deltaGibbs = initialGibbs - finalGibbs
    return(deltaGibbs, initialGibbs, finalGibbs, total)


if __name__ == "__main__":
    sequence = argv[1]
    #print('qualidade\tinicio\tfim\ttotal')
    print('%.2f\t%.2f\t%.2f\t%.2f' % (deltaGibbCalculator(sequence)))
