#!/usr/bin/env python3
from sys import argv


def calculate_tm(x, w, y, z):
    return 64.9 + 41 * ((y+z-16.4)/(x+y+w+z))


if __name__ == "__main__":
    sequence = argv[1]
    sequence = sequence.replace('U', 'T')
    x, w, y, z = [sequence.count(nucl) for nucl in 'TAGC']
    print("%.2f" % calculate_tm(x, w, y, z))
