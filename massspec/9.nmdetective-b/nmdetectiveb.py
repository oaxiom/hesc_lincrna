

from glbase3 import *

# TODO: Port into glbase gtf

def nmdb_score(self, inlastexon, orflength, exonlength, within_50nt_of_lastEJ):
    if inlastexon:
        return 0.0

    if orflength < 150:
        return 0.12

    if exonlength >407:
        return 0.41

    if within_50nt_of_lastEJ:
        nmd = 0.20

    return 0.65

for transcript in assembly:




