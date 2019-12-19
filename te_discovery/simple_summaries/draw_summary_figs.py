
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

import pies
sys.path.append('../../')
import shared
from glbase3 import glload, utils, expression, genelist

gtf = shared.get_pickle('current_gtf/results.pickle')
gencode = shared.get_pickle('gencode/results.pickle')

print(gtf)
print(gencode)

data = {
    'hPSC': gtf['ncrna'],
    'GENCODE': gencode['ncrna']
    }

data = {
    'hPSC': gtf['pc'],

    'GENCODE': gencode['pc']
    }
