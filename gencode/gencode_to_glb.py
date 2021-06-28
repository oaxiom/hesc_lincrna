
'''

Pack the GTF into a gneome_sql for chipFish

'''

import sys, os
from glbase3 import genelist, format

gencode = genelist('gencode.v32.annotation.gtf', format=format.gtf)
gencode.save('gencode.glb')
