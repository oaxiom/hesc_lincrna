
import glob, sys, os
from glbase3 import genelist, glload

for f in glob.glob('../*.glb'):
    stub = os.path.split(f)[1].replace('.glb', '')
    gl = glload(f)

    for item in gl:
        item['ensg'] = item['ensg'].split('.')[0]

    gl._optimiseData()
    gl = gl.removeDuplicates('ensg')
    gl.saveTSV('{0}.tsv'.format(stub), key_order=['ensg'])
