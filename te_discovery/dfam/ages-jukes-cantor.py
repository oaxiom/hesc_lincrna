import math
from glbase3 import *

# Jukes-Cantor ages for TEs, a simple estiamte of TE age;

res = {}

oh = open('hg38_rmsk.tsv', 'rt')

for l in oh:
    if 'repName' in l:
        continue
    t = l.split('\t')
    full_name = '{0}:{1}:{2}'.format(t[11], t[12], t[10])
    if full_name not in res:
        res[full_name] = []
    res[full_name].append(float(t[2])) # milliDiv score
print('Found {0:,} TEs'.format(len(res)))

# get averages, jk:
jk = {}
for k in res:
    if res[k]:
        #res[k] = max(res[k]) / 1000.0
        res[k] = (sum(res[k]) / len(res[k])) / 1000.0 # milliDiv is parts per thousand, this should be in range 0 .. 1
    substitutions_fraction = res[k]
    if substitutions_fraction:
        jk[k] = -3.0/4.0 * math.log(1 - 4.0/3.0 * substitutions_fraction)
        #jk[k] /= 4.5e-9 # Mm substitution rate
        jk[k] /= 2.2e-9 # Hs substitution rate

print('Ages:')
for k in jk:
    print('%s: \t%s' % (k, "{:,}".format(int(jk[k]))))

gl = genelist()
ll = [{'TE': k, 'age': jk[k]} for k in jk]
gl.load_list(ll)
gl.sort('TE')
gl.save('te_ages.glb')
gl.saveTSV('te_ages.tsv', key_order=['TE', 'age'])

# Just for speed;
# Note that this version has not been filtered; For that, see the one in
dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam.save('dfam_annotation.glb')
