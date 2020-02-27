
from glbase3 import glload, genelist

gl = genelist('cell_data.csv', format={'de_clusters': 16})
print(gl)

res = {}

for item in gl:
    if item['de_clusters'] not in res:
        res[item['de_clusters']] = 0
    res[item['de_clusters']] += 1

total = len(gl)

for k in res:
    print(k, res[k], res[k]/total * 100)
print(res)
