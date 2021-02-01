
'''

# Filter HMMS for only relevant human/primate ones;

'''

import glob, sys, os, gzip
from operator import itemgetter
from glbase3 import utils, expression, genelist, glload, config

dfam = genelist('dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'HMMlength': 1, 'species': 2, 'type': 3, 'subtype': 4})

reject_species = set([
    'Caenorhabditis_elegans',
    'Chrysemys', # Turtles
    'Chrysemys_picta_bellii',
    'Chrysochloris_asiatica', # Golden mole;
    'Cyprinidae',
    'Danio',
    'Danio_rerio',
    'Drosophila_fruit_fly_genus',
    'Drosophila_melanogaster',
    'Durocryptodira',
    'Ficedula_albicollis',
    'Halyomorpha_halys',
    'Heliconiini',
    'Heliconius',
    'Muridae',
    'Murinae',
    'Mus_mouse_genus',
    'Mus_musculus',
    'Protostomia', # Humans are Deuterostomia
    'Rodentia',
    'Sauropsida',
    'Sciuromorpha',
    'Teleostei',
    'Testudines',
    'Testudinoidea',
    'Theria_Mammalia',
    'Uraeginthus_cyanocephalus',
    ])

all_species = sorted(list(set(dfam['species'])))

print('\n'.join(all_species))

new_dfam = []

for hmm in dfam:
    if hmm['species'] in reject_species: # Only keep Hominid TEs
        continue
    else:
        new_dfam.append(hmm)

gl = genelist()
gl.load_list(new_dfam)
gl.saveTSV('filtered.dfams.tsv')
