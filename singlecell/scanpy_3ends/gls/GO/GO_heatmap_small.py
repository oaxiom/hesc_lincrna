
import glob, sys, os, math
from glbase3 import *
import matplotlib.cm as cm
config.draw_mode = 'pdf'

clus_order = [0, 1, 2, 3, 4]

format = {'force_tsv': True, 'pvalue': 1, 'name': 0}

main_cluster_membership = {}

go_to_keep = set([
    # 0 Stem cells;
    'GO:0001825 blastocyst formation',
    'GO:0001708 cell fate specification',
    'GO:0040019 positive regulation of embryonic development',
    'GO:0007369 gastrulation',
    # 1: Cell cycle;
    'GO:0045740 positive regulation of DNA replication',
    'GO:0090068 positive regulation of cell cycle process',
    # 2:
    'GO:0080154 regulation of fertilization',
    'GO:1901992 positive regulation of mitotic cell cycle phase transition',
    # 3: differnetiation;
    'GO:0050807 regulation of synapse organization',
    'GO:0072073 kidney epithelium development',
    'GO:0045669 positive regulation of osteoblast differentiation',
    # 4: Quiescence?
    'GO:0017148 negative regulation of translation',
    ])

for ont in ('BP', ):
    go_store = {}
    main_cluster_membership = {}
    for filename in glob.glob('tabs/tab{0}_de_genes-grp*.tsv'.format(ont)):
        print(filename)
        go = glgo(filename=filename, format=format)
        if not go:
            continue

        clus_number = int(os.path.split(filename)[1].split('.')[0].split('-')[-1].replace('grp', ''))

        #go.sort('1')
        top5 = go

        for item in top5:
            if item['pvalue'] < 0.01:
                if item['name'] not in go_to_keep:
                    continue

                print('Found')
                if item['name'] not in go_store:
                    go_store[item['name']] = [-1] * len(clus_order)

                go_store[item['name']][clus_number] = -math.log10(item['pvalue'])
                #main_cluster_membership[item['name']] = clus_number-1

    # fill in the holes:
    for filename in glob.glob('tabs/tab{0}_de_genes-grp*.tsv'.format(ont)):
        go = glgo(filename=filename, format=format)
        if not go:
            continue

        clus_number = int(os.path.split(filename)[1].split('.')[0].split('-')[-1].replace('grp', ''))
        for k in go_store:
            this_k = go.get(key='name', value=k, mode='lazy') # by default
            if this_k:
                print(k, clus_number)
                go_store[k][clus_number] = -math.log10(float(this_k[0]['pvalue']))

    newe = []

    print(go_store)

    for k in go_store:
        newe.append({'name': k, 'conditions': go_store[k]})

    goex = expression(loadable_list=newe, cond_names=clus_order)

    #goex = goex.sliceConditions(clus_order)
    goex = goex.filter_low_expressed(1.9, 1)

    res = goex.heatmap(filename='atmap_small_%s.png' % ont, size=[9, 8], bracket=[1.0,10],
        row_cluster=True, col_cluster=False, imshow=False,
        heat_wid=0.06, cmap=cm.Reds, border=True,
        row_font_size=7, heat_hei=0.011*len(goex), grid=True,
        draw_numbers=True, draw_numbers_fmt='*',
        draw_numbers_threshold=2.0,
        draw_numbers_font_size=6) # 1.30 = 0.05

    print('\n'.join(reversed(res['reordered_rows'])))
