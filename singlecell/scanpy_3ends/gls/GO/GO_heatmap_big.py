
import glob, sys, os, math
from glbase3 import *
import matplotlib.cm as cm
config.draw_mode = 'pdf'

clus_order = [0, 1, 2, 3, 4, 5]

format = {'force_tsv': True, 'pvalue': 1, 'name': 0}

main_cluster_membership = {}

for ont in ('BP', 'CC', 'MF'):
    go_store = {}
    main_cluster_membership = {}
    for filename in glob.glob('tabs/tab{0}_de_genes-grp*.tsv'.format(ont)):
        go = glgo(filename=filename, format=format)
        if not go:
            continue

        clus_number = int(os.path.split(filename)[1].split('.')[0].split('-')[-1].replace('grp', ''))

        #go.sort('1')
        top5 = go[0:20] # Make sure to add blastocyst and stem cell terms:

        for item in top5:
            if item['pvalue'] < 0.01:
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

    goex.heatmap(filename='atmap_big_%s.png' % ont,
        size=[9, 12], bracket=[1.0,10],
        row_cluster=True, col_cluster=False, imshow=False,
        heat_wid=0.06, cmap=cm.Reds, border=True,
        row_font_size=7, heat_hei=0.007*len(goex), grid=True,
        draw_numbers=True, draw_numbers_fmt='*',
        draw_numbers_threshold=2.0,
        draw_numbers_font_size=6) # 1.30 = 0.05


