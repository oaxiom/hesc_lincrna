import os, sys, glob, numpy
import matplotlib.pyplot as plot
from matplotlib import gridspec
import matplotlib.cm as cm
from glbase3 import genelist, glload, expression, config
config.draw_mode = 'pdf'

TE_lookup = glload('../../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
TE_lookup = {i['transcript_id']: i for i in TE_lookup}

dfam = glload('../../../te_discovery/dfam/dfam_annotation.glb') #, format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam = {d['name']: '{0}:{1}:{2}'.format(d['type'], d['subtype'], d['name']) for d in dfam}

res = {}
res_coding = {}
res_TEs = {}
res_TEs_type = {}
res_esc_expn = {}

for filename in sorted(glob.glob('../gls/*.glb')):
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')

    print(stub)

    genes = glload(filename)
    num_genes = len(genes)

    res[stub] = {'noTE': 0, 'TE': 0}
    res_TEs[stub] = {}
    res_TEs_type[stub] = {}
    res_esc_expn[stub] = {'ES-': 0, 'ES:': 0,  'ES+': 0,}
    res_coding[stub] = {'C': 0, 'NC': 0}

    for gene in genes:
        tags = gene['name'].split(' ')[1].split(';')
        res_esc_expn[stub][tags[2]] += 1
        if tags[1] == 'U':
            continue
        res_coding[stub][tags[1]] += 1

        if gene['transcript_id'] not in TE_lookup:
            res[stub]['noTE'] += 1
        else:
            res[stub]['TE'] += 1
            tes = TE_lookup[gene['transcript_id']]
            for TE in tes['doms']:
                if TE['dom'] not in dfam:
                    continue
                TE = dfam[TE['dom']] # Get fullname

                if True not in [i in TE for i in ('LTR', 'LINE', 'SINE', 'SVA', 'DNA')]:
                    continue
                if '.' in TE:
                    continue

                type_subtype = ':'.join(TE.split(':')[0:2])
                if type_subtype not in res_TEs_type[stub]:
                    res_TEs_type[stub][type_subtype] = 0
                res_TEs_type[stub][type_subtype] += 1

                if TE not in res_TEs[stub]:
                    res_TEs[stub][TE] = 0
                res_TEs[stub][TE] += 1

    # normalise to numebr of genes:
    for TE in res_TEs[stub]:
        res_TEs[stub][TE] /= num_genes

    for TE in res_TEs_type[stub]:
        res_TEs_type[stub][TE] /= num_genes

tab10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

heat_data = {}

for gid, grp in enumerate(res):
    fig = plot.figure(figsize=[4,1.5])
    fig.subplots_adjust(left=0.07, right=0.98, top=0.99, bottom=0.12)#, hspace=0.2, wspace=0.2)

    #ax = plot.subplot2grid((2,2), (0,0), rowspan=2)

    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=1.5, hspace=0.9)

    ax = plot.subplot(gs[0, 0])
    ax.barh([0, 1], list(res[grp].values()), color=tab10[0:2])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(res[grp].keys())
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    ax = plot.subplot(gs[1, 0])
    ax.barh([0, 1, 2], list(res_esc_expn[grp].values()), color=tab10[1:4])
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(res_esc_expn[grp].keys())
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    ax = plot.subplot(gs[2, 0])
    ax.barh([0, 1], list(res_coding[grp].values()), color=tab10[0:2])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(res_coding[grp].keys())
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    ax = plot.subplot(gs[0:, 1])
    res_TEs_type[grp] = {k: res_TEs_type[grp][k] for k in sorted(res_TEs_type[grp], key=res_TEs_type[grp].get, reverse=True)}
    res_TEs_type[grp] = {k: res_TEs_type[grp][k] for k in list(res_TEs_type[grp].keys())[0:10]}
    res_TEs_type[grp] = {k: res_TEs_type[grp][k] for k in sorted(res_TEs_type[grp], key=res_TEs_type[grp].get, reverse=False)}
    ys = numpy.arange(len(res_TEs_type[grp]))
    ax.barh(ys, list(res_TEs_type[grp].values()))
    ax.set_yticks(ys)
    ax.set_yticklabels(res_TEs_type[grp].keys())
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    for TE in res_TEs[grp]:
        if TE not in heat_data:
            heat_data[TE] = [0] * 5 # grps;
        heat_data[TE][gid] = res_TEs[grp][TE]

    # Top10:
    ax = plot.subplot(gs[0:, 2])
    res_TEs[grp] = {k: res_TEs[grp][k] for k in sorted(res_TEs[grp], key=res_TEs[grp].get, reverse=True)} # Top 10:
    res_TEs[grp] = {k: res_TEs[grp][k] for k in list(res_TEs[grp].keys())[0:10]}
    res_TEs[grp] = {k: res_TEs[grp][k] for k in sorted(res_TEs[grp], key=res_TEs[grp].get, reverse=False)} # resrot

    ys = numpy.arange(len(res_TEs[grp]))
    ax.barh(ys, list(res_TEs[grp].values()))
    ax.set_yticks(ys)
    ax.set_yticklabels(res_TEs[grp].keys())
    print('\n', grp)
    print('\n'.join(reversed(list(res_TEs[grp].keys()))))
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig('bp-{0}.png'.format(grp))
    fig.savefig('bp-{0}.pdf'.format(grp))

tes_to_keep = [
    'DNA:TcMar-Tigger',
    'DNA:hAT-Charlie',
    'DNA:hAT-Tip100',
    #'LINE:CR1',
    'LINE:L1',
    'LINE:L2',
    'LTR:ERV1',
    'LTR:ERVK',
    'LTR:ERVL',
    'LTR:ERVL-MaLR',
    #'LINE:RTE-X',
    'Retroposon:SVA',
    'SINE:Alu',
    'SINE:MIR',
    ]

tes_to_keep.reverse()

fig = plot.figure(figsize=[5,1.6])
fig.subplots_adjust(left=0.2, right=0.99)

for axidx, k in enumerate(sorted(res_TEs_type)):
    print(res_TEs_type[k])
    ax = fig.add_subplot(1, len(res_TEs_type), axidx+1)
    vals = []
    for i in tes_to_keep:
        if i in res_TEs_type[k]:
            vals.append(res_TEs_type[k][i])
        else:
            vals.append(0)

    ys = numpy.arange(len(vals))
    ax.barh(ys, vals)
    ax.set_yticks(ys)
    if axidx == 0:
        ax.set_yticklabels(tes_to_keep)
    else:
        ax.set_yticklabels('')
    ax.set_title(k)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    ax.set_xlim([0, 1.5])

fig.savefig('te_class.pdf'.format(grp))

# Heat;
e = expression(loadable_list=[{'name': k, 'conditions': heat_data[k]} for k in sorted(heat_data)], cond_names=[0, 1,2,3,4])
e = e.filter_low_expressed(0.04, 1)
e.heatmap(filename='all.pdf', col_cluster=False, row_cluster=False, heat_wid=0.12,
    heat_hei=0.011*len(e),
    bracket=[0, 0.2], border=True,
    cmap=cm.plasma)
