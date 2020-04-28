import os, sys, glob, numpy
import matplotlib.pyplot as plot
from matplotlib import gridspec
from glbase3 import genelist, glload

TE_lookup = glload('../../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
TE_lookup = {i['transcript_id']: i for i in TE_lookup}

dfam = glload('../../../te_discovery/dfam/dfam_annotation.glb') #, format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam = {d['name']: '{0}:{1}:{2}'.format(d['type'], d['subtype'], d['name']) for d in dfam}

res = {}
res_coding = {}
res_TEs = {}
res_TEs_type = {}
res_esc_expn = {}

for filename in glob.glob('../gls/*.glb'):
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')

    genes = glload(filename)

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

                if True not in [i in TE for i in ('LTR', 'LINE', 'SINE', 'SVA')]:
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

tab10 = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for grp in res:
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
