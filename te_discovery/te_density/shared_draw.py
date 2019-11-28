
import math, numpy, sys
import matplotlib.pyplot as plot
import matplotlib.cm as cm

sys.path.append('../../')
import shared

def draw_density(filename, selected_genes, selected_genes2, TE=None, threshold=20):

    res = [numpy.zeros(1000), numpy.zeros(1000)] # scaled bins to put in;
    for ridx, sl in enumerate([selected_genes, selected_genes2]):
        for n, gene in enumerate(sl):
            tlen = shared.convert_genocode_to_local(gene)[1]

            for d in gene['doms']:
                if TE:
                    if d['dom'] not in TE:
                        continue

                s = math.floor(d['span'][0] / tlen * 1000)
                e = math.ceil(d['span'][1] / tlen * 1000)

                res[ridx][s:e] += 1

    if max(res[0]) < threshold and max(res[1]) < threshold:
        return None # stop from drawing

    # # norm values to the number of transcripts in total set;
    res[0] /= len(selected_genes)
    res[1] /= len(selected_genes2)

    fig = plot.figure(figsize=[2,1.3])
    fig.subplots_adjust(left=0.30, bottom=0.35,)
    ax = fig.add_subplot(111)

    ax.plot(res[0], label=selected_genes.name)
    ax.plot(res[1], label=selected_genes2.name)
    ax.tick_params(labelsize=6)
    ax.set_xticks([0, 1000])
    ax.set_xticklabels(['TSS', 'TTS'])
    #fig.savefig(filename)
    #ax.legend()
    fig.savefig(filename.replace('.png', '.pdf'))
    plot.close(fig)

    print('Saved %s' % filename)
    return res

def draw_heatmap(filename, res):
    # should be just one set of labels...

    res_labels = sorted(res.keys(), reverse=True)

    # normalise each row;
    for k in res_labels:
        s = max(res[k][0].std(), res[k][1].std())
        m = max(res[k][0].mean(), res[k][1].mean())
        res[k][0] -= m
        res[k][0] /= s

        res[k][1] -= m
        res[k][1] /= s

    fig = plot.figure(figsize=[8,9])
    heat_hei = 0.010*len(res_labels)
    fig.subplots_adjust(left=0.2, right=0.75, bottom=0.95-heat_hei, top=0.95, wspace=0.1)

    for i in [0,1]: # only 2 heatmaps supported at the moment;
        res_tab = numpy.array([res[k][i] for k in res_labels])

        ax = fig.add_subplot(1,2,i+1)
        #ax.set_position([0.5, 0.1, 0.48, 0.8])

        ax.imshow(res_tab, cmap=cm.PuOr_r, aspect="auto",
            origin='lower',
            extent=[0, res_tab.shape[1], 0, res_tab.shape[0]])

        ax.set_yticklabels(res_labels)
        if i>0:
            ax.set_yticklabels('')
        ax.set_yticks(numpy.arange(len(res_labels))+0.5)

        ax.set_xticklabels('')
        ax.tick_params(labelsize=6, bottom=False)
        #plot.colorbar(cax=ax)

    #fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))

def draw_density_utrs(filename, known, novel, bkgd, TE=None):
    data = {
        'known': {
            'utr5': numpy.zeros(1000), # scaled bins to put in;
            'cds': numpy.zeros(1000),
            'utr3': numpy.zeros(1000),
            },
        'novel': {
            'utr5': numpy.zeros(1000), # scaled bins to put in;
            'cds': numpy.zeros(1000),
            'utr3': numpy.zeros(1000),
            },
        'bkgd': {
            'utr5': numpy.zeros(1000), # scaled bins to put in;
            'cds': numpy.zeros(1000),
            'utr3': numpy.zeros(1000),
            }
        }

    genes = {'known': selected_genes,
        'novel': bkgd1,
        'bkgd': bkgd2}

    for gl in genes:
        for n, gene in enumerate(genes[gl]):
            #print(gene)
            # scale the TE to the mRNA
            if gl != 'novel':
                pos = shared.convert_genocode_to_local(gene) # return 0, tlength, cdsl, cdsr, splice_sites
            else:
                # Have to get the pos from the gene entry
                pos = (0, 100, 20, 50, [])
                print(gene)
            print(pos)

            if pos[2] == pos[3]:
                # Bad CDS, skip this one;
                continue

            utr5_l = 0
            utr5_r = pos[2]

            cds_l = pos[2]
            cds_r = pos[3]
            cds_len  = pos[3] - pos[2]

            utr3_l = pos[3]
            utr3_r = pos[1]
            utr3_len = (pos[1] - pos[3]) +1 # in case some utr = 0

            tlen = len(gene['loc'])

            for d in gene['doms']:
                if TE:
                    if d['dom'] not in TE:
                        continue

                s = d['span'][0]
                e = d['span'][1]

                #print(s, e, 'utr', utr5_l, utr5_r, 'cds', cds_l, cds_r, cds_len, 'utr3', utr3_l, utr3_r, utr3_len)

                if s <= utr5_r:
                    ls = max([math.floor(s / utr5_r * 1000), 0])
                    le = min([math.ceil(e / utr5_r * 1000), 1000])
                    data[gl]['utr5'][ls:le] += 1
                    #print("Add 5'UTR")
                if e >= cds_l and s <= cds_r:
                    ls = max([math.floor((s-cds_l) / cds_len * 1000), 0])
                    le = min([math.ceil((e-cds_l) / cds_len * 1000), 1000])
                    data[gl]['cds'][ls:le] += 1
                    #print('Add CDS')
                if e > utr3_l:
                    ls = max([math.floor((s-utr3_l) / utr3_len * 1000), 0])
                    le = min([math.ceil((e-utr3_l) / utr3_len * 1000), 1000])
                    data[gl]['utr3'][ls:le] += 1
                    #print("Add 3'UTR")

                #print(ls, le)
                #print()

            if (n+1) % 1000 == 0:
                print('Processed: {:,} transcripts'.format(n+1))
                if n> 20: break

    fig = plot.figure(figsize=[2.2,1.4])

    ymax = max([utr5.max(), cds.max(), utr3.max(), 1])

    ax1 = fig.add_subplot(131)
    ax1.plot(utr5)
    ax1.set_ylim([0, ymax])
    ax1.tick_params(labelsize=6)
    ax1.set_xticklabels('')

    ax2 = fig.add_subplot(132)
    ax2.plot(cds)
    ax2.set_ylim([0, ymax])
    ax2.set_yticklabels('', fontsize=6)
    ax2.set_xticklabels('')

    ax3 = fig.add_subplot(133)
    ax3.plot(utr3)
    ax3.set_ylim([0, ymax])
    ax3 .set_yticklabels('')
    ax3.set_xticklabels('')


    #fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    plot.close(fig)

    print('Saved %s' % filename)
