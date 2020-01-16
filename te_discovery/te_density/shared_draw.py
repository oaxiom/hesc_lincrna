
import math, numpy, sys
import matplotlib.pyplot as plot
import matplotlib.cm as cm

sys.path.append('../../')
import shared

def draw_density(filename, data_dict, TE=None, threshold=5):

    res = {k: numpy.zeros(1000) for k in data_dict} # scaled bins to put in;
    for k in data_dict:
        for gene in data_dict[k]:
            tlen = shared.get_transcript_length(gene)

            for d in gene['doms']:
                if TE:
                    if d['dom'] not in TE:
                        continue

                s = math.floor(d['span'][0] / tlen * 1000)
                e = math.ceil(d['span'][1] / tlen * 1000)

                res[k][s:e] += 1

    if max([max(i) for i in res.values()]) < threshold:
        return None # stop from drawing

    # # norm values to the number of transcripts in total set;
    for k in res:
        res[k] /= len(data_dict[k])

    fig = plot.figure(figsize=[2.8,1.3])
    fig.subplots_adjust(left=0.15, right=0.65, bottom=0.35,)
    ax = fig.add_subplot(111)

    for k in res:
        ax.plot(res[k], label=k)

    ax.tick_params(labelsize=6)
    ax.set_xticks([0, 1000])
    ax.set_xticklabels(['TSS', 'TTS'])
    #fig.savefig(filename)
    ax.legend()
    plot.legend(loc='upper left', bbox_to_anchor=(1.1, 0.8), prop={'size': 6})
    fig.savefig(filename.replace('.png', '.pdf'))
    plot.close(fig)

    print('Saved %s' % filename)
    return res

def draw_heatmap(filename, res, dataset):
    # should be just one set of labels...

    te_labels = sorted(res.keys())

    # normalise each row, for the max s and max mean in any dataset
    for te in te_labels:
        s = max([res[te][k].std() for k in dataset])
        m = max([res[te][k].mean() for k in dataset])
        for k in dataset.keys():
            res[te][k] -= m
            res[te][k] /= s

    number_of_heatmaps = len(dataset.keys())
    fig = plot.figure(figsize=[1+((number_of_heatmaps-1)*3),9])
    heat_hei = 0.010*len(te_labels)
    fig.subplots_adjust(left=0.2, right=0.80, bottom=0.95-heat_hei, top=0.95, wspace=0.1)

    for i, k in enumerate(dataset.keys()):
        res_tab = numpy.array([res[te][k] for te in te_labels])

        ax = fig.add_subplot(1,number_of_heatmaps,i+1)
        #ax.set_position([0.5, 0.1, 0.48, 0.8])

        ax.imshow(res_tab, cmap=cm.PuOr_r, aspect="auto",
            origin='lower', extent=[0, res_tab.shape[1], 0, res_tab.shape[0]])

        ax.set_yticklabels(te_labels)
        if i>0:
            ax.set_yticklabels('')
        ax.set_yticks(numpy.arange(len(te_labels))+0.5)
        ax.set_title(k, {'fontsize': 6})
        ax.set_xticklabels('')
        ax.tick_params(labelsize=6, bottom=False)
        #plot.colorbar(cax=ax)

    #fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))

def draw_density_utrs(filename, dataset_dict, TE=None):

    data = {k: {'utr5': numpy.zeros(1000),'cds': numpy.zeros(1000), 'utr3': numpy.zeros(1000)} for k in dataset_dict}

    for glk in dataset_dict:
        for n, gene in enumerate(dataset_dict[glk]):
            # scale the TE to the mRNA
            pos = gene['cds_local_locs'] # return 0, tlength, cdsl, cdsr, splice_sites

            if pos[0] == pos[1]:
                # Bad CDS, skip this one;
                continue

            #if 'tlength' not in gene:
            #    gene['tlength'] = shared.get_transcript_length(gene)

            utr5_l = 0
            utr5_r = pos[0]

            cds_l = pos[0]
            cds_r = pos[1]
            cds_len  = pos[1] - pos[0]

            utr3_l = pos[1]
            utr3_r = gene['tlength']
            utr3_len = (utr3_r - utr3_l) # in case some utr = 0
            #print(utr5_l, utr5_r, pos, utr3_l, utr3_r, utr3_len)

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
                    data[glk]['utr5'][ls:le] += 1
                    #print("Add 5'UTR")
                if e >= cds_l and s <= cds_r:
                    ls = max([math.floor((s-cds_l) / cds_len * 1000), 0])
                    le = min([math.ceil((e-cds_l) / cds_len * 1000), 1000])
                    data[glk]['cds'][ls:le] += 1
                    #print('Add CDS')
                if utr3_len > 1 and e > utr3_l: # there are a bunch of messages with UTR3' = 1
                    ls = max([math.floor((s-utr3_l) / utr3_len * 1000), 0])
                    le = min([math.ceil((e-utr3_l) / utr3_len * 1000), 1000])
                    data[glk]['utr3'][ls:le] += 1
                    #print("Add 3'UTR")

                #print(ls, le)
                #print()

            if (n+1) % 10000 == 0:
                print('Processed: {:,} transcripts'.format(n+1))

    fig = plot.figure(figsize=[2.8,1.3])
    fig.subplots_adjust(left=0.15, right=0.65, bottom=0.35,)

    # # norm values to the number of transcripts in total set;
    for k in data:
        data[k]['utr5'] /= len(dataset_dict[k])
        data[k]['cds'] /= len(dataset_dict[k])
        data[k]['utr3'] /= len(dataset_dict[k])

    ymax = max([max(data[k]['utr5'].max(), data[k]['cds'].max(), data[k]['utr3'].max()) for k in data])

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    for glk in dataset_dict:
        ax1.plot(data[glk]['utr5'], label=glk)
        ax2.plot(data[glk]['cds'], label=glk)
        ax3.plot(data[glk]['utr3'], label=glk)

    ax1.set_ylim([0, ymax])
    ax1.tick_params(labelsize=6)
    ax1.set_xticklabels('')

    ax2.set_ylim([0, ymax])
    ax2.set_yticklabels('', fontsize=6)
    ax2.set_xticklabels('')

    ax3.set_ylim([0, ymax])
    ax3.set_yticklabels('')
    ax3.set_xticklabels('')

    ax3.legend()
    plot.legend(loc='upper left', bbox_to_anchor=(1.1, 0.8), prop={'size': 6})
    #fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.pdf'))
    plot.close(fig)

    print('Saved %s' % filename)
