
import math, numpy, sys
import matplotlib.pyplot as plot

sys.path.append('../../')
import shared

def draw_density(filename, selected_genes, TE=None):

    arr = numpy.zeros(1000) # scaled bins to put in;
    for n, gene in enumerate(selected_genes):
        # scale the TE to the mRNA
        #print(gene)

        tlen = len(gene['loc'])

        for d in gene['doms']:
            if TE:
                if d['dom'] not in TE:
                    continue

            s = math.floor(d['span'][0] / tlen * 1000)
            e = math.ceil(d['span'][1] / tlen * 1000)

            #print(s,e)

            arr[s:e] += 1

    fig = plot.figure()
    ax = fig.add_subplot(111)

    ax.plot(arr)

    fig.savefig(filename)
    plot.close(fig)
    print('Saved %s' % filename)

def draw_density_utrs(filename, selected_genes, TE=None):

    utr5 = numpy.zeros(1000) # scaled bins to put in;
    cds = numpy.zeros(1000)
    utr3 = numpy.zeros(1000)

    for n, gene in enumerate(selected_genes):
        # scale the TE to the mRNA
        pos = shared.convert_genocode_to_local(gene)

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

            #print(s, e, utr5_l, utr5_r, cds_l, cds_r, cds_len, utr3_l, utr3_r, utr3_len)

            if s <= utr5_r:
                ls = max(math.floor(s / utr5_r * 1000), 0)
                le = min(math.ceil(e / utr5_r * 1000), 1000)
                utr5[ls:le] += 1
                #print("Add 5'UTR")
            if e >= cds_l and s <= cds_r:
                ls = max(math.floor((s-cds_l) / cds_len * 1000), 0)
                le = min(math.ceil((e-cds_l) / cds_len * 1000), 1000)
                cds[s:e] += 1
                #print('Add CDS')
            if e > utr3_l:
                ls = max(math.floor((s-utr3_l) / utr3_len * 1000), 0)
                le = min(math.ceil((e-utr3_l) / utr3_len * 1000), 1000)
                utr3[s:e] += 1
                #print("Add 3'UTR")

        if (n+1) % 1000 == 0:
            print('Processed: {:,} transcripts'.format(n+1))
            #if n> 2000: break

    fig = plot.figure(figsize=[7,5])

    ymax = max([utr5.max(), cds.max(), utr3.max(), 1])

    ax1 = fig.add_subplot(131)
    ax1.plot(utr5)
    ax1.set_ylim([0, ymax])

    ax2 = fig.add_subplot(132)
    ax2.plot(cds)
    ax2.set_ylim([0, ymax])

    ax3 = fig.add_subplot(133)
    ax3.plot(utr3)
    ax3.set_ylim([0, ymax])


    fig.savefig(filename)
    plot.close(fig)
    print('Saved %s' % filename)
