

def pretify_sample_names(sams):
    sams = [i.replace('_rp1', '').replace('_', ' ').replace('Hs ', '') for i in sams]
    sams = [i.replace('rnaseq ', '').replace('hesc', 'hESC').replace('esc', 'hESC').replace('ipsc', 'iPSC') for i in sams]
    sams = [i.replace(' day', ' Troph day').replace('hde d', 'DE day').replace('esc', 'hESC').replace('ipsc', 'iPSC') for i in sams]
    sams = [i.replace('embryo te', 'Trophoblast').replace('embryo ', '').replace('2c', '2C').replace('4c', '4C').replace('8c', '8C') for i in sams]
    sams = [i.replace('blastocyst', 'Blastocyst').replace('icm', 'ICM').replace('tsc', 'TSC').replace('epsc', 'EPSC') for i in sams]
    sams = [i.replace('zygote', 'Zygote').replace('oocyte', 'Oocyte').replace('morula', 'Morula').replace('unknown ', '') for i in sams]
    return(sams)

def get_cols(cond_names):
    cols = []
    for s in cond_names:
        if '2C' in s or '4C' in s or '8C' in s or 'ICM' in s or 'Zygote' in s or 'Oocyte' in s or 'Blastocyst' in s or 'Trophoblast' in s:
            cols.append('blue')
        elif 'TSC' in s:
            cols.append('purple')
        elif 'de' in s.lower(): # Primed cells;
            cols.append('orange')
        elif 'tsc' in s.lower(): # Primed cells;
            cols.append('orange')
        elif 'primed' in s:
            cols.append('green')
        elif 'naive' in s:
            cols.append('red')
        else:
            print("No colour found for sample '%s'" % s)
            cols.append('#cccccc')
    return(cols)
