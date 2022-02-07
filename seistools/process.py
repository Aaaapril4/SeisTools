def get_factor(n):
    '''
    Find factors of n
    Return:
        a list of factors
    '''
    tempn = n
    f = 2
    factor = []
    while tempn > 16:
        if tempn%f == 0 and f <= 16:
            factor.append(int(f))
            tempn /= f
            f = f+1
        else:
            f = f+1
    factor.append(int(tempn))
    return factor



def downsample(tr, sr=1):
    '''
    Downsample a trace
    Return:
        trace
    '''
    freq = round(tr.stats.sampling_rate,0)

    if freq == sr:
        return tr

    elif freq > 16*sr:
        factor = get_factor(freq/sr)
        for f in factor:
            tr.decimate(f,strict_length=False)

    elif freq > sr and freq <= 16*sr:
        tr.decimate(freq/sr,strict_length=False)
    return tr
