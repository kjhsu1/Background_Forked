def giza(genome, rov=2, sbv=1):
    values = [.001,.125,.25,.5,1,2,4,8,1000]
    count = 0
    v = 0   # value
    giza = genome.copy()
    for n in range(len(giza)):
        if count == rov+sbv:
            v += 1
            if v > len(values)-1: v = 0
            count = 0
        if count < rov:
            giza[n] = values[v]
            count += 1
        elif count < rov+sbv:
            count += 1
            giza[n] = 1 
    return giza