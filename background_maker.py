def giza(genome, rov=20, sbv=100):
    values = [.001,.125,.25,.5,0,2,4,8,1000]
    count = 0
    v = 0   # value
    for n in range(len(genome)):
        if count == rov+sbv:
            v += 1
            if v > len(values)-1: v = 0
            count = 0
        if count < rov:
            genome[n] = values[v]
            count += 1
        elif count < rov+sbv:
            count += 1 
    return genome