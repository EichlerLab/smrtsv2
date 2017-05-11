
def TSDAlign(query, target, side):
    qlen = len(query)
    tlen = len(target)
    
    score = [ [0]*(tlen+1) for i in range(qlen+1)]
    
    if (side == 'suffix'):
        query  = query[::-1]
        target = target[::-1]
    # The TSD is an exact match,
    maxScore = 0;
    maxi = 0;
    maxj = 0;

    for i in range(qlen):
        for j in range(tlen):
            if (query[i] == target[j]):
                score[i+1][j+1] += score[i][j] + 1
                if (score[i+1][j+1] >= maxScore):
                    maxScore = score[i+1][j+1]
                    maxi = i
                    maxj = j
    qs = query[maxi+1-maxScore:maxi]
    ts = target[maxj+1-maxScore:maxj]

    if (side == 'suffix'):
        qs = qs[::-1]
        ts = ts[::-1]

    return (qs, ts, maxScore)
    



def SWAlign(query, target, match=1,mismatch=-1,indel=-1):
    qlen = len(query)
    tlen = len(target)
    
    scoremat = [ [0]*(tlen+1) for i in range(qlen+1)]
    pathmat  = [ [0]*(tlen+1) for i in range(qlen+1)]

    endalign = 0
    up = 1
    left = 2
    diag = 3
    maxScore = 0
    maxI = 0
    maxJ = 0
    for i in range(qlen):
        for j in range(tlen):
            
            if query[i] == target[j]:
                op = match
            else:
                op = mismatch
            insScore = scoremat[i+1][j] + indel
            delScore = scoremat[i][j+1] + indel
            mmScore  = scoremat[i][j] + op
            
            scoremat[i+1][j+1] = max(insScore, max(delScore, max(mmScore, 0)))
            if (scoremat[i+1][j+1] == 0):
                pathmat[i+1][j+1] = endalign
            elif (scoremat[i+1][j+1] == insScore):
                pathmat[i+1][j+1] = left
            elif (scoremat[i+1][j+1] == delScore):
                pathmat[i+1][j+1] = up
            elif (scoremat[i+1][j+1] == mmScore):
                pathmat[i+1][j+1] = diag
            if (scoremat[i+1][j+1] > maxScore):
                maxScore = scoremat[i+1][j+1]
                maxI = i+1
                maxJ = j+1
    # do backgrack
    i = maxI
    j = maxJ
    optT = []
    optQ = []
    while (pathmat[i][j] != endalign and i > 0 and j > 0):
        if (pathmat[i][j] == diag):
            optT.append(target[j-1])
            optQ.append(query[i-1])
            i -= 1
            j -= 1
        elif (pathmat[i][j] == left):
            optQ.append(query[i-1])
            i -= 1
        elif (pathmat[i][j] == up):
            optT.append(target[j-1])
            j -= 1
    optT.reverse()
    optQ.reverse()
    optQs = ''.join(optQ)
    optTs = ''.join(optT)
    return (optQs, optTs, maxScore)
        
            
    
                
                
