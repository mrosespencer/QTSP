def printmatrix(m,l):
    line = []

    for i in range(l):
        for j in range(l):
            line.append(m[i,j])
        print(line)
        line =[]