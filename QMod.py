

# Make Q symmetric
def half(q, e):
    qone = {}

    for i in range(e):
        for j in range(e):
            qone[i, j] = (q[i, j] + q[j, i]) * 0.5

    return qone


# Make Q into an upper triangular matrix
def triangular(q, e):
    qtwo = {}

    for i in range(e):
        for j in range(i + 1, e):
            qtwo[i, j] = q[i, j] + q[j, i]
            qtwo[j, i] = 0
        qtwo[i, i] = q[i, i]

    return qtwo


# Make Q positive semi-definite
def plusm(q, e, m):
    qthree = {}

    for i in range(e):
        for j in range(e):
            qthree[i, j] = q[i, j]
        qthree[i, i] = q[i, i] + m

    return qthree


# Make Q negative semi-definite
def minusm(q, e, m):
    qfour = {}

    for i in range(e):
        for j in range(e):
            qfour[i, j] = q[i, j]
        qfour[i, i] = (q[i, i] - m)

    return qfour



def getval(a, b, n):
    val = 0

    if a == (n):
        val = (n - 1) * (n - 1) + b

    elif b == (n):
        val = a * (n - 1) + (n - 2)

    elif a > b:
        val = a * (n - 1) + b

    elif b > a:
        val = a * (n - 1) + (b - 1)

    return val


def quadred2(q, e, n):
    qfive2 = {}

    for i in range(e):
        for j in range(e):
            qfive2[i, j] = 0

    # Creating q~

    for i in range(n - 1):
        for j in range(n - 1):
            if i !=j:
                for k in range(n - 1):
                    for l in range(n - 1):
                        if k!= l:
                            ij = getval(i, j, n)
                            kl = getval(k, l, n)
                            kn = getval(k, n, n)
                            nl = getval(n, l, n)
                            inn = getval(i, n, n)
                            nj = getval(n, j, n)



                            qfive2[ij, kl] = (
                                q[ij, kl] - q[ij, kn] - q[ij, nl] - q[inn, kl] + q[inn, kn] + q[inn, nl] - q[nj, kl] +
                                q[nj, kn] + q[nj, nl])
                    # print(qfive2[ij,kl])



    # Creating d

    d = {}
    for i in range(e):
        for j in range(e):
            d[i, j] = 0

    for i in range(n):
        for j in range(n):
            ij = getval(i, j, n)
            d[ij, ij] = qfive2[ij, ij]

    for i in range(e):
        for j in range(e):
            qfive2[i,j] = qfive2[i,j] - d[i,j]


    # Creating Q^R
    qr = {}

    for i in range(e):
        for j in range(e):
            qr[i, j] = qfive2[i, j]

    # Creating l~
    lt = {}
    for i in range(n):
        for j in range(n):
            s = 0
            if i != j:
                if i == (n - 1):
                    for k in range(n - 1):
                        ij = getval(i, j, n)
                        kn = getval(k, n, n)
                        nk = getval(n, k, n)
                        nj = getval(n, j, n)

                        s += (q[ij, kn] + q[ij, nk] + q[kn, ij] - 0 - q[kn, nj] + q[nk, ij] - 0 - q[nk, nj])

                elif j == (n - 1):
                    for k in range(n - 1):
                        ij = getval(i, j, n)
                        kn = getval(k, n, n)
                        nk = getval(n, k, n)
                        inn = getval(i, n, n)

                        s += (q[ij, kn] + q[ij, nk] + q[kn, ij] - q[kn, inn] - 0 + q[nk, ij] - q[nk, inn] - 0)

                else:

                    for k in range(n - 1):
                        ij = getval(i, j, n)
                        kn = getval(k, n, n)
                        nk = getval(n, k, n)
                        inn = getval(i, n, n)
                        nj = getval(n, j, n)

                        s += (q[ij, kn] + q[ij, nk] + q[kn, ij] - q[kn, inn] - q[kn, nj] + q[nk, ij] - q[nk, inn] - q[nk, nj])
            lt[i, j] = s

    # Creating L

    l = {}

    for i in range(n):
        for j in range(n):
            l[i, j] = 0
            if i != j:
                ij = getval(i, j, n)
                l[i, j] = lt[i, j] + d[ij, ij]

    return qr, l
