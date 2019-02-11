from random import *
import QMod

minsize = 5
maxsize = 31
s = False

properties = ["nonneg", "negskew", "posskew", "balanced", "psd", "rankone", "ranktwo", "nonnegpsd"]


def Makecost(n, s):
    e = n * (n - 1)

    # random.seed(0)

    q = {}

    symsize = 0
    for i in range(n - 1):
        symsize += i

    if s == False:
        for i in range(e):
            for j in range(e):
                q[i, j] = randint(1, 10)
                # q[i,j] = 1
            # q[i,i] = 0
    else:
        for i in range(symsize):
            for j in range(e):
                q[i, j] = randint(0, 10)

    return q


def MakeQ(n, s, p, t):
    name = "Q" + str(n) + str(p) + "-" + str(t)

    filename = "%s.txt" % name
    file = open(filename, 'w')
    e = n * (n - 1)

    # random.seed(0)

    q = {}

    symsize = 0
    for i in range(n - 1):
        symsize += i

    B = {}
    D = {}
    col = {}
    col1 = {}
    col2 = {}
    row = {}
    row1 = {}
    row2 = {}

    # print(p)

    for i in range(e):
        for j in range(e):
            B[i, j] = randint(-5, 5)
            D[i, j] = randint(5, 10)

        # print(i)
        col[i] = randint(-10, 10)
        col1[i] = randint(-10, 10)
        col2[i] = randint(-10, 10)
        row[i] = randint(-5, 5)
        row1[i] = randint(-5, 5)
        row2[i] = randint(-5, 5)

    if s == False:
        for i in range(e):
            for j in range(e):
                if p == "nonneg":
                    q[i, j] = int(randint(5, 10))
                if p == "negskew":
                    q[i, j] = randint(-10, 5)
                if p == "posskew":
                    q[i, j] = randint(-5, 10)
                if p == "balanced":
                    q[i, j] = randint(-5, 5)
                if p == "psd":
                    q[i, j] = (B[i, j] + B[j, i])
                if p == "rankone":
                    q[i, j] = col[i] * row[j]
                if p == "ranktwo":
                    q[i, j] = col1[i] * row1[j] + col2[i] * row2[j]
                if p == "nonnegpsd":
                    q[i, j] = (D[i, j] + D[j, i])


    else:
        for i in range(symsize):
            for j in range(e):
                q[i, j] = randint(0, 10)

    a = [[q[i, j] for i in range(e)] for j in range(e)]

    for i in range(len(a)):
        print()
        file.write("\n")
        for j in a[i]:
            # print("%f " % j, end='')
            file.write("%f " % j)

    file.close()

    return q


def MakeD(n, s):
    d = {}

    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    # d[i,j,k,l] = randint(0,10)
                    d[i, j, k, l] = 1
            d[i, j, i, j] = 0
        d[i, i, i, i] = 0

    return d


def MakeQ2(n, s):
    q2 = {}

    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    ij = QMod.getval(i, j, n)
                    kl = QMod.getval(k, l, n)

                    if i == k:
                        if j != l:
                            q2[ij, kl] = 0
                        else:
                            q2[ij, kl] = randint(1, 10)
                    else:
                        if j == l:
                            q2[ij, kl] = 0
                        else:
                            q2[ij, kl] = randint(1, 10)
            q2[ij, ij] = 0

    return q2


def MakeC(n, s):
    c = {}
    for i in range(n):
        for j in range(i, n):
            c[i, j] = randint(1, 10)
            # c[i,j] = 1
            if s == True:
                c[j, i] = c[i, j]
            else:
                c[j, i] = randint(1, 10)
                # c[j, i] = 1

        c[i, i] = 0
    return c


def MakeCmat(n, t):
    name = "C" + str(n) + "-" + str(t)

    filename = "%s.txt" % name
    file = open(filename, 'w')

    c = {}
    for i in range(n):
        for j in range(n):
            c[i, j] = randint(0, 1)

    a = [[c[i, j] for i in range(n)] for j in range(n)]

    for i in range(len(a)):
        print()
        file.write("\n")
        for j in a[i]:
            # print("%f " % j, end='')
            file.write("%f " % j)

    file.close()


for n in range(minsize, maxsize, 5):
    if n < 15:
        trials = 5
    else:
        trials = 1

    for t in range(trials):
        for p in range(7, 8):
            MakeQ(n, s, properties[p], t)
        MakeCmat(n, t)
