def check(x, n):
    tour = []
    tour.append(0)

    end = 1
    stop = 0

    for j in range(n):
        if x[0, j] > 0:
            tour.append(j)
            stop = j

    for l in range(2*n):
        # while end >0:
        for i in range(n):
            if x[stop, i] > 0:
                tour.append(i)
                stop = i

            if stop < 1:
                end = 0
                break


    full = False

    if len(tour) > n:
        full = True


    return tour, full
