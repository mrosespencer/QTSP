
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