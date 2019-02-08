from gurobipy import *
import time


def SolveTSP(n, c, q, qname):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "scf_" + qname + "-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0
    # Set time limit to 3 hours

    m.setParam(GRB.Param.TimeLimit, 10800.0)

    m.setParam("logfile", "%s.txt" % logname)

    # Create variables
    x = {}
    y = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='e ' + str(i) + '_' + str(j))
            y[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='flow-' + str(i) + '_' + str(j))



        # Create edge matrix

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]
    f = len(e)

    g = {}

    for i in range(f):
        for j in range(f):
            g[i, j] = (e[i] * e[j])

    # Set objective

    objective = QuadExpr()

    for i in range(f):
        for j in range(f):
            objective.addTerms(q[i, j], e[i], e[j])

    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j], x[i, j])

    m.setObjective(objective, GRB.MINIMIZE)

    # Constraints

    for i in range(n):
        x[i, i].ub = 0
        y[i, i].ub = 0
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))


    m.addConstr(quicksum(y[0, j] for j in range(1, n)) == (n - 1), "c3-" + str(j))

    for i in range(n):
        for j in range(n):
            m.addConstr(y[i, j] <= (n - 1) * x[i, j], "c4-" + str(i) + '_' + str(j))

    for j in range(1, n):
        m.addConstr((quicksum(y[i, j] for i in range(n)) - quicksum(y[j, k] for k in range(n))) == 1, "c5=" + str(j))

    # m.setParam('OutputFlag', False)

    m.update()
    m.optimize()
    finalx = {}
    varlist = []

    for v in m.getVars():
        if v.VType != GRB.CONTINUOUS:
            varlist.append(v.x)

    for i in range(n):
        for j in range(n):
            finalx[i, j] = varlist[(n * i) + j]

    gap = m.MIPGAP



    t1 = time.time()
    totaltime = t1 - t0

    return m.objVal, totaltime, finalx, gap
