from gurobipy import *
import time
import GetVal


def SolveTSP(n, c, q, qname, adj, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "scf_" + qname + "-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0
    # Set time limit to 3 hours

    m.setParam(GRB.Param.TimeLimit, 10800.0)

    m.setParam("logfile", "%s.txt" % logname)


    #Turn off presolve
    m.setParam("Presolve", presolve)
    m.setParam("PreQLinearize", presolve)

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


    # Set objective

    objective = QuadExpr()
    if adj == False:

        for i in range(n):
            for j in range(n):
                if i != j:
                    for k in range(n):
                        for l in range(n):
                            if k != l:
                                ij = GetVal.getval(i, j, n)
                                kl = GetVal.getval(k, l, n)
                                objective.addTerms(q[ij, kl], x[i, j], x[k, l])

    else:
        for i in range(n):
            for j in range(n):
                if i != j:
                    for k in range(n):
                        if j != k:
                            ij = GetVal.getval(i, j, n)
                            jk = GetVal.getval(j, k, n)
                            objective.addTerms(q[ij, jk], x[i, j], x[j, k])
    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])

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

    status = m.status
    return m.objVal, totaltime, finalx, gap, status
