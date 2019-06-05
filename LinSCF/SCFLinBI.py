from gurobipy import *
import time
import GetVal


def SolveTSP(n, c, q, qname, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "scf_linbi-" + qname + "-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0
    # Set time limit to 3 hours

    m.setParam(GRB.Param.TimeLimit, 10800.0)
    if n <10:
        m.setParam("logfile","")
    else:
        m.setParam("logfile", "%s.txt" % logname)

    #Turn off presolve
    m.setParam("Presolve", presolve)


    # Create variables
    x = {}
    w = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='e ' + str(i) + '_' + str(j))
            w[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='flow-' + str(i) + '_' + str(j))

        # Create edge matrix

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]
    f = len(e)

    y = {}

    for i in range(f):
        for j in range(f):
            y[i, j] = m.addVar(vtype=GRB.BINARY, name="y " + str(i) + '_' + str(j))

    # Set objective

    objective = LinExpr()

    for i in range(f):
        for j in range(f):
            objective.addTerms(q[i, j], y[i, j])

    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])

    m.setObjective(objective, GRB.MINIMIZE)

    # Constraints

    for i in range(n):
        x[i, i].ub = 0
        w[i, i].ub = 0
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))

    m.addConstr(quicksum(w[0, j] for j in range(1, n)) == (n - 1), "c3-" + str(j))

    for i in range(n):
        for j in range(n):
            m.addConstr(w[i, j] <= (n - 1) * x[i, j], "c4-" + str(i) + '_' + str(j))

    for j in range(1, n):
        m.addConstr((quicksum(w[i, j] for i in range(n)) - quicksum(w[j, k] for k in range(n))) == 1, "c5=" + str(j))

    for i in range(n):
        for j in range(n):
            if i != j:
                for k in range(n):
                    for l in range(n):
                        if k != l:
                            ij = GetVal.getval(i,j,n)
                            kl = GetVal.getval(k,l,n)
                            m.addConstr(x[i,j] +x[k,l] <= 1+ y[ij,kl], "BI1-" +str(ij)+"-"+str(kl))
                            m.addConstr(-x[i, j] - x[k, l] <= -2* y[ij, kl], "BI2-" + str(ij) + "-" + str(kl))

    # m.setParam('OutputFlag', False)

    m.update()
    m.optimize()
    finalx = {}
    varlist = []

    # m.computeIIS()
    # m.write("model.ilp")
    # print('\nThe following constraint(s) cannot be satisfied:')
    # for c in m.getConstrs():
    #     if c.IISConstr:
    #         print('%s' % c.constrName)

    for v in m.getVars():
        if v.VType == GRB.BINARY:
            varlist.append(v.x)

    for i in range(n):
        for j in range(n):
            finalx[i, j] = varlist[(n * i) + j]


    gap = m.MIPGAP



    t1 = time.time()
    totaltime = t1 - t0

    status = m.status

    return m.objVal, totaltime, finalx, gap, status
