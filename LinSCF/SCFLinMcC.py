from gurobipy import *
import time
import GetVal


def SolveTSP(n, c, q, qname, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "scf_linmcc-" + qname + "-log"

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
    m.setParam("PreQLinearize", presolve)


    # Define constants
    f = n*(n-1)
    upperv = {}
    lowerv={}
    for i in range(f):
        upperv[i] = 0
        lowerv[i] = 0
        for j in range(f):
            upperv[i] += max(q[i,j],0)
            lowerv[i] += min(q[i,j],0)

    # Create variables
    x = {}
    w = {}
    t = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))
            if i != j:
                ij = GetVal.getval(i, j, n)
                t[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='t-' + str(i) + '_' + str(j), lb = lowerv[ij]-1)   # determine appropriate lower bound
            else:
                t[i,j] = m.addVar(vtype=GRB.CONTINUOUS, name='t-' + str(i) + '_' + str(j))
            w[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='flow-' + str(i) + '_' + str(j))

        # Create edge matrix

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]

    # y = {}
    #
    # for i in range(f):
    #     for j in range(f):
    #         y[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name="y " + str(i) + '_' + str(j))

    # Set objective

    objective = LinExpr()


    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])
            objective.add(t[i, j])

    m.setObjective(objective, GRB.MINIMIZE)



    # Constraints

    for i in range(n):
        x[i, i].ub = 0
        w[i, i].ub = 0
        m.addConstr(t[i, i] == 0)
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))

    m.addConstr(quicksum(w[0, j] for j in range(1, n)) == (n - 1), "c3-" + str(j))

    for i in range(n):
        for j in range(n):
            if i != j:
                m.addConstr(w[i, j] <= (n - 1) * x[i, j], "c4-" + str(i) + '_' + str(j))

    for j in range(1, n):
        m.addConstr((quicksum(w[i, j] for i in range(n)) - quicksum(w[j, k] for k in range(n))) == 1, "c5-" + str(j))


    for i in range(n):
        for j in range(n):
            if i != j:
                ij = GetVal.getval(i, j, n)
                m.addConstr(t[i,j] <= x[i, j]*upperv[ij], "MCL1-" + str(ij))
                m.addConstr(t[i,j] >= x[i, j]*lowerv[ij], "MCL2-" + str(ij))

                con1 = LinExpr()
                con2 = LinExpr()

                for k in range(n):
                    for l in range(n):
                        if k != l:
                            kl = GetVal.getval(k,l,n)

                            con1.addTerms(q[ij, kl], x[k, l])
                            con2.addTerms(q[ij, kl], x[k, l])

                con1.addTerms(lowerv[ij], x[i,j])
                con2.addTerms(upperv[ij], x[i,j])

                con1.addConstant(-lowerv[ij])
                con2.addConstant(-upperv[ij])


                m.addConstr(con1 >= t[i,j] , "MCL3-" + str(ij))
                m.addConstr(con2 <= t[i,j], "MCL4" + str(ij))

                con1.clear()
                con2.clear()



    # m.setParam('OutputFlag', False)

    m.update()
    m.optimize()
    finalx = {}
    varlist = []
    if m.status != GRB.Status.INFEASIBLE:




        for v in m.getVars():
            if v.VType == GRB.BINARY:
                varlist.append(v.x)

        for i in range(n):
            for j in range(n):
                finalx[i, j] = varlist[(n * i) + j]


        gap = m.MIPGAP



    else:
        m.computeIIS()
        m.write("model.ilp")
        print('\nThe following constraint(s) cannot be satisfied:')
        for c in m.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
        gap = 0

    status = m.status
    t1 = time.time()
    totaltime = t1 - t0
    return m.objVal, totaltime, finalx, gap, status
