from gurobipy import *
import time
import GetVal
import PrintM


def SolveTSP(n, c, q, qname, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "mtz_linmcc-"+qname +"-log"

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


    # Define constants
    f = n*(n-1)

    upperv = {}
    lowerv = {}

    for i in range(f):
        upperv[i] = 0
        lowerv[i] = 0

        for j in range(f):
            upperv[i] += max(q[i,j],0)
            lowerv[i] += min(q[i,j],0)

    # Create variables
    x = {}
    u = {}
    t = {}

    for i in range(n):
        for j in range(n):
            # x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))
            x[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='x ' + str(i) + '_' + str(j), lb=0, ub=1)
            if i != j:
                ij = GetVal.getval(i, j, n)
                t[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='t-' + str(i) + '_' + str(j), lb = lowerv[ij]-1)   # determine appropriate lower bound
            else:
                t[i,j] = m.addVar(vtype=GRB.CONTINUOUS, name='t-' + str(i) + '_' + str(j))

    for i in range(1, n):
        u[i] = m.addVar(vtype=GRB.CONTINUOUS, name='u' + str(i))

   # Create edge matrix

    # m.update()

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]


    # PrintM.printmatrix(q,f)


    # Set objective

    objective = LinExpr()


    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])
            objective.add(t[i, j])

    # print(objective.size())

    m.setObjective(objective, GRB.MINIMIZE)



    # Constraints

    for i in range(n):
        x[i, i].ub = 0
        m.addConstr(t[i, i] == 0)
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))

    for i in range(1, n):
        m.addConstr(1 <= u[i], "u-" + str(i))
        m.addConstr(u[i] <= (n - 1), "u2-" + str(i))

    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                m.addConstr(((u[i] - u[j] + ((n - 1) * x[i, j])) <= (n - 2)), "u3-" + str(i) + "_" + str(j))



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

    m.update()
    m.optimize()



    finalx = {}
    varlist = []

    if m.status != GRB.Status.INFEASIBLE:
        print("McC: %f", m.objVal)
        for v in m.getVars():
            if v.VarName.find('x ') != -1:
                varlist.append(v.x)
                print(v.varName, v.x)

        for i in range(n):
            for j in range(n):
                finalx[i, j] = varlist[(n * i) + j]

        # print(PrintM.printmatrix(finalx, n))

        # gap = m.MIPGAP
        gap = 0
        obj= m.objVal

        # for v in m.getVars():
        #     if v.x > 0:
        #         print(v.varName, v.x)
        # print('Obj:', obj)



    else:
        m.computeIIS()
        m.write("model.ilp")
        print('\nThe following constraint(s) cannot be satisfied:')
        for c in m.getConstrs():
            if c.IISConstr:
                print('%s' % c.constrName)
        gap = 0
        obj = 0

    # gap = m.MIPGAP

    # for v in m.getVars():
    #     if v.x > 0:
    #         print(v.varName, v.x)
    # print('Obj:', m.objVal)

    # logfile = open('%s.log' % (qname), 'w')

    # logfile = m._logfile

    # m.write("%s.log" %qname)

    t1 = time.time()
    totaltime = t1 - t0

    status = m.status


    return obj, totaltime, finalx, gap, status
