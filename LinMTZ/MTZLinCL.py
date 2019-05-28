from gurobipy import *
import time
import GetVal


def SolveTSP(n, c, q, qname, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "mtz_lincl-"+qname +"-log"

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

    # Create variables
    x = {}
    u = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))

    for i in range(1, n):
        u[i] = m.addVar(vtype=GRB.CONTINUOUS, name='u' + str(i))

   # Create edge matrix

    # m.update()

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]
    f = len(e)

    # print(e)
    y = {}

    for i in range(f):
        for j in range(f):
            y[i,j] = m.addVar(vtype=GRB.CONTINUOUS, name="y "+ str(i) + '_'+str(j))

    # Set objective

    objective = LinExpr()


    for i in range(f):
        for j in range(f):
            objective.addTerms(q[i,j], y[i,j])

    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])
    m.setObjective(objective, GRB.MINIMIZE)

    # m.update()
    # print(objective)
    # Constraints

    for i in range(n):
        x[i, i].ub = 0
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
                for k in range(n):
                    for l in range(n):
                        if k != l:
                            ij = GetVal.getval(i,j,n)
                            kl = GetVal.getval(k,l,n)
                            m.addConstr(x[i,j] +x[k,l] <= 1+ y[ij,kl], "CL1-" +str(ij)+"-"+str(kl))
                            m.addConstr(y[ij,kl] <= x[i,j], "CL2-"+str(ij))
                            m.addConstr(y[ij, kl] <= x[k,l], "CL3-"+str(kl))

    m.addConstrs((y[i,j] >=0 for i in range(f) for j in range(f)), "CL4-"+str(i)+"-"+str(j))

    m.update()
    m.optimize()


    finalx = {}
    varlist = []

    for v in m.getVars():
        if v.VarName.find('x ') != -1:
            varlist.append(v.x)


    for i in range(n):
        for j in range(n):
            finalx[i,j] = varlist[(n*i)+j]

    gap = m.MIPGAP

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

    return m.objVal, totaltime, finalx, gap, status
