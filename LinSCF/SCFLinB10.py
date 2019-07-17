from gurobipy import *
import time
import GetVal
import math


def SolveTSP(n, c, q, qname, presolve, relax):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "scf_linb10-" + qname + "-log"

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
    f = n * (n - 1)
    uppery = {}
    lowery = {}
    for i in range(f):
        uppery[i] = 0
        lowery[i] = 0
        for j in range(f):
            uppery[i] += max(q[i, j], 0)
            lowery[i] += min(q[i, j], 0)

    upperp = {}
    lowerp = {}

    upperv = 9
    lowerv = 0

    for i in range(f):
        if uppery[i] > 0:
            upperp[i] = math.floor(math.log(uppery[i], 10)) + 1
        else:
            upperp[i] = 0
        if lowery[i] < 0:
            lowerp[i] = math.floor(math.log(abs(lowery[i]), 10)) + 1
        else:
            lowerp[i] = 0

    # Create variables
    x = {}
    w = {}
    rone = {}
    rtwo = {}
    vone = {}
    vtwo = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))

            w[i, j] = m.addVar(vtype=GRB.CONTINUOUS, name='flow-' + str(i) + '_' + str(j))
            if i != j:
                ij = GetVal.getval(i, j, n)
                for p in range(upperp[ij]):
                    rone[ij,p] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'r1 ' + str(i)+str(j) + '_' +str(p))
                    vone[ij,p] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'v1 ' + str(i)+str(j) + '_' +str(p))

                for o in range(lowerp[ij]):
                    rtwo[ij, o] = m.addVar(vtype=GRB.INTEGER, ub=9, lb=0, name='r2 ' + str(i) + str(j) + '_' + str(o))
                    vtwo[ij, o] = m.addVar(vtype=GRB.INTEGER, ub=9, lb=0, name='v2 ' + str(i) + str(j) + '_' + str(o))


    # Set objective

    objective = LinExpr()

    costone = {}
    costtwo = {}

    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])
            if i != j:
                ij = GetVal.getval(i, j, n)
                for p in range(upperp[ij]):
                    objective.addTerms(10**p, rone[ij,p])
                for o in range(lowerp[ij]):
                    objective.addTerms(-(10**o), rtwo[ij, o])

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
                ij = GetVal.getval(i, j, n)

                qsum = LinExpr()
                r = LinExpr()


                for k in range(n):
                    for l in range(n):
                        if k != l:
                            kl = GetVal.getval(k, l, n)

                            qsum.addTerms(q[ij,kl], x[k,l])


                for p in range(upperp[ij]):
                    r.addTerms(10**p, vone[ij,p])
                for o in range(lowerp[ij]):
                    r.addTerms(-(10**o),vtwo[ij,o])

                m.addConstr(qsum == r, "b1-"+ str(ij))

                for p in range(upperp[ij]):
                    m.addConstr(upperv*x[i,j] + vone[ij,p]-rone[ij,p] <= upperv, "b2-"+str(ij)+"_"+str(p))
                    m.addConstr(lowerv*x[i, j] + vone[ij, p] - rone[ij, p] >= lowerv, "b9-" + str(ij) + "_" + str(p))
                    m.addConstr(rone[ij,p] <= x[i,j]*upperv, "b3-"+ str(ij)+"_"+str(p))
                    # m.addConstr(rone[ij, p] >= x[i, j] * lowerv, "b8-" + str(ij) + "_" + str(p))
                for o in range(lowerp[ij]):
                    m.addConstr(upperv*x[i,j] + vtwo[ij,o]-rtwo[ij,o] <= upperv, "b5-"+str(ij)+"_"+str(o))
                    m.addConstr(lowerv*x[i, j] + vtwo[ij, o] - rtwo[ij, o] >= lowerv, "b10-" + str(ij) + "_" + str(o))
                    m.addConstr(rtwo[ij,o] <= x[i,j]*upperv, "b6-"+ str(ij)+"_"+str(o))
                    # m.addConstr(rtwo[ij, o] >= x[i, j] * lowerv, "b11-" + str(ij) + "_" + str(o))

                qsum.clear()
                r.clear()

    # m.setParam('OutputFlag', False)

    m.update()
    if relax == True:
        m = m.relax()
        gap = 0
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
        if v.VarName.find('x ') != -1:
            varlist.append(v.x)

    for i in range(n):
        for j in range(n):
            finalx[i, j] = varlist[(n * i) + j]


    if relax == False:
        gap = m.MIPGAP

    t1 = time.time()
    totaltime = t1 - t0

    status = m.status

    return m.objVal, totaltime, finalx, gap, status
