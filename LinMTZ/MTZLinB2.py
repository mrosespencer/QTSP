from gurobipy import *
import time
import GetVal
import math


def SolveTSP(n, c, q, qname):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "mtz_linb2-"+qname +"-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0

    # Set time limit to 3 hours

    m.setParam(GRB.Param.TimeLimit, 10800.0)

    m.setParam("logfile", "%s.txt" % logname)

    # Define constants
    f = n * (n - 1)
    upperv = {}
    lowerv = {}
    for i in range(f):
        upperv[i] = 0
        lowerv[i] = 0
        for j in range(f):
            upperv[i] += max(q[i, j], 0)
            lowerv[i] += min(q[i, j], 0)

    upperp = {}
    lowerp = {}

    for i in range(f):
        if upperv[i] > 0:
            upperp[i] = math.floor(math.log(upperv[i], 2)) + 1
        else:
            upperp[i] = 0
        if lowerv[i] < 0:
            lowerp[i] = math.floor(math.log(abs(lowerv[i]), 2)) + 1
        else:
            lowerp[i] = 0

    # Create variables
    x = {}
    wone = {}
    wtwo = {}
    tone = {}
    ttwo = {}
    u = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))
            if i != j:
                ij = GetVal.getval(i, j, n)
                for p in range(upperp[ij]):
                    wone[ij, p] = m.addVar(vtype=GRB.BINARY, name='w1 ' + str(i) + str(j) + '_' + str(p))
                    tone[ij, p] = m.addVar(vtype=GRB.BINARY, name='t1 ' + str(i) + str(j) + '_' + str(p))
                for o in range(lowerp[ij]):
                    wtwo[ij, o] = m.addVar(vtype=GRB.BINARY, name='w1 ' + str(i) + str(j) + '_' + str(o))
                    ttwo[ij, o] = m.addVar(vtype=GRB.BINARY, name='t1 ' + str(i) + str(j) + '_' + str(o))

    for i in range(1, n):
        u[i] = m.addVar(vtype=GRB.CONTINUOUS, name='u' + str(i))

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
                    costone[ij,p] = math.pow(2,p)
                    objective.addTerms(costone[ij,p], wone[ij,p])
                for o in range(lowerp[ij]):
                    costtwo[ij,o] = -math.pow(2,o)
                    objective.addTerms(costtwo[ij,o], wtwo[ij, o])
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
                ij = GetVal.getval(i, j, n)

                qsum = LinExpr()
                r = LinExpr()


                for k in range(n):
                    for l in range(n):
                        if k != l:
                            kl = GetVal.getval(k, l, n)

                            qsum.addTerms(q[ij,kl], x[k,l])


                for p in range(upperp[ij]):
                    r.addTerms(costone[ij,p], tone[ij,p])
                for o in range(lowerp[ij]):
                    r.addTerms(costtwo[ij,o],ttwo[ij,o])

                m.addConstr(qsum == r, "b1-"+ str(ij))

                for p in range(upperp[ij]):
                    m.addConstr(x[i,j] + tone[ij,p]-wone[ij,p] <= 1, "b2-"+str(ij)+"_"+str(p))
                    m.addConstr(wone[ij,p] <= x[i,j], "b3-"+ str(ij)+"_"+str(p))
                    m.addConstr(wone[ij,p] <= tone[ij,p], "b4-" + str(ij)+"_"+str(p))
                for o in range(lowerp[ij]):
                    m.addConstr(x[i,j] + ttwo[ij,o]-wtwo[ij,o] <= 1, "b5-"+str(ij)+"_"+str(o))
                    m.addConstr(wtwo[ij,o] <= x[i,j], "b6-"+ str(ij)+"_"+str(o))
                    m.addConstr(wtwo[ij,o] <= ttwo[ij,o], "b7-" + str(ij)+"_"+str(o))

                qsum.clear()
                r.clear()

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
    # print('Obj MTZ:', m.objVal)

    # logfile = open('%s.log' % (qname), 'w')

    # logfile = m._logfile

    # m.write("%s.log" %qname)

    t1 = time.time()
    totaltime = t1 - t0

    status = m.status

    return m.objVal, totaltime, finalx, gap, status
