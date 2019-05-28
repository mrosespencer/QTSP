from gurobipy import *
import time
import GetVal
import math
import PrintM


def SolveTSP(n, c, q, qname, presolve):
    t0 = time.time()
    # Create model
    m = Model()

    logname = "mtz_linb10-"+qname +"-log"

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
    f = n * (n - 1)
    uppery = {}
    lowery={}
    for i in range(f):
        uppery[i] = 0
        lowery[i] = 0
        for j in range(f):
            uppery[i] += max(q[i,j],0)
            lowery[i] += min(q[i,j],0)


    upperp = {}
    lowerp = {}

    upperv = 9
    lowerv = 0


    for i in range(f):
        if uppery[i] > 0:
            upperp[i] = math.floor(math.log(uppery[i],10))+1
        else:
            upperp[i] = 0
        if lowery[i] < 0:
            lowerp[i] = math.floor(math.log(abs(lowery[i]),10))+1
        else:
            lowerp[i] = 0

    # Create variables
    x = {}
    rone = {}
    rtwo = {}
    vone = {}
    vtwo = {}
    u = {}

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))
            if i != j:
                ij = GetVal.getval(i, j, n)
                for p in range(upperp[ij]):
                    rone[ij,p] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'r1 ' + str(i)+str(j) + '_' +str(p))
                    vone[ij,p] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'v1 ' + str(i)+str(j) + '_' +str(p))
                for o in range(lowerp[ij]):
                    rtwo[ij,o] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'r2 ' + str(i)+str(j) + '_' +str(o))
                    vtwo[ij,o] = m.addVar(vtype = GRB.INTEGER, ub = 9,lb = 0, name = 'v2 ' + str(i)+str(j) + '_' +str(o))

    for i in range(1, n):
        u[i] = m.addVar(vtype=GRB.CONTINUOUS, name='u' + str(i))

    # Set objective

    objective = LinExpr()



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
