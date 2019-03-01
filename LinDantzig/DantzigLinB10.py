from gurobipy import *
import time
import math
import PrintM

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



def SolveTSP(n, c, q, qname):
    t0 = time.time()

    # Callback - use lazy constraints to eliminate sub-tours

    def subtourelim(model, where):
        if where == GRB.callback.MIPSOL:
            # print("Added constraint")
            selected = []
            # make a list of edges selected in the solution
            for i in range(n):
                sol = model.cbGetSolution([model._vars[i, j] for j in range(n)])
                selected += [(i, j) for j in range(n) if sol[j] > 0.5]
            # find the shortest cycle in the selected edge list
            tour = subtour(selected)
            if len(tour) < n:
                # add a subtour elimination constraint

                model.cbLazy( quicksum(model._vars[i, j] for j  in tour for i in tour) <= (len(tour)-1) )


    # Given a list of edges, finds the shortest subtour

    def subtour(edges):
        visited = [False] * n
        cycles = []
        lengths = []
        selected = [[] for i in range(n)]
        for x, y in edges:
            selected[x].append(y)
        while True:
            current = visited.index(False)
            thiscycle = [current]
            while True:
                visited[current] = True
                neighbors = [x for x in selected[current] if not visited[x]]
                if len(neighbors) == 0:
                    break
                current = neighbors[0]
                thiscycle.append(current)
            cycles.append(thiscycle)
            lengths.append(len(thiscycle))
            if sum(lengths) == n:
                break
        return cycles[lengths.index(min(lengths))]


    # Create model
    m = Model()

    logname = "dantzig_linb10-"+qname +"-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0

    # Set time limit to 3 hours

    m.setParam(GRB.Param.TimeLimit, 10800.0)
    # m.setParam(GRB.Param.TimeLimit, 300.0)

    m.setParam("logfile", "%s.txt" % logname)

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

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))
            if i != j:
                ij = getval(i, j, n)
                for p in range(upperp[ij]):
                    rone[ij,p] = m.addVar(vtype=GRB.INTEGER, ub = 9, lb = 0, name = 'r1 ' + str(i)+str(j) + '_' +str(p))
                    vone[ij,p] = m.addVar(vtype = GRB.INTEGER, ub = 9, lb = 0, name = 'v1 ' + str(i)+str(j) + '_' +str(p))
                for o in range(lowerp[ij]):
                    rtwo[ij,o] = m.addVar(vtype=GRB.INTEGER, ub = 9, lb = 0, name = 'r2 ' + str(i)+str(j) + '_' +str(o))
                    vtwo[ij,o] = m.addVar(vtype = GRB.INTEGER, ub = 9,lb = 0, name = 'v2 ' + str(i)+str(j) + '_' +str(o))


    # Set objective

    objective = LinExpr()


    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])
            if i != j:
                ij = getval(i, j, n)
                for p in range(upperp[ij]):
                    objective.addTerms(10**p, rone[ij,p])
                for o in range(lowerp[ij]):
                    objective.addTerms(-(10**o), rtwo[ij, o])


    m.setObjective(objective, GRB.MINIMIZE)

    for i in range(n):
        x[i, i].ub = 0
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))


    for i in range(n):
        for j in range(n):
            if i != j:
                ij = getval(i, j, n)

                qsum = LinExpr()
                r = LinExpr()


                for k in range(n):
                    for l in range(n):
                        if k != l:
                            kl = getval(k, l, n)

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

                for o in range(lowerp[ij]):
                    m.addConstr(upperv*x[i,j] + vtwo[ij,o]-rtwo[ij,o] <= upperv, "b5-"+str(ij)+"_"+str(o))
                    m.addConstr(lowerv*x[i, j] + vtwo[ij, o] - rtwo[ij, o] >= lowerv, "b10-" + str(ij) + "_" + str(o))
                    m.addConstr(rtwo[ij,o] <= x[i,j]*upperv, "b6-"+ str(ij)+"_"+str(o))


                qsum.clear()
                r.clear()


    m.update()



    m._vars = x
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)
    gap = m.MIPGAP
    status = m.status
    # obj = m.getObjective()
    solcnt = m.SolCount



    finalx = {}

    if solcnt > 0:

        vals = m.getAttr('x', x)
        selected = tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

        tour = subtour(selected)
        assert len(tour) == n, "The tour length is: %d" %len(tour)




        varlist = []
        for v in m.getVars():
            if v.VType == GRB.BINARY:
                varlist.append(v.x)
        # selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5]



        for i in range(n):
            for j in range(n):
                finalx[i, j] = varlist[(n * i) + j]


    else:
        for i in range(n):
            for j in range(n):
                finalx[i, j] = 0


    t1 = time.time()
    totaltime = t1 - t0

    return m.objVal, totaltime, finalx, gap, status