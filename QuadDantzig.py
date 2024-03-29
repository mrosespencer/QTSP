from gurobipy import *
import time
import GetVal


def SolveTSP(n, c, q, qname, adj, presolve):
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

                model.cbLazy(quicksum(model._vars[i, j] for j in tour for i in tour) <= (len(tour) - 1))

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

    logname = "Dantzig_" + qname + "-log"

    # m.setParam('OutputFlag', False)
    m.Params.logtoconsole = 0.0
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

    for i in range(n):
        for j in range(n):
            x[i, j] = m.addVar(vtype=GRB.BINARY, name='x ' + str(i) + '_' + str(j))

    e = []
    e = [x[i, j] for i in range(n) for j in range(n) if i != j]
    f = len(e)

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
                            objective.addTerms(q[ij, jk], x[i, j], x[j,k])
    for i in range(n):
        for j in range(n):
            objective.addTerms(c[i, j], x[i, j])

    m.setObjective(objective, GRB.MINIMIZE)

    for i in range(n):
        x[i, i].ub = 0
        m.addConstr((quicksum(x[i, j] for j in range(n))) == 1, "c1-" + str(i))
        m.addConstr((quicksum(x[j, i] for j in range(n))) == 1, "c2-" + str(i))

    m.update()
    # m.optimize()

    m._vars = x
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)

    vals = m.getAttr('x', x)
    selected = tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

    tour = subtour(selected)
    assert len(tour) == n, "The tour length is: %d" % len(tour)

    gap = m.MIPGAP

    varlist = []
    for v in m.getVars():
        varlist.append(v.x)
    # selected = [(i,j) for i in range(n) for j in range(n) if solution[i,j] > 0.5]

    finalx = {}

    for i in range(n):
        for j in range(n):
            finalx[i, j] = varlist[(n * i) + j]

    t1 = time.time()
    totaltime = t1 - t0
    status = m.status
    return m.objVal, totaltime, finalx, gap, status
