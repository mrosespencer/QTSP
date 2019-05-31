# This program runs the experimentation for the Dantzig formulation of the linearized quadratic traveling salesman problem.
# We test one version of the quadratic cost matrix (balanced Q), and one modifcation (upper triangular). We then test
# methods to linearize the quadratic objective function. Those linearizations are: replacing the quadratic term with a
# binary variable, the classical linearization using continuous variables, the McCormick envelopes, base 2 and base 10
# linearizations. All linearized models are found in the sub folder LinDantzig. We test these linearizations against the
# non-linear form. The outputs for this program are the latex file that produces the tables in my report, as well as
#  simple text files of the objective value, gap value, Gurobi status, and the tour of the TSP. It also outputs the
#  Gurobi log files.



import VerifyTour
import QMod

import GetVal
import os

import QuadDantzig
import LinDantzig.DantzigLinCL
import LinDantzig.DantzigLinBI
import LinDantzig.DantzigLinMcC
import LinDantzig.DantzigLinB2
import LinDantzig.DantzigLinB10
#
# import QuadSCF
# import LinSCF.SCFLinCL
# import LinSCF.SCFLinMcC
# import LinSCF.SCFLinBI
# import LinSCF.SCFLinB2
# import LinSCF.SCFLinB10

# import QuadMTZ
# import LinMTZ.MTZLinCL
# import LinMTZ.MTZLinBI
# import  LinMTZ.MTZLinMcC
# import LinMTZ.MTZLinB2
# import LinMTZ.MTZLinB10


minsize = 5
maxsize = 11


# Presolve value: 0 = off, -1 = default
presolve = 0

# Used for testing purposes only
skip = False
adj = False
s = False

# Set number of trials of 5 to average (can do up to 100 size 5, and 5 size 10)
fivetrials = 100
tentrials = 5


# p = 3 # we use only balanced Q for this experiment
m = 10000
for p in range(8):

    properties = ["nonneg", "negskew", "posskew", "balanced", "psd", "rankone", "ranktwo", "nonnegpsd",  "other"]

    prop = properties[p]

    filename = "DantzigLin%s" % prop
    file = open(filename+".txt", 'w')

    file.write("\\documentclass[11pt]{article}\n")
    file.write(
        "\\usepackage{amsmath, amssymb, amsthm, amsfonts,multirow,booktabs,siunitx, lscape, multirow, rotating, booktabs}\n")
    file.write("\\begin{document}\n")
    file.write("\\begin{table}[h] \n")
    file.write("\\centering\n \def\\arraystretch{1.3}\n")
    file.write("\\resizebox{\\textwidth}{!}{\n")
    file.write("\\begin{tabular}{ @{} ccccc @{}} \\toprule\n")
    file.write(
        "Size & Quadratic & Binary & Classic & McCormick & Base 2 & Base 10 \\\\ \\midrule \n")

    objfilename = "DantzigLin%s-obj.txt" % prop
    timefilename = "DantzigLin%s-time.txt" % prop
    tourfilename = "DantzigLin%s-tour.txt" % prop
    gapfilename = "DantzigLin%s-gap.txt" % prop
    statusfilename = "DantzigLin%s-status.txt" % prop

    objfile = open(objfilename, "w")
    timefile = open(timefilename, "w")
    tourfile = open(tourfilename, "w")
    gapfile = open(gapfilename, "w")
    statusfile = open(statusfilename, "w")

    obj = {}
    time = {}
    tour = {}
    full = {}
    gap = {}
    status = {}

    count = 0

    for n in range(minsize, maxsize, 5):

        if p == 8:
            trials = 1
        elif n == 5:
            trials = fivetrials
        elif n == 10:
            trials = tentrials
        else:
            trials = 1

        for t in range(trials):

            e = n * (n - 1)

            qname = "Q" + str(n) + str(prop) + "-" + str(t)

            qcname = "%s.txt" % qname

            f = open(os.path.join('Cost', qcname), "r")

            # Read quadratic cost matrix from file
            arr = []
            for line in f:
                line = line.split()

                if line:
                    line = [int(float(i)) for i in line]
                    arr.append(line)

            f.close()

            # Convert to array
            q = {}
            for i in range(e):
                for j in range(e):
                    q[i, j] = arr[i][j]

            name = "C" + str(n) + "-" + str(t)


            cname = "%s.txt" % name
            fc = open(os.path.join('Cost', cname), "r")



            # Read linear cost from file

            arrc = []
            for line in fc:
                line = line.split()

                if line:
                    line = [int(float(i)) for i in line]
                    arrc.append(line)

            fc.close()

            c = {}
            for i in range(n):
                for j in range(n):
                    c[i, j] = arrc[i][j]

        # We choose to add the quadratic cost associated with a single edge to the linear cost so that the diagonals are 0

            for i in range(n):
                for j in range(n):
                    ij = GetVal.getval(i,j,n)
                    c[i, j] = c[i, j] + q[ij,ij]
                    q[ij,ij] = 0




            # Triangular Q
            name = qname + "-" + str(count)
            print(count)

            # Dantzig Formulation
            obj[0, count], time[0, count], x, gap[0, count], status[0, count] = QuadDantzig.SolveTSP(n, c, q, name, adj, presolve)
            tour[0, count], full[0, count] = VerifyTour.check(x, n)

            # Binary replacement
            obj[1, count], time[1, count], x, gap[1, count], status[1, count] = LinDantzig.DantzigLinBI.SolveTSP(n, c, q,name, presolve)
            tour[1, count], full[1, count] = VerifyTour.check(x, n)

            # Classic
            obj[2, count], time[2, count], x, gap[2, count], status[2, count] = LinDantzig.DantzigLinCL.SolveTSP(n, c, q, name, presolve)
            tour[2, count], full[2, count] = VerifyTour.check(x, n)

            # McCormick Envelopes
            obj[3, count], time[3, count], x, gap[3, count], status[3, count] = LinDantzig.DantzigLinMcC.SolveTSP(n, c, q, name, presolve)
            tour[3, count], full[3, count] = VerifyTour.check(x, n)

            # Base 2 Formulation
            obj[4, count], time[4, count], x, gap[4, count], status[4, count] = LinDantzig.DantzigLinB2.SolveTSP(n, c, q, name, presolve)
            tour[4, count], full[4, count] = VerifyTour.check(x, n)

            # Base 10 Formulation
            obj[5, count], time[5, count], x, gap[5, count], status[5, count] = LinDantzig.DantzigLinB10.SolveTSP(n, c, q, name, presolve)
            tour[5, count], full[5, count] = VerifyTour.check(x, n)


            objline = []
            timeline = []
            tours = []
            gapline = []
            statusline = []
            for i in range(6):
                objline.append(obj[i, count])
                a = time[i, count]
                time[i, count] = round(a, 3)
                timeline.append(time[i, count])
                tours.append(full[i, count])
                tourfile.write(str(tour[i, count]) + "\n")
                gapline.append(gap[i, count])
                statusline.append(status[i,count])
            print(objline)
            print(timeline)
            print(tours)
            print(statusline)

            objfile.write(str(objline) + "\n")
            timefile.write(str(timeline) + "\n")
            tourfile.write(str(tours) + "\n")
            gapfile.write(str(gapline) + "\n")
            statusfile.write(str(statusline) + "\n")

            count += 1
            # print(str(objline))
    fiveavg = [0] * 6
    for i in range(6):
        fiveavg[i] = round(sum(time[i, j] for j in range(fivetrials)) / fivetrials, 4)

    print(fiveavg)

    file.write("%d & %g & %g & %g & %g & %g & %g  \\\\  \n" % (
        5, fiveavg[0], fiveavg[1], fiveavg[2], fiveavg[3], fiveavg[4], fiveavg[5]))

    for i in range(fivetrials, fivetrials+tentrials):
        file.write("%d & %g & %g & %g & %g & %g & %g \\\\  \n" % (
            10, time[0, i], time[1, i], time[2, i], time[3, i], time[4, i], time[5, i]))


    file.write("\\end{tabular}\n } \n")
    file.write("\\end{table}\n")
    file.write("\\end{document}")

    file.close()
    objfile.close()
    timefile.close()
    tourfile.close()
    gapfile.close()
