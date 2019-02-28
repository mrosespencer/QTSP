# This program runs the experimentation for the MTZ formulation of the linearized quadratic traveling salesman problem.
# We test one version of the quadratic cost matrix (balanced Q), and one modifcation (upper triangular). We then test
# methods to linearize the quadratic objective function. Those linearizations are: replacing the quadratic term with a
# binary variable, the classical linearization using continuous variables, the McCormick envelopes, base 2 and base 10
# linearizations. All linearized models are found in the sub folder LiinMTZ. We test these linearization against the
# non-linear form. The outputs for this program are the latex file that produces the tables in my report, as well as
#  simple text files of the objective value, gap value, Gurobi status, and the tour of the TSP. It also outputs the
#  Gurobi log files.



import VerifyTour
import QMod

import GetVal
import os

# import QuadDantzig
# import DantzigLinCL
# import DantzigLinBI
# import DantzigLinMcC
# import DantzigLinB2
# import DantzigLinB10
#
# import QuadSCF
# import SCFLinCL
# import SCFLinMcC
# import SCFLinBI
# import SCFLinB2
# import SCFLinB10

import QuadMTZ
import LinMTZ.MTZLinCL
import LinMTZ.MTZLinBI
import  LinMTZ.MTZLinMcC
import LinMTZ.MTZLinB2
import LinMTZ.MTZLinB10


minsize = 5
maxsize = 31

s = False
p = 3 # we use only balanced Q for this experiment
m = 1000
skip = False

properties = ["nonneg", "negskew", "posskew", "balanced", "psd", "rankone", "ranktwo", "other", "nonnegpsd"]

prop = properties[p]

filename = "MTZLin%s.txt" % prop
file = open(filename, 'w')

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

objfilename = "MTZLin%s-obj.txt" % prop
timefilename = "MTZLin%s-time.txt" % prop
tourfilename = "MTZLin%s-tour.txt" % prop
gapfilename = "MTZLin%s-gap.txt" % prop
statusfilename = "MTZLin%s-status.txt" % prop

objfile = open(objfilename, "w")
timefile = open(timefilename, "w")
tourfile = open(tourfilename, "w")
gapfile = open(gapfilename, "w")
statusfile = open(statusfilename, "w")


for n in range(minsize, maxsize, 5):

    if n < 15:
        trials = 5
    else:
        trials = 1

    if p ==7:
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
        obj = {}
        time = {}
        tour = {}
        full = {}
        gap = {}
        status = {}



        # Triangular Q
        # name = qname + "-" + str(2)
        q2 = QMod.triangular(q, e)

        # MTZ Formulation
        obj[0, n / 5 - 1 + t], time[0, n / 5 - 1 + t], x, gap[0, n / 5 - 1 + t], status[0, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, c, q2, name)
        tour[0, n / 5 - 1 + t], full[0, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Binary replacement
        obj[1, n / 5 - 1 + t], time[1, n / 5 - 1 + t], x, gap[1, n / 5 - 1 + t], status[1, n / 5 - 1 + t] = LinMTZ.MTZLinBI.SolveTSP(n, c, q2, name)
        tour[1, n / 5 - 1 + t], full[1, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Classic
        obj[2, n / 5 - 1 + t], time[2, n / 5 - 1 + t], x, gap[2, n / 5 - 1 + t], status[2, n / 5 - 1 + t] = LinMTZ.MTZLinCL.SolveTSP(n, c, q2, name)
        tour[2, n / 5 - 1 + t], full[2, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # McCormick Envelopes
        obj[3, n / 5 - 1 + t], time[3, n / 5 - 1 + t], x, gap[3, n / 5 - 1 + t], status[3, n / 5 - 1 + t] = LinMTZ.MTZLinMcC.SolveTSP(n, c, q2, name)
        tour[3, n / 5 - 1 + t], full[3, n / 5 - 1 + t] = VerifyTour.check(x, n)
        
        # Base 2 Formulation
        obj[4, n / 5 - 1 + t], time[4, n / 5 - 1 + t], x, gap[4, n / 5 - 1 + t], status[4, n / 5 - 1 + t] = LinMTZ.MTZLinB2.SolveTSP(n, c, q2, name)
        tour[4, n / 5 - 1 + t], full[4, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Base 10 Formulation
        obj[5, n / 5 - 1 + t], time[5, n / 5 - 1 + t], x, gap[5, n / 5 - 1 + t], status[5, n / 5 - 1 + t] = LinMTZ.MTZLinB10.SolveTSP(n, c, q2, name)
        tour[5, n / 5 - 1 + t], full[5, n / 5 - 1 + t] = VerifyTour.check(x, n)


        objline = []
        timeline = []
        tours = []
        gapline = []
        statusline = []
        for i in range(6):
            objline.append(obj[i, n / 5 - 1 + t])
            a = time[i, n / 5 - 1 + t]
            time[i, n / 5 - 1 + t] = round(a, 3)
            timeline.append(time[i, n / 5 - 1 + t])
            tours.append(full[i, n / 5 - 1 + t])
            tourfile.write(str(tour[i, n / 5 - 1 + t]) + "\n")
            gapline.append(gap[i, n / 5 - 1 + t])
            statusline.append(status[i, n/5 -1 +t])
        print(objline)
        print(timeline)
        print(tours)
        print(statusline)

        file.write("%d & %g & %g & %g & %g & %g &%g \\\\  \n" % (
            n, time[0, n / 5 - 1 + t], time[1, n / 5 - 1 + t], time[2, n / 5 - 1 + t], time[3, n / 5 - 1 + t], time[4, n / 5 - 1 + t], time[5, n/5-1+t]))

        objfile.write(str(objline) + "\n")
        timefile.write(str(timeline) + "\n")
        tourfile.write(str(tours) + "\n")
        gapfile.write(str(gapline) + "\n")
        statusfile.write(str(statusline) + "\n")

        # print(str(objline))

file.write("\\end{tabular}\n } \n")
file.write("\\end{table}\n")
file.write("\\end{document}")

file.close()
objfile.close()
timefile.close()
tourfile.close()
gapfile.close()
