# This program runs the experimentation for the Single Commodity Flow formulation of the quadratic traveling salesman problem. We
# provide quadratic cost files with different properties (e.g. nonnegative entries, balanced positive/negative
# entries, etc.), which are read from other files. We then modify the quadratic costs in the QMod file, and then send
#  the resulting modified cost to the appropriate TSP model. The outputs for this program are the latex file that
# produces the tables in my report, as well as simple text files of the objective value, gap value, Gurobi status,
# and the tour of the TSP. It also outputs the Gurobi log files.

import QuadSCF
import VerifyTour
import QMod
from pathlib import Path
import os

minsize = 5
maxsize =6
# We provide 100 size 5 problems (averaged), 10 size 8 problems (averaged), 5 size 10 problems, 5 size 12 problems (one solved), one of each sizes 15, 20, 25, and 30 (none solved).
# size = [5,8,10,15]
size = [12]

# Presolve value: 0 = off, -1 = default
presolve = 0

# Used for testing purposes only
skip = False
adj = False

# Set number of trials to (can do up to 100 size 5, and 5 size 10, five each size 10 and 12)
fivetrials = 0
eighttrials = 0
tentrials = 1
twelvetrials = 0
totaltrials = fivetrials + tentrials + eighttrials +twelvetrials

s = False
# p = 0  # change within 0, ..., 7 for the different properties of randomly generated quadratic cost files
# details on creation of the quadratic costs can be found in MakeTSP

for p in range(8):

    #Parameters

    # Plus/minus M value
    m = 10000

    properties = ["nonneg", "negskew", "posskew", "balanced", "psd", "rankone", "ranktwo", "nonnegpsd", "other"]

    prop = properties[p]

    filename = "SCF_%s" % prop
    file = open(filename+".txt", 'w')

    file.write("\\documentclass[11pt]{article}\n")
    file.write(
        "\\usepackage{amsmath, amssymb, amsthm, amsfonts,multirow,booktabs,siunitx, lscape, multirow, rotating, booktabs}\n")
    file.write("\\begin{document}\n")
    file.write("\\begin{table}[h] \n")
    file.write("\\centering\n ")
    file.write("\\begin{tabular}{ @{} ccccccccc @{}} \toprule\n")
    file.write(
        "Size &Original & Sym & Upper Triangular & PSD & NSD & QR & Sym QR & Upper Triangular QR \\\\ \\midrule \n")

    objfilename = "%s-obj.txt" % filename
    timefilename = "%s-time.txt" % filename
    tourfilename = "%s-tour.txt" % filename
    gapfilename = "%s-gap.txt" % filename
    statusfilename = "%s-status.txt" % filename

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

    for n in size:

        if p == 8:
            trials = 1
        elif n == 5:
            trials = fivetrials
        elif n == 8:
            trials = eighttrials
        elif n == 10:
            trials = tentrials
        elif n == 12:
            trials = twelvetrials
        else:
            trials = 1

        for t in range(trials):

            e = n * (n - 1)

            qname = "Q" + str(n) + str(prop) + "-" + str(t)

            qcname = "%s.txt" % qname

            data_q = Path("Cost/" + qcname)

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
            data_c = Path("Cost/" + cname)
            # fc = open(data_c, "r")
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


            name = qname + "-" + str(0)

            # We now go through each of the modifications, beginning with no modification, then symmeterizing, then upper triangular, etc.

            # Original Q
            obj[0,count], time[0,count], x, gap[0,count], status[0,count] = QuadSCF.SolveTSP(n, c, q, name, adj, presolve)
            tour[0,count], full[0,count] = VerifyTour.check(x, n)

            # Symmetric Q
            name = qname + "-" + str(1)
            q1 = QMod.half(q, e)
            obj[1,count], time[1,count], x, gap[1,count], status[1,count] = QuadSCF.SolveTSP(n, c, q1, name, adj, presolve)
            tour[1,count], full[1,count] = VerifyTour.check(x, n)
            # print("Here")

            # Triangular Q
            name = qname + "-" + str(2)
            q2 = QMod.triangular(q, e)
            obj[2,count], time[2,count], x, gap[2,count], status[2,count] = QuadSCF.SolveTSP(n, c, q2, name, adj, presolve)
            tour[2,count], full[2,count] = VerifyTour.check(x, n)

            # Make Q positive semi-definite
            name = qname + "-" + str(3)
            q3 = QMod.plusm(q, e, m)
            obj[3,count], time[3,count], x, gap[3,count], status[3,count]  = QuadSCF.SolveTSP(n, c, q3, name, adj, presolve)
            tour[3,count], full[3,count] = VerifyTour.check(x, n)
            obj[3,count] += -n * m

            # Make Q negative semi-definite
            name = qname + "-" + str(4)
            q4 = QMod.minusm(q, e, m)
            obj[4,count], time[4,count], x, gap[4,count], status[4,count]  = QuadSCF.SolveTSP(n, c, q4, name, adj, presolve)
            tour[4,count], full[4,count] = VerifyTour.check(x, n)
            obj[4,count] += n * m

            # Node N Removal
            name = qname + "-" + str(5)
            qr, lr = QMod.quadred2(q, e, n)

            cr = {}
            for i in range(n):
                for j in range(n):
                    cr[i, j] = c[i, j] + lr[i, j]  # Add L cost to non-quadratic cost

            obj[5,count], time[5,count], x, gap[5,count], status[5,count] = QuadSCF.SolveTSP(n, cr, qr, name, adj, presolve)
            tour[5,count], full[5,count] = VerifyTour.check(x, n)

            # Make QR into symmetric matrix
            name = qname + "-" + str(6)
            q6 = QMod.half(qr, e)
            obj[6,count], time[6,count], x, gap[6,count], status[6,count] = QuadSCF.SolveTSP(n, cr, q6, name, adj, presolve)
            tour[6,count], full[6,count] = VerifyTour.check(x, n)

            # Make QR into upper triangular
            name = qname + "-" + str(7)

            q7 = QMod.triangular(qr, e)
            obj[7,count], time[7,count], x, gap[7,count], status[7,count] = QuadSCF.SolveTSP(n, cr, q7, name, adj, presolve)
            tour[7,count], full[7,count] = VerifyTour.check(x, n)

            # Print results to files

            objline = []
            timeline = []
            tours = []
            gapline = []
            statusline = []
            for i in range(8):
                objline.append(obj[i,count])
                a = time[i,count]
                time[i,count] = round(a, 3)
                timeline.append(time[i,count])
                tours.append(full[i,count])
                tourfile.write(str(tour[i,count]) + "\n")
                gapline.append(gap[i,count])
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

    # Average all size 5 problems
    if fivetrials > 0:
        fiveavg = [0] * 8
        for i in range(8):
            fiveavg[i] = round(sum(time[i, j] for j in range(fivetrials)) / fivetrials, 3)

        print(fiveavg)

        file.write("%d & %g & %g & %g & %g & %g & %g & %g & %g \\\\  \n" % (
            5, fiveavg[0], fiveavg[1], fiveavg[2], fiveavg[3], fiveavg[4], fiveavg[5], fiveavg[6], fiveavg[7]))

    # Average all size 8 problems
    if eighttrials > 1:
        eightavg = [0] * 8
        for i in range(8):
            eightavg[i] = round(sum(time[i, j] for j in range(eighttrials)) / eighttrials, 4)

        print(eightavg)

        file.write("%d & %g & %g & %g & %g & %g & %g & %g & %g \\\\  \n" % (
            8, eightavg[0], eightavg[1], eightavg[2], eightavg[3], eightavg[4], eightavg[5], eightavg[6], eightavg[7]))

    # Report remaining trials
    for i in range(fivetrials+eighttrials, totaltrials):
        file.write("%d & %g & %g & %g & %g & %g & %g & %g & %g \\\\  \n" % (
                n, time[0, i], time[1, i], time[2,i], time[3, i], time[4, i], time[5, i],time[6, i], time[7,i]))

    file.write("\\end{tabular}\n } \n")
    file.write("\\end{table}\n")
    file.write("\\end{document}")

    file.close()
    objfile.close()
    timefile.close()
    tourfile.close()
    gapfile.close()
    statusfile.close()
