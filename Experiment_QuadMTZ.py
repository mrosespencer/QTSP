# This program runs the experimentation for the MTZ formulation of the quadratic traveling salesman problem. We
# provide quadratic cost files with different properties (e.g. nonnegative entries, balanced positive/negative
# entries, etc.), which are read from other files. the outputs for this program are the latex file that produces the
# tables in my report, as well as simple text files of the objective value, gap value, Gurobi status, and the tour of
#  the TSP. It also outputs the Gurobi log files.

import QuadMTZ
import VerifyTour
import QMod
from pathlib import Path
import os

minsize = 5
maxsize = 31

s = False
p = 0  # change within 0, ..., 7 for the different properties of randomly generated quadratic cost files
# details on creation of the quadratic costs can be found in MakeTSP

m = 1000
skip = False

properties = ["nonneg", "negskew", "posskew", "balanced", "psd", "rankone", "ranktwo", "other", "nonnegpsd"]

prop = properties[p]

filename = "MTZ_%s.txt" % prop
file = open(filename, 'w')

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

for n in range(minsize, maxsize, 5):

    if n < 15:
        trials = 5
    else:
        trials = 1

    if p == 7:
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

        obj = {}
        time = {}
        tour = {}
        full = {}
        gap = {}
        status = {}

        name = qname + "-" + str(0)

        # Original Q
        obj[0, n / 5 - 1 + t], time[0, n / 5 - 1 + t], x, gap[0, n / 5 - 1 + t], status[0, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, c, q, name)
        tour[0, n / 5 - 1 + t], full[0, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Symmetric Q
        name = qname + "-" + str(1)
        q1 = QMod.half(q, e)
        obj[1, n / 5 - 1 + t], time[1, n / 5 - 1 + t], x, gap[1, n / 5 - 1 + t], status[1, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, c, q1, name)
        tour[1, n / 5 - 1 + t], full[1, n / 5 - 1 + t] = VerifyTour.check(x, n)
        # print("Here")

        # Triangular Q
        name = qname + "-" + str(2)
        q2 = QMod.triangular(q, e)
        obj[2, n / 5 - 1 + t], time[2, n / 5 - 1 + t], x, gap[2, n / 5 - 1 + t], status[2, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, c, q2, name)
        tour[2, n / 5 - 1 + t], full[2, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Make Q positive semi-definite
        name = qname + "-" + str(3)
        q3 = QMod.plusm(q, e, m)
        obj[3, n / 5 - 1 + t], time[3, n / 5 - 1 + t], x, gap[3, n / 5 - 1 + t], status[3, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, c, q3, name)
        tour[3, n / 5 - 1 + t], full[3, n / 5 - 1 + t] = VerifyTour.check(x, n)
        obj[3, n / 5 - 1 + t] = obj[3, n / 5 - 1 + t] - n * m

        # Make Q negative semi-definite
        name = qname + "-" + str(4)
        q4 = QMod.minusm(q, e, m)
        obj[4, n / 5 - 1 + t], time[4, n / 5 - 1 + t], x, gap[4, n / 5 - 1 + t], status[4, n / 5 - 1 + t]  = QuadMTZ.SolveTSP(n, c, q4, name)
        tour[4, n / 5 - 1 + t], full[4, n / 5 - 1 + t] = VerifyTour.check(x, n)
        obj[4, n / 5 - 1 + t] = obj[4, n / 5 - 1 + t] + n * m

        # Node N Removal
        name = qname + "-" + str(5)
        qr, lr = QMod.quadred2(q, e, n)

        cr = {}
        for i in range(n):
            for j in range(n):
                cr[i, j] = c[i, j] + lr[i, j]  # Add L cost to non-quadratic cost

        obj[5, n / 5 - 1 + t], time[5, n / 5 - 1 + t], x, gap[5, n / 5 - 1 + t], status[5, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, cr, qr, name)
        tour[5, n / 5 - 1 + t], full[5, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Make QR into symmetric matrix
        name = qname + "-" + str(6)
        q6 = QMod.half(qr, e)
        obj[6, n / 5 - 1 + t], time[6, n / 5 - 1 + t], x, gap[6, n / 5 - 1 + t], status[6, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, cr, q6, name)
        tour[6, n / 5 - 1 + t], full[6, n / 5 - 1 + t] = VerifyTour.check(x, n)

        # Make QR into upper triangular
        name = qname + "-" + str(7)

        q7 = QMod.triangular(qr, e)
        obj[7, n / 5 - 1 + t], time[7, n / 5 - 1 + t], x, gap[7, n / 5 - 1 + t], status[7, n / 5 - 1 + t] = QuadMTZ.SolveTSP(n, cr, q7, name)
        tour[7, n / 5 - 1 + t], full[7, n / 5 - 1 + t] = VerifyTour.check(x, n)

        objline = []
        timeline = []
        tours = []
        gapline = []
        statusline = []
        for i in range(8):
            objline.append(obj[i, n / 5 - 1 + t])
            a = time[i, n / 5 - 1 + t]
            time[i, n / 5 - 1 + t] = round(a, 3)
            timeline.append(time[i, n / 5 - 1 + t])
            tours.append(full[i, n / 5 - 1 + t])
            tourfile.write(str(tour[i, n / 5 - 1 + t]) + "\n")
            gapline.append(gap[i, n / 5 - 1 + t])
            statusline.append(status[i, n / 5 - 1 + t])
        print(objline)
        print(timeline)
        print(tours)
        print(statusline)

        file.write("%d & %g & %g & %g & %g & %g & %g & %g & %g \\\\  \n" % (
            n, obj[0, n / 5 - 1 + t], time[1, n / 5 - 1 + t], time[2, n / 5 - 1 + t],
            time[3, n / 5 - 1 + t],
            time[4, n / 5 - 1 + t], time[5, n / 5 - 1 + t],
            time[6, n / 5 - 1 + t], time[7, n / 5 - 1 + t]))

        objfile.write(str(objline) + "\n")
        timefile.write(str(timeline) + "\n")
        tourfile.write(str(tours) + "\n")
        gapfile.write(str(gapline) + "\n")
        statusfile.write(str(statusline) + "\n")

file.write("\\end{tabular}\n } \n")
file.write("\\end{sidewaystable}\n")
file.write("\\end{document}")

file.close()
objfile.close()
timefile.close()
tourfile.close()
gapfile.close()
statusfile.close()