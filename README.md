# QTSP
Python code to support QTSP research

The Cost subfolder contains the randomly generated cost files and models for the quadratic TSP. Details on creating the different random classes of quadratic costs can be found in the "MakeTSP.py" file. These are read as the cost inputs for the different quadratic models.

To run the models, select the appropriate Experiment file (QuadMTZ, QuadSCF, QuadDantzig). You can change the property of the quadratic matrix. The model will then read a different set of costs from the Costs folder. The model will test the different modifications within the QMod file.

The outputs will be a textfile that creates the Latex files in my report, as well as simple text files of the objective values, time to run, gurobi status code, tour information, and the gurobi log files.
