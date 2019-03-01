# QTSP
Python code to support QTSP research

This repo contains python code to recreate the experiments that I completed as part of my research project work in 2018-2019. The work was overseen by Dr. Abraham Punnen, and later by Dr. Tamon Stephen, both of Simon Fraser University.

We test the effect of modifying the structure of the quadratic cost matrix on computation time for the quadratic form of the traveling salesman problem. We use randomly generated cost matrices with different properties (details below), and three formulations (Dantzig subtour elimination, Miller Tucker Zemlin - MTZ, and single commodity flow - SCF). 

We then test the effect on computation time of different linearizations on the quadratic models. We test five linearizations: binary replacement, classic method, McCormick envelopes, base-2 and base-10. Details on these linearizations can be found in my written documentation. We restrict to one quadratic cost property (balanced) and one modification (upper triangular).

To run the quadratic models, select the appropriate Experiment file (QuadMTZ, QuadSCF, QuadDantzig). You can change the property of the quadratic cost matrix by changing the property value in the model. The model will then read a different set of costs from the Costs folder. The model will test the different modifications within the QMod file.

To run the linearized models, select the appropriate Experiment file (LinMTZ, LinSCF, LinDantzig). We restrict to one property/modification, though the cost matrices are provided for the others, and the experiment file can access the QMod file if so desired. The code will test the different linearizations.

The outputs will be a textfile that creates the Latex files in my report, as well as simple text files of the objective values, time to run, gurobi status code, tour information, and the gurobi log files.


The subfolders are organized as follows:
The Cost subfolder contains the randomly generated cost files and models for the quadratic TSP. Details on creating the different random classes of quadratic costs can be found in the "MakeTSP.py" file. These are read as the cost inputs for the different quadratic models.

The LinDantzig folder contains the linearized quadratic models for the Dantzig subtour elimination formulation. 

The LinSCF folder contains the linearized quadratic models for the SCF formulation.

The LinMTZ folder contains the linearized quadratic models for the MTZ formulation. 
