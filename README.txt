This ModelDB entry contains the code for the model investigated in
 
Savelli F., "Spontaneous dynamics of hippocampal place fields in a model of
combinatorial competition among stable inputs", Journal of Neuroscience, 2024
(in press).

The code, an example of input data, and a simulation parameters/folder are
organized as follows:

The folder "ModelCode" contains all C++ code necessary to compile the
executable for running a simulation and a LICENSE file (MIT license). The
environment is GNU/Unix. Currently, the code is set up to compile and execute
in Cygwin. Porting it to another GNU environment (e.g., Linux) is
straightforward if the difference in slash-backslash conventions for building
paths is addressed (Currently, the "\" backslash MSDOS convention is used; see
below.)

The Makefile handles compilation. To compile, navigate to the folder in a
terminal and run "make". To clean compiled files before recompiling, run "make
clean".

To execute the simulation, run from the terminal: 

simulator [relative path to simulation folder] 

For example (with backslash convention):

simulator ..\SimulationPlat137\

The simulation folder is assumed to contain the file "parameters.txt"
specifying the simulation parameters. At the end of the simulation, this same
folder will be populated with several output files:

"spike_file": main output file containing which cells spiked at each timestep 

"weights_[cell_id]" (if enabled): samples of all the synaptic weights of one
place cell at regular intervals

"grid_parameters_file": geometrical parameters of every grid cell created at
the beginning of the simulation

"connectivity": which grid cells are connected to which place cells

Note that the actual trajectory data file constituting the simulation input is
stored in an external folder. The (relative) path to this folder is specified
in the parameter file in the simulation folder (see above) to avoid redundant
copies of the data if multiple simulations use the same data.

Questions can be sent to fsavelli.research@gmail.com (the same contact address
as in the article).

