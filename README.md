The code is associated to the paper with title: 
Collective Influence Maximization for Multiple CompetingProducts with an Awareness-to-Influence Model

## Code details
The code is written in C++ and there are a few python scripts to run experiments and plot the outputs. All the python modules can be downloaded through pip and there are no additional libraries needed for C++ code. The compiler for C++ is g++ and the version is c++11.

## Compile the code
First run:
```bash
make clean
```
Then:
```bash
make
```
## Parameters
The code has multiple input parameters. The lines below contain the flag and the potential inputs. The order of the paraeters does not matter.
#### Algoritms:
Example:
```bash
-alg EGNA
```
EG refers to GCW and NA refers to the Naive algorithm. You can run individually each algorithm as well 

#### Dataset:
Example:
```bash
-dt Gowalla
```
It will use the Gowalla dataset. There is a variety of datasets: Foursquare, Gowalla, LastFM, Dallas, test

#### Clustering:
Example:
```bash
-dt kmeans
```
For the spatial datasets we applied diferent clustering algorithms. The ones that are available are: original, kmeans, dgcd, ff

#### Box:
Example:
```bash
-box true
```
For the spatial datasets in order to randomize the locations of the events you could select to pick the locations of the events from a random box. (true, false)

#### Epsilon:
Example:
```bash
-e 0.1
```
This sets the epsilon paramter for the DSSA algorithm. It defaults to 0.1.

#### Delta:
Example:
```bash
-d 0.01
```
This sets the delta paramter for the DSSA algorithm. It defaults to 0.01.

#### Budget:
Example:
```bash
-k 10
```
This sets the budget for each competitor. It defaults to 10.

#### Competitors:
Example:
```bash
-c 16
```
This sets the number of competitors. It defaults to 16.

#### Number of MC for updating AP tables:
Example:
```bash
-mc1 100
```
This sets the number of Monte Carlo simulations to be performed at the end of every best response in order to update the Awareness Probability tables. It defaults to 100.

#### Number of MC for calculating total similarity:
Example:
```bash
-mc2 100
```
This sets the number of Monte Carlo simulations to estimate the total similarity at the end of every round. It defaults to 100.

#### Variable:
Example:
```bash
-v k
```
This sets the if the budget or number of competitors is the variable (k or c respectively). By default neither is a variable.

#### Repetitions:
Example:
```bash
-r 10
```
This sets the number of repetitions to be performed on the given experiment. It defaults to 1 where there is no varible, 5 otherwise.

#### Condition:
Example:
```bash
-cnd rounds
```
This sets the condition to stop successive rounds ("rounds", "perc" or "seeds"). "rounds" ends the loop when the given number of rounds is met. "perc" ends the loop when successive rounds do not increase the total similarity more than the given percent. "seeds" ends the loop when there is no change is the seedsets. It defaults to "rounds".

#### Rounds:
Example:
```bash
-rnd 1
```
This condition sets the number of rounds to be performed to e.g. 1. It defaults to 1.

#### Quality percent:
Example:
```bash
-p 0.03
```
This condition forces the rounds to stop if the increase in total similarity is smaller than the parameter e.g. 3%. It defaults to 0.03.

#### Propagation Model:
Example:
```bash
-m 1
```
This condition sets the propagation model either to IC ("1") or LT ("0"). It defaults to 1.

#### Experiments:
Example:
```bash
-exp 1
```
This condition sets the type of experiment to be performed. It defaults to 1.

#### Output:
Example:
```bash
-o def
```
This condition sets the path of the output file. It defaults to "def" which creates files automatically and stores them in the results folder.

## Running the code
To run the code you can set the desired parameters accordingly. Below you can find an example where, dataset is Gowalla on which kmeans was performed, the algorithms run are GCW (EG) and Naive (NA), for budget k=10 and competitors c=16 under IC.
```bash
.GCW -dt Gowalla -clust kmeans -alg EGNA -k 10 -c 16 -m 1
```
## Additional scripts
The file plot_graph.py is a script that contains a lot of helper functions to plot and run experiments easily. The file is very easy to understand and you can set the desired parameters at the top.



