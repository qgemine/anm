# ANM Benchmark

In order to promote the development of computational techniques for **Active Network Management (ANM)**, we propose a test instance of a generic ANM problem. It can be used as a test bed for comparing existing computational techniques, as well as for developing new ones.

![Screenshot of a simulation run](http://www.montefiore.ulg.ac.be/~anm/anm_simulation.png)

## Problem formulation

The generic formulation of the ANM problem is motivated and described in [this paper](http://arxiv.org/pdf/1405.2806.pdf), which details also the procedure used for building the test instance.

## Installation

This code runs on Matlab&reg; and requires the Statistics and Optimization toolboxes as well as the third-party [YALMIP toolbox](http://users.isy.liu.se/johanl/yalmip/) (free and easy to install). In order to run the illustrative control policy, you will also need to intall Gurobi&reg; (free academic license).


## Manual

The implementation and interface of the test instance are both contained in the `ANM_System` class. To initialize it, just create a new object:

    inst = ANM_System();

Documentation for this class is available within Matlab&reg; by typing `doc ANM_System`. In addition, we present here some use cases that illustrate how to use the interface.

### Parameters of the test instance

TODO

### Simulation of the system

TODO

### Curtailement instructions

TODO

### Activation of flexible loads

TODO

### Using a control policy

TODO