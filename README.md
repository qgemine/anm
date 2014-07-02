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

The class has several public properties that describe the characteristics of the test instance. You get obtain a list of them by using the command `properties ANM_System`. Among other information, you can find the number of buses in the distribution network, its admittance matrix, the number of distributed generators and loads, the voltage and current limits, the duration of modulation signals, ...

    >> inst = ANM_System();
    >> display(['This test system relies on a ' num2str(inst.N_buses) '-bus network,' ...
                ' to which ' num2str(inst.N_gens) ' generators and ' num2str(inst.N_loads) ...
                ' loads are connected.']);
    This test system relies on a 77-bus network, to which 59 generators and 53 loads are connected.
    
    >> display(['The average duration of modulation signals ' ...
                'is ' num2str(mean(inst.Tflex)) ' periods.']);
    The average duration of modulation signals is 14.6415 periods.


### Simulation of the system

The main method that is responsible for simulating the system is `transition()`. Every time it is called, it uses the stochastic model of the system (such as described in [this paper](http://arxiv.org/pdf/1405.2806.pdf)) to produce a transition to the next time step. It is then possible the obtain the new electrical quantities using methods `getV()`, `getI()`, `getPGens()`, `getPLoads()`, ...

    >> inst = ANM_System();
    >> for i = 1:96 % one-day long simulation
    >>     maxV = max(inst.getV()); % get the highest voltage magnitude
    >>     minV = min(inst.getV()); % get the lowest voltage magnitude
    >>     maxIuse = max(inst.getI()./inst.ratings); % get the highest link usage
    >>     PgenTot = sum(inst.getPGens()); % get the total active production
    >>     fprintf(['t=%d: maximum and minimum voltage magnitude are %1.2f and %1.2f, ' ...
    >>              'the most used link is used at %2.2f%% of its capacity and the total ' ...
    >>              'amount of distributed production is %2.2f MW.\n'], i, maxV, minV, maxIuse, PgenTot);
    >>     inst.transition(); % generate the transition of the system to the next time step
    >> end
    t=1: maximum and minimum voltage magnitude are 1.02 and 1.01, the most used link is used at 0.08% of its capacity and the total amount of distributed production is 6.31 MW.
    t=2: maximum and minimum voltage magnitude are 1.02 and 1.01, the most used link is used at 0.17% of its capacity and the total amount of distributed production is 7.78 MW.
    t=3: maximum and minimum voltage magnitude are 1.03 and 1.01, the most used link is used at 0.25% of its capacity and the total amount of distributed production is 9.29 MW.
    t=4: maximum and minimum voltage magnitude are 1.03 and 1.01, the most used link is used at 0.32% of its capacity and the total amount of distributed production is 10.49 MW.
    ...
    t=95: maximum and minimum voltage magnitude are 1.01 and 0.99, the most used link is used at 0.35% of its capacity and the total amount of distributed production is 0.14 MW.
    t=96: maximum and minimum voltage magnitude are 1.01 and 1.00, the most used link is used at 0.28% of its capacity and the total amount of distributed production is 0.40 MW.

### Curtailement instructions

TODO

### Activation of flexible loads

TODO

### Using a control policy

TODO