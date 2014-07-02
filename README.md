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
    >>     fprintf(['t=%d: maximum and minimum voltage magnitudes are %1.2f and %1.2f, ' ...
    >>              'the most used link is used at %2.2f%% of its capacity and the total ' ...
    >>              'amount of distributed production is %2.2f MW.\n'], i, maxV, minV, 100*maxIuse, PgenTot);
    >>     inst.transition(); % generate the transition of the system to the next time step
    >> end
    t=1: maximum and minimum voltage magnitudes are 1.01 and 1.00, the most used link is used at 15.73% of its capacity and the total amount of distributed production is 0.73 MW.
    t=2: maximum and minimum voltage magnitudes are 1.01 and 1.01, the most used link is used at 10.80% of its capacity and the total amount of distributed production is 1.45 MW.
    ...
    t=39: maximum and minimum voltage magnitudes are 1.03 and 1.01, the most used link is used at 35.59% of its capacity and the total amount of distributed production is 17.93 MW.
    t=40: maximum and minimum voltage magnitudes are 1.03 and 1.01, the most used link is used at 38.64% of its capacity and the total amount of distributed production is 18.33 MW.
    ...
    t=95: maximum and minimum voltage magnitudes are 1.01 and 0.99, the most used link is used at 29.48% of its capacity and the total amount of distributed production is 1.17 MW.
    t=96: maximum and minimum voltage magnitudes are 1.01 and 1.00, the most used link is used at 26.50% of its capacity and the total amount of distributed production is 0.56 MW.

### Curtailment instructions

Before calling the method `transition()`, it is possible to provide curtailment instructions by using the method `setCurtailment()`. It expects as an argument a column vector with as many components as the number of generators. The i-th component corresponds to the upper limit in MW for i-th generator for the next time step.

    >> inst = ANM_System();
    >> for i = 1:96
    >>     % Percentage of curtailed production.
    >>     curt = (sum(inst.getPGens())-sum(inst.getPCurtGens()))/sum(inst.getPGens());
    >>     if sum(inst.getPGens()) == 0
    >>         curt = 0;
    >>     end
    >>     fprintf('%2.2f%% of production has been curtailed at time step %d.\n', 100*curt, i);
    >>     % Random curtailment instructions.
    >>     inst.setCurtailment(10*rand(inst.N_gens,1));
    >>     % Generate next time step.
    >>     inst.transition();
    >> end
    0.00% of production has been curtailed at time step 1.
    0.00% of production has been curtailed at time step 2.
    ...
    45.77% of production has been curtailed at time step 45.
    12.79% of production has been curtailed at time step 46.
    ...
    0.00% of production has been curtailed at time step 95.
    10.80% of production has been curtailed at time step 96.

### Activation of flexible loads

TODO

### Using a control policy

TODO