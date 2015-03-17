from __future__ import division
import time
import sys
import operator
import functools
from numpy import *
from scipy import *
from sklearn.cluster import *
import copy as cpy
import coopr.environ
from coopr.pyomo import *
from coopr.opt import SolverFactory
import collections

from ANM import *

opt = SolverFactory("cplex")

def get_model(sim, T, Pgens, Ploads, dev2bus, A_hyperplane, b_hyperplane):
	model = ConcreteModel()

	# Loading modulation signals from simulator
	Pflex = {}
	for i in range(sim.N_loads):
		signal = sim.getFlexSignal(i)
		for j in range(min(sim.Tflex[i],T)):
			Pflex[(j+1,i+1)] = asscalar(signal[j])
	# Loading modulation offsets currently active
	Pflex_init = {}
	for i in range(sim.N_loads):
		flexState = sim.getFlexState()[i]
		signal = sim.getFlexSignal(i)
		for j in range(asscalar(sim.Tflex[i]-flexState),min(sim.Tflex[i],T)):
			Pflex_init[(j+1-sim.Tflex[i]+flexState,i+1)] = asscalar(signal[j])

	### Periods within the optimization horizon
	model.periods = RangeSet(1,T,name="periods")

	### Generators
	model.generators = RangeSet(1, sim.N_gens,name="generators")
	# Active power injections of generators
	model.Pgens = Param(model.periods, model.generators, name="Pgen", within=Reals, initialize=dict(zip(((i+1,j+1) for i in range(Pgens.shape[0]) for j in range(Pgens.shape[1])), (asscalar(x) for x in Pgens.reshape(Pgens.size)))))
	# Curtailment cost
	model.CurtCost = Param(model.periods, name="CurtCost", within=Reals, initialize=dict(zip(range(1,T+1),(asscalar(sim.getCurtPrice(((sim.getQuarter()+i) % 96)+1))/4.0 for i in range(T)))))

	### Loads
	model.loads = RangeSet(1, sim.N_loads, name="loads")
	# Active power injections of loads
	model.Ploads = Param(model.periods, model.loads, name="Pload", within=Reals, initialize=dict(zip(((i+1,j+1) for i in range(Ploads.shape[0]) for j in range(Ploads.shape[1])), (asscalar(x) for x in Ploads.reshape(Ploads.size)))))
	# Cost of flexibility services
	model.FlexCost = Param(model.loads, within=NonNegativeReals, initialize=dict(zip(range(1,sim.N_loads+1),(asscalar(sim.getFlexCost(i)) for i in range(sim.N_loads)))))
	# Duration of flexibility services.
	model.FlexDuration = Param(model.loads, within=PositiveIntegers, initialize=dict(zip(range(1,sim.N_loads+1),(asscalar(sim.Tflex[i]) for i in range(sim.N_loads)))))
	# Modulation signal of flexibility services.
	model.FlexSignal = Param(model.periods, model.loads, within=Reals, default=0.0, initialize=Pflex)
	# State of flexibility services at the time of the first stage.
	model.InitFlexState = Param(model.loads, within=NonNegativeIntegers, initialize=dict(zip(range(1,sim.N_loads+1),(asscalar(sim.getFlexState()[i]) for i in range(sim.N_loads)))))
	# Modulation signals decided beforehand.
	model.InitFlexSignal = Param(model.periods, model.loads, within=Reals, default=0.0, initialize=Pflex_init)

	### Buses
	model.buses = RangeSet(1, sim.N_buses-1)

	### Parameters of hyperplane A*x+b = 0 approximating operational constraints.
	model.A_hyperplane = Param(model.buses, within=Reals, initialize=dict(zip((i+1 for i in range(sim.N_buses-1)),(asscalar(x) for x in A_hyperplane))))#Param(model.buses, within=Reals)
	model.b_hyperplane = Param(initialize=b_hyperplane, within=Reals)

	### Variables
	# Curtailment instructions
	model.Pmax = Var(model.periods, model.generators, within=NonNegativeReals, bounds=(0,10))
	# Curtailed production
	model.DeltaProd = Var(model.periods, model.generators, within=NonNegativeReals)
	# Cost
	model.Cost = Var(model.periods, within=NonNegativeReals)
	# Active power injection at buses
	model.Pbuses = Var(model.periods, model.buses, within=Reals)
	# Activation of flexibility services
	model.ActFlex = Var(model.periods, model.loads, within=Binary)
	# State of flexible loads (nb of remaining active periods)
	model.FlexState = Var(model.periods, model.loads, within=NonNegativeIntegers)
	# Effective modulation of loads
	model.DeltaCons = Var(model.periods, model.loads, within=Reals)
	# Binary variables used to model max(.,.) operator
	model.MaxBin = Var(RangeSet(2,T), model.loads, within=Binary)
	model.MaxVar = Var(RangeSet(2,T), model.loads, within=NonNegativeIntegers)

	### Constraints
	# Cost assignement per time period
	def cost(model, t):
		return model.Cost[t] == sum(model.CurtCost[t]*model.DeltaProd[t,gen] for gen in model.generators)
		+ sum(model.FlexCost[load]*model.ActFlex[load,t] for load in model.loads)
	model.CostBinding = Constraint(model.periods, rule=cost)

	# Determine the amount of power curtailed
	def amount_curt(model, t, gen):
		return model.DeltaProd[t,gen] >= model.Pgens[t,gen] - model.Pmax[t,gen]
	model.AmountCurt = Constraint(model.periods, model.generators, rule=amount_curt)

	# Dertermine the power injected in buses
	def inj_in_buses(model, t, bus):
		return sum(asscalar(dev2bus[bus-1,sim.N_loads+gen-1])*(model.Pmax[t,gen]) for gen in model.generators) \
			   + sum(asscalar(dev2bus[bus-1,load-1])*(model.Ploads[t,load]+model.DeltaCons[t,load]) for load in model.loads) \
			 == model.Pbuses[t,bus]
	model.Mapping = Constraint(model.periods, model.buses, rule=inj_in_buses)

	# Force initial state of flexible loads
	def init_flex_state(model, load):
		return model.FlexState[1,load] == model.InitFlexState[load]
	model.DefInitFlexState = Constraint(model.loads, rule=init_flex_state)

	# Model MaxVar = max(0,FlexState-1)
	def max_op_rule(model, rule_num, t, load):
		if rule_num == 1:
			return model.MaxVar[t,load] <= model.FlexState[t-1,load]-1 + model.MaxBin[t,load]
		if rule_num == 2:
			return model.MaxVar[t,load] <= model.FlexDuration[load] * (1 - model.MaxBin[t,load])
		if rule_num == 3:
			return model.MaxVar[t,load] >= model.FlexState[t-1,load]-1
	model.MaxOpRule = Constraint(RangeSet(1,3), RangeSet(2,T), model.loads, rule=max_op_rule)

	# Transition rule of flexible states
	def trans_flex(model, t, load):
		return model.FlexState[t,load] == model.MaxVar[t,load] + model.ActFlex[t,load]*model.FlexDuration[load]
	model.TransFlex = Constraint(RangeSet(2,T), model.loads, rule=trans_flex)

	# Forbit double activation of a flexibility service
	def enforce_single_activation(model, t, load):
		return model.ActFlex[t,load] + (1.0/model.FlexDuration[load]) * model.FlexState[t,load] <= 1
	model.SingleAct = Constraint(model.periods, model.loads, rule=enforce_single_activation)

	# Determine effective modulation of loads
	def eff_mod_loads(model, t, load):
		return model.DeltaCons[t,load] == sum(model.ActFlex[k,load]*model.FlexSignal[1+t-k,load] for k in range(1,t+1)) + model.InitFlexSignal[t,load]
	model.EffModLoads = Constraint(model.periods, model.loads, rule=eff_mod_loads)

	# Approximation of operational constraints
	def op_consts(model, t):
		return sum(model.A_hyperplane[i]*model.Pbuses[t,i] for i in model.buses) + model.b_hyperplane <= 0
	model.OperationalConsts = Constraint(model.periods, rule=op_consts)

	### Objective
	def ObjRule(model):
		return summation(model.Cost)

	model.costFct = Objective(rule=ObjRule)

	return model

def Tree():
    return collections.defaultdict(Tree)

# Data structure used to populate list of nodes of a scenario tree
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value
# Dict-like data structre using "." to access keys.
def bunch(**kw):
    return type('bunch',(), kw)

def enum_ancestors(parent, branching, ancestors):
	if branching == []:
		return
	depth = parent[0]+1
	prior = ancestors[depth].keys()
	offset = 0
	if prior:
		offset = max(prior)
	for child in range(1,branching[0]+1):
		ancestors[depth][child+offset] = parent
		enum_ancestors((depth,child+offset), branching[1::], ancestors)
		

def policy(param, simu):
	# Number of devices
	N_dev = simu.N_gens + simu.N_loads
	# Shape and depth of the scenario tree (excluding root node)

	# Check the size of the parameter passed as argument.
	if asarray(param).shape != (simu.N_buses,):
		raise "Policy's paramater must be a vector of size " + str(simu.N_buses)

	# Split the parameter into A and b, as in "A*x + b <= 0".
	A_hyperplane = asarray(param)[0:-1]
	b_hyperplane = asscalar(asarray(param)[-1])

	# Generate 'N_samples' 'T'-long trajectories of the active power injection of the devices.
	timeseries = zeros((policy.N_samples, policy.T*N_dev))
	for k in range(0,policy.N_samples):
		inst = cpy.copy(simu)
		for t in range(0,policy.T):
			inst.transition()
			timeseries[k, (t*N_dev):((t+1)*N_dev)] = concatenate((inst.getPLoads(),inst.getPGens()))
	return timeseries
	#hc = AgglomerativeClustering(n_clusters=policy.branching[0],\
	#							 affinity='euclidean', linkage='ward')#functools.reduce(operator.mul,policy.branching, 1)\
	#hc.fit(timeseries)

	#return hc

	# Cluster the timeseries into 'N_scens' scenarios
	#km = KMeans(n_clusters=policy.N_scens)
	#km.fit(timeseries)
	#Ps = km.cluster_centers_.reshape((policy.T,N_dev))
	#Ploads = Ps[:,0:simu.N_loads:]
	#Pgens = Ps[:,simu.N_loads::]

	# Get the mathematical program instanciated with paramters' value and foreseen power injections
	# model = get_model(simu, policy.T, Pgens, Ploads, simu.dev2bus, A_hyperplane, b_hyperplane)
	# instance = model.create()

	# # Solve the mathematical program and check the termination condition of the solver.
	# results = opt.solve(instance)
	# if str(results.solver.termination_condition) != "optimal":
	# 	raise Exception("Unexpected solver status: " + str(results.solver.termination_condition))

	# # Load the results from the solver
	# instance.load(results)

	# # Gather first stage's action variables into 'action.curt' and 'action.flex'
	# actions = bunch(curt=[],flex=[])
	# actions.curt = ones(simu.N_gens)
	# for i in range(0,simu.N_gens):
	# 	actions.curt[i] = instance.Pmax[1,i+1].value
	# actions.flex = ones(simu.N_loads, dtype=integer)
	# for i in range(0,simu.N_loads):
	# 	actions.flex[i] = instance.ActFlex[1,i+1].value

	# return actions

# Helper for the user to know the expected size of the policy's parameter.
policy.param_size = Simulator.N_buses
# Policy's static parameters
policy.N_scens = 1
policy.N_samples = 100
# Shape and depth of the scenario tree (excluding root node)
policy.branching = [3,2] + [1] * 13
policy.T = len(policy.branching)
# Mapping of nodes to their ancestor
_raw_ancs = Vividict()
enum_ancestors((0,1),policy.branching,_raw_ancs)
policy.ancestors = dict([((i,j),_raw_ancs[i][j]) \
				 	   for i in _raw_ancs.keys() \
				 	   for j in _raw_ancs[i].keys()])

def evaluate(param, N_runs, L_run, discount):
	# Array of weighted instantaneous rewards
	r = zeros((N_runs,L_run))

	# Compute the N_runs * L_run rewards
	for n in range(0, N_runs):
		simu = Simulator()
		for t in range(0, L_run):
			sys.stdout.write("\r%2.2f%% completed..." % (100*(n*L_run+t)/(N_runs*L_run)))
			sys.stdout.flush()
			action = policy(param, simu)
        	simu.setCurtailment(action.curt)
        	simu.setFlexAct(action.flex)
        	simu.transition()
        	r[n,t] = (discount**t) * simu.getReward()

	sys.stdout.write("\r100% completed...\r\n")
	sys.stdout.flush()

    # return the mean and standard dev. of the policy's return
	return (r.sum(1).mean(),r.sum(1).std())

