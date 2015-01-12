from __future__ import division
from numpy import *
from scipy import *
from sklearn.cluster import *
import copy as cpy
import coopr.environ
from coopr.pyomo import *
from coopr.opt import SolverFactory

from ANM import *

opt = SolverFactory("glpk")

def get_model(sim, T, Pgens, Ploads, dev2bus, A_hyperplane, b_hyperplane):
	T = 1
	model = ConcreteModel()

	### Periods within the optimization horizon
	model.periods = RangeSet(1,T)

	### Generators
	model.generators = RangeSet(1, sim.N_gens)
	# Active power injections of generators
	model.Pgens = Param(model.periods, model.generators, within=Reals, initialize=dict(zip(((i+1,j+1) for i in range(Pgens.shape[0]) for j in range(Pgens.shape[1])), (asscalar(x) for x in Pgens.reshape(Pgens.size)))))
	# Curtailment cost
	model.CurtCost = Param(model.periods, within=Reals, initialize=dict(zip(range(1,T+1),(asscalar(sim.getCurtPrice(sim.getQuarter()+i+1)) for i in range(T)))))

	### Loads
	model.loads = RangeSet(1, sim.N_loads)
	# Active power injections of loads
	model.Ploads = Param(model.periods, model.loads, within=Reals, initialize=dict(zip(((i+1,j+1) for i in range(Ploads.shape[0]) for j in range(Ploads.shape[1])), (asscalar(x) for x in Ploads.reshape(Ploads.size)))))
	# State of flexibility services at the time of the first stage
	#model.PriorFlexState = Param(model.loads, within=NonNegativeIntegers)

	### Buses
	model.buses = RangeSet(1, sim.N_buses-1)
	# Active power injections at buses
	#model.Pbuses = Param(model.periods, model.buses, initialize=dict(zip(((i+1,j+1) for i in range(Pbuses.shape[0]) for j in range(Pbuses.shape[1])), Pbuses.reshape(Pbuses.size))))

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

	### Constraints
	# Cost assignement per time period
	def cost(model, t):
		return model.Cost[t] == sum(model.CurtCost[t]*model.DeltaProd[t,gen] for gen in model.generators)
	model.CostBinding = Constraint(model.periods, rule=cost)

	# Determine the amount of power curtailed
	def amount_curt(model, t, gen):
		return model.DeltaProd[t,gen] >= model.Pgens[t,gen] - model.Pmax[t,gen]
	model.AmountCurt = Constraint(model.periods, model.generators, rule=amount_curt)

	# Dertermine the power injected in buses
	def inj_in_buses(model, t, bus):
		return sum(asscalar(dev2bus[bus-1,gen-1])*(model.Pgens[t,gen]+model.DeltaProd[t,gen]) for gen in model.generators) \
			   + sum(asscalar(dev2bus[bus-1,load-1+sim.N_gens])*(model.Ploads[t,load]) for load in model.loads) \
			 == model.Pbuses[t,bus]
	model.Mapping = Constraint(model.periods, model.buses, rule=inj_in_buses)

	# Approximation of operational constraints
	def op_consts(model, t):
		return sum(model.A_hyperplane[i]*model.Pbuses[t,i] for i in model.buses) + model.b_hyperplane <= 0
	model.OperationalConsts = Constraint(model.periods, rule=op_consts)

	### Objective
	def ObjRule(model):
		return summation(model.Cost)

	model.costFct = Objective(rule=ObjRule)

	return model

def policy(A_hyperplane, b_hyperplane, simu):
	T = 1
	N_scens = 1
	N_samples = 100
	N_dev = simu.N_gens+simu.N_loads

	timeseries = zeros((N_samples, T*N_dev))

	for k in range(0,N_samples):
		inst = cpy.copy(simu)
		for t in range(0,T):
			inst.transition()
			timeseries[k, (t*N_dev):((t+1)*N_dev)] = concatenate((inst.getPLoads(), inst.getPGens()))

	km = KMeans(n_clusters=N_scens)
	km.fit(timeseries)
	Ps = km.cluster_centers_.reshape((T,N_dev))#.sum(2)

	Pgens = Ps[:,0:simu.N_gens]
	Ploads = Ps[:,simu.N_gens::]

	model = get_model(simu, T, Pgens, Ploads, simu.dev2bus, A_hyperplane, b_hyperplane)#ones(simu.N_buses-1), -15)
	
	instance = model.create()
	results = opt.solve(instance)
	instance.load(results)

	if str(results.solver.termination_condition) != "optimal":
		raise Exception("Unexpected solver status: " + str(results.solver.termination_condition))

	actions = ones(simu.N_gens)

	for i in range(0,simu.N_gens):
		actions[i] = instance.Pmax[1,i+1].value


	return actions
