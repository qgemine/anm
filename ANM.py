#!/usr/bin/env python

from numpy import *
from numpy.random import RandomState
from scipy import *
from scipy import optimize
from scipy import sparse
from models import WindSampler, SunSampler

class Simulator:
    # ANM.Simulator Simulate the system dynamics of the Active Network
    # Management benchmark available at http://www.montefiore.ulg.ac.be/~anm/ .
    # It takes into account the control actions provided by the user and
    # computes the reward associated to every transition.

    def __init__(self, case, rng=None, wind=WindSampler(), sun=SunSampler()):
        # __init__ Initialize the data members of the class and 
        # generate the initial state of the system.

        # Load the test case
        self.buses = case["bus"]
        self.branches = case["branch"]
        self.devices = case["gen"]
        self.services = case["flex"]
        self.baseMVA = case["baseMVA"]
        # Dictionary of parameters that require a call to a user-defined function to update their value at each transition.
        self.callables = dict()

        # Number of buses of the network.
        self.N_buses = self.buses.shape[0]
        # Number of links in the network.
        self.N_branches = self.branches.shape[0]
        # 0-based indices of the loads and generators within the list of devices.
        loads = (self.devices[:,1]<0.0).nonzero()[0].tolist()
        gens = (self.devices[:,1]>0.0).nonzero()[0].tolist()
        # Number of loads.
        self.N_loads = len(loads)
        # Number of generators.
        self.N_gens = len(gens)
        # Q-P ratio of loads and generators.
        self.qp = self.devices[loads+gens,2]/self.devices[loads+gens,1]
        # Bounds of voltage magnitude at buses.
        self.Vmin = asarray(self.buses[1::,12])
        self.Vmax = asarray(self.buses[1::,11])

        # The IDs of buses are edited to be an increasing list of integers (1,2,3,...,N_buses)
        bus_idx = dict()
        for i in range(self.N_buses):
            bus_idx[int(self.buses[i,0])] = i+1
            self.buses[i,0] = i+1

        # Update the bus ids in the array of the devices
        for i in range(self.devices.shape[0]):
            self.devices[i,0] = bus_idx[self.devices[i,0]]

        # Check the correctness of the slack bus specifications and get its voltage magnitude.
        if self.buses[0,1] == 3:
            PV_dev = [dev for dev in range(self.devices.shape[0]) if dev not in (loads+gens)]
            if len(PV_dev) > 1:
                raise NotImplementedError("A single PV device is supported at the moment (at the slack).")
            if len(PV_dev) < 1:
                raise ValueError("The test system does not have a PV device at the slack bus.")
            PV_id = PV_dev[0]
            if PV_id != 0:
                raise ValueError("The PV device at the slack bus must be specified as the first device in the input file.")
            if int(self.devices[PV_id,0]) != 1:
                raise ValueError("The test system must have a single PV device, connected at the slack bus.")
            self.V_slack = self.devices[PV_id,5]
            if (self.V_slack < 0.5) or (self.V_slack > 1.5):
                print "Warning: voltage magnitude ("+str(self.V_slack)+") at the slack bus does not seem coherent."
        else:
            raise ValueError("The slack bus of the test case must be specified as the first bus in the input file.")

        # Build the mapping matrix M such that M_ij = 1 if device j is at bus i, and 0 otherwise.
        # Bus indices are shifted by 2 to have a 0-based indexing and because the slack bus is not considered.
        # Device indices are shifted by 1 because the PV-like device at the slack bus is not considered.
        d2b_cols = []
        d2b_rows = []
        for load in loads:
            d2b_cols += [load-1]
            d2b_rows += [int(self.devices[load,0])-2]
        for gen in gens:
            d2b_cols += [gen-1]
            d2b_rows += [int(self.devices[gen,0])-2]
        d2b_data = array([1]*(self.N_loads+self.N_gens))
        self.dev2bus = sparse.coo_matrix((d2b_data, (d2b_rows, d2b_cols)),
                                              shape=(self.N_buses-1, self.N_loads+self.N_gens)).toarray()

        # Update bus indices in branch specifications to match the new indexing (1,...,N_buses)
        for i in range(self.branches.shape[0]):
            self.branches[i,0] = bus_idx[self.branches[i,0]]
            self.branches[i,1] = bus_idx[self.branches[i,1]]
        # Build the nodal admittance matrix from the branch specifications.
        self.Y_bus = zeros((self.N_buses,self.N_buses), dtype=complex)
        self.g_ij = empty((self.N_branches,))
        self.b_ij = empty((self.N_branches,))
        self.V_to_I  = []
        self.ratings = empty((self.N_branches,))
        self.links = []
        for i in range(self.N_buses):
            # Shunt admittance at buses
            self.Y_bus[i,i] = (self.buses[i,4]+1.0*self.buses[i,5])/self.baseMVA
        for i in range(self.N_branches):
            k = int(self.branches[i,0]-1) # from bus, 0-based id
            m = int(self.branches[i,1]-1) # to bus, 0-based id
            # Serie admittance
            r = self.branches[i,2]
            x = self.branches[i,3]
            y_line = 1.0/(r+1.0j*x)
            # Charging admittance
            b = self.branches[i,4]
            y_shunt = 0.5j*b
            # Tap ratios
            tap = self.branches[i,8] if self.branches[i,8] > 0.0 else 1.0
            shift = self.branches[i,9]*pi/180.0
            tap = tap * exp(1.0j*shift)

            self.links += [[k,m]]
            self.g_ij[i] = y_line.real
            self.b_ij[i] = y_line.imag
            self.ratings[i] = self.branches[i,5]/self.baseMVA
            self.Y_bus[k,k] += (y_line + y_shunt) / (tap*conj(tap))
            self.Y_bus[m,m] += (y_line + y_shunt)
            self.Y_bus[k,m] = -y_line / conj(tap)
            self.Y_bus[m,k] = -y_line / tap
            self.V_to_I += [lambda V_buses,k=k,m=m,tap=tap,y_line=y_line: ((tap*conj(tap))*V_buses[k]-conj(tap)*V_buses[m])*y_line]

        self.links = array(self.links)

        # Initial value of reals and imaginary parts of voltages for power flow computations.
        self.V_init = array([self.V_slack] + [1]*(self.N_buses-1) + [0]*self.N_buses)

        # Get functions to compute power levels of generators and loads.
        self.Ploads_fcts = [None]*self.N_loads
        self.Pgens_fcts = [None]*self.N_gens
        for id, fct in case["power"]:
            if id in loads:
                self.Ploads_fcts[loads.index(id)] = fct
            elif id in gens:
                self.Pgens_fcts[gens.index(id)] = fct
            else:
                raise ValueError("Specifying a function for an unknown device.")

        if any([x is None for x in self.Ploads_fcts]):
            raise ValueError("A power level function is missing for at least one load.")

        if any([x is None for x in self.Pgens_fcts]):
            raise ValueError("A power level function is missing for at least one generator.")

        # Load services of flexible loads from the input file.
        self.flexSignals = []
        self.flexT = []
        self.flexPrice = []
        self.flexIdInLoads = []
        flexLoads = [flex[1::] for flex in self.services if int(flex[0])==1]
        for loadID, price, signal in flexLoads:
            try:
                id = loads.index(int(loadID))
            except ValueError:
                raise ValueError("A flexible service is attached to a non-existing load.")
            self.flexIdInLoads += [id]
            self.flexSignals += [signal]
            self.flexT += [len(signal)]
            if callable(price):
                self.callables[("flexPrice",len(self.flexPrice))] = price
                self.flexPrice += [None]
            else:
                self.flexPrice += [price]
        self.N_flex = len(self.flexIdInLoads)
        self.flexT = asarray(self.flexT)

        # Load curtailment services from input.
        self.curtPrice = []
        self.curtIdInGens = []
        curtGens = [curt[1::] for curt in self.services if int(curt[0])==2]
        for genID, price in curtGens:
            try:
                id = gens.index(int(genID))
            except ValueError:
                raise ValueError("A curtailment service is attached to a non-existing generator.")
            self.curtIdInGens += [id]
            if callable(price):
                self.callables[("curtPrice",len(self.curtPrice))] = price
                self.curtPrice += [None]
            else:
                self.curtPrice += [price]
        self.N_curt = len(self.curtIdInGens)

        # Initialize random generator.
        self.rng = RandomState() if rng is None else rng

        # Get weather samplers.
        self.sun = sun
        self.wind = wind

        # Stochastic initialization of weather's state.
        self.ws = self.wind(self.rng)
        self.ir = self.sun(self.rng)

        # Initialize consumption of loads
        self.Ploads = array([f(self.rng) for f in self.Ploads_fcts])

        # Deduce initial distributed generation from weather data
        self.Pgens = array([f(self.ir, self.ws) for f in self.Pgens_fcts])

        # Electrical quantities are not computed initially.
        self.computed = False

        # Initialize state and action vectors.
        self.flex_state = zeros(self.N_flex, dtype=integer)
        self.dP_flex = zeros(self.N_flex)
        self.last_flex_act = zeros(self.N_flex, dtype=integer)
        self.flex_act = zeros(self.N_flex, dtype=integer)
        self.prod_limit = 1e2*ones(self.N_curt)
        self.curt_insts = 1e2*ones(self.N_curt)
        self.q = 1

    def getCurtPrice(self, gen, timesteps=0):
        # GETCURTPRICE Get the cost of curtailing generation per MWh 
        # for the specified generator and for time period 't+TIMESTEPS',
        # 't' being the current period.
        #   C = GETCURTPRICE(OBJ, GEN, TIMESTEPS) or
        #   C = OBJ.GETCURTPRICE(GEN, TIMESTEPS) returns in C the cost of
        #   curtailment per MWh for period t+TIMESTEPS that belongs to
        #   {1,...,96} and for GEN that belongs to {0,...,N_gens-1}.
        try:
            if self.curtPrice[self.curtIdInGens.index(gen)] is None:
                return self.callables[("curtPrice",self.curtIdInGens.index(gen))](self.getQuarter(), self.getQuarter(timesteps))
            else:
                return self.curtPrice[self.curtIdInGens.index(gen)]
        except ValueError:
            raise ValueError("No curtailment service for generator "+str(gen)+".")

    def getProdLimit(self, gen):
        # GETCURTSTATE Get the last curtailment instructions, i.e. the
        # power production limits that are currently active.
        try:
            return self.prod_limit[self.curtIdInGens.index(gen)]
        except ValueError:
            raise ValueError("No curtailment service for generator "+str(gen)+".")

    def getFlexCost(self, load, timesteps=0):
        # GETFLEXCOST Get the activation cost of the flexibility service
        # associated to the specified load.
        # LOAD must belong to {0,...,N_loads-1} and have a flexibility
        # service.
        try:
            if self.flexPrice[self.flexIdInLoads.index(load)] is None:
                return self.callables[("flexPrice",self.flexIdInLoads.index(load))](self.getQuarter(), self.getQuarter(timesteps))
            else:
                return self.flexPrice[self.flexIdInLoads.index(load)]
        except ValueError:
            raise ValueError("No flexibility service for load "+str(load)+".")

    def getFlexOffset(self, load):
        # GETFLEXOFFSETS Get the amount of power demand shift currently provided by the
        # flexible service of the load.
        try:
            return self.dP_flex[self.flexIdInLoads.index(load)]
        except ValueError:
            raise ValueError("No flexibility service for load "+str(load)+".")

    def getFlexSignal(self, load):
        # GETFLXSIGNAL Get a list of the offset values for the whole
        # modulation signal of the specified load.
        # LOAD must belong to {0,...,N_loads-1} and have a flexibility
        # service.
        try:
            return self.flexSignals[self.flexIdInLoads.index(load)]
        except ValueError:
            raise ValueError("No flexibility service for load "+str(load)+".")

    def getFlexState(self, load):
        # GETFLEXSTATE Get the state (number of active periods left) of
        # the flexible load. Zero means that the flexible service is inactive.
        try:
            return self.flex_state[self.flexIdInLoads.index(load)]
        except ValueError:
            raise ValueError("No flexibility service for load "+str(load)+".")
    
    def getI(self, branch):
        # GETI Get the magnitude of the current magnitude in the branch, in per unit.
        
        if not self.computed:
            self.comp_elec_state()

        # Get current magnitude [p.u.] in links.
        try:
            return self.I[branch]
        except IndexError:
            raise IndexError("Wrong branch id '"+str(branch)+"'.")
    
    def getPCurtGen(self, gen):
        # GETPCURTGENS Get the active power injection of the generator, in
        # MW, accounting for prospective curtailment instructions.
        try:
            return minimum(self.Pgens[gen], self.prod_limit[self.curtIdInGens.index(gen)])
        except ValueError:
             return self.Pgens[gen]
        except IndexError:
            raise IndexError("Wrong generator id '"+str(gen)+"'.")
    
    def getPGen(self, gen):
        # GETPGEN Get the potential (without curtailment) active power
        # injection of the generator, in MW.
        try:
            return self.Pgens[gen]
        except IndexError:
            raise IndexError("Wrong generator id '"+str(gen)+"'.")
            
    def getPLoad(self, load):
        # GETPLOAD Get the baseline (without modulation) active power
        # injection of loads, in MW.
        try:
            return -self.Ploads[load]
        except IndexError:
            raise IndexError("Wrong load id '"+str(load)+"'.")
    
    def getPModLoad(self, load):
        # GETPMODLOAD Get the active power injection of loads, in MW,
        # accounting for prospective modulations.
        try:
            return -self.Ploads[load] + self.dP_flex[self.flexIdInLoads.index(load)]
        except ValueError:
            return -self.Ploads[load]
        except IndexError:
            raise IndexError("Wrong load id '"+str(load)+"'.")

    def getQuarter(self, timesteps=0):
        # GETQUARTER Get the quarter of an hour in day of the current
        # time period.

        return ((self.q+timesteps-1)%96)+1

    def getCost(self):
        if not self.computed:
            self.comp_elec_state()

        # Compute last transition's operational costs
        curt_cost = sum([self.getCurtPrice(gen)*(self.getPGen(gen)-self.getPCurtGen(gen))/4.0 for gen in self.curtIdInGens])
        flex_cost = sum([self.getFlexCost(load)*activated for load, activated in zip(self.flexIdInLoads, self.last_flex_act.tolist())])
        return curt_cost + flex_cost

    def isSafe(self):
        if not self.computed:
            self.comp_elec_state()

        return self.isCurrentSafe() and self.isVoltageSafe()

    def isVoltageSafe(self):
        if not self.computed:
            self.comp_elec_state()

        return not any(greater(self.Vmagn, self.Vmax).tolist() + less(self.Vmagn, self.Vmin).tolist())

    def isCurrentSafe(self):
        if not self.computed:
            self.comp_elec_state()

        return not any(greater(self.I, self.ratings).tolist())

    def getViolationMagn(self):
        Vu = self.Vmagn-self.Vmax
        Vl = self.Vmin-self.Vmagn
        Iu = self.I-self.ratings
        Vu = Vu[greater(Vu,0.0).nonzero()[0]]
        Vl = Vl[greater(Vl,0.0).nonzero()[0]]
        Iu = Iu[greater(Iu,0.0).nonzero()[0]]
        return Vu.sum()+Vl.sum()+Iu.sum()

    def getReward(self):
        # GETREWARD Return the instantaneous reward associeted to the
        # last transition of the system.
        
        # Compute last transition's reward
        if self.isSafe():
            return -self.getCost()
        else:
            return -self.getCost() - 1.0e3*min(exp(self.getViolationMagn())-1.0,1.0e3)

    def getSolarIr(self):
        # GETSOLARIR Get the current solar irradiance, in W per m^2.

        return self.ir
    
    def getV(self, bus):
       # GETV Get the magnitude of the voltages at the bus, in per
       # unit.
        
        if not self.computed:
            self.comp_elec_state()

        if bus==0:
            return self.V_slack

        try:
            return self.Vmagn[bus-1]
        except IndexError:
            raise IndexError("Wrong bus id '"+str(bus)+"'.")

    def getWindSpeed(self):
        # GETWINDSPEED Get the current wind speed, in meters per second.
        
        return self.ws

    def setPmax(self, gen, p_max):
        try:
            self.curt_insts[self.curtIdInGens.index(gen)] = p_max
        except ValueError:
            raise ValueError("No curtailment service for generator "+str(gen)+".")

    def activateFlex(self, load):
        try:
            self.flex_act[self.flexIdInLoads.index(load)] = 1
        except ValueError:
            raise ValueError("No flexibility service for load "+str(load)+".")

    def transition(self):
        # TRANSITION Simulate a transistion of the system, accounting
        # for the control actions provided through SETCURTAILMENT and 
        # SETFLEXACT functions. The reward associated to this
        # transition can be obtained by calling GETREWARD function.
            
        # Electrical state is now out of date.
        self.computed = False

        # Update state of flexible loads.
        self.flex_state = maximum(self.flex_state-1, 0) + self.flex_act*self.flexT

        # Update curtailment state.
        self.prod_limit = self.curt_insts
        
        # Compute value of load modulation signals.
        self.dP_flex = greater(self.flex_state,0)*array([signal[T-s-1] for signal, T, s in zip(self.flexSignals,self.flexT,self.flex_state)])

        # Increment the quarter of hour of the day
        self.q = (self.q % 96) + 1

        # Stochastic transition of the consumption of loads.
        self.Ploads = array([f(self.rng) for f in self.Ploads_fcts])

        # Stochastic transition of the weather (irradiance and wind speed)
        self.ws = self.wind(self.rng)
        self.ir = self.sun(self.rng)

        # Distributed generation from weather data
        self.Pgens = array([f(self.ir, self.ws) for f in self.Pgens_fcts])

        # Reset instructions for next period.
        self.last_flex_act = self.flex_act
        self.flex_act = zeros(self.N_flex, dtype=integer)
        self.curt_insts = 1e2*ones(self.N_curt)

    def pf_eqs(self,v,y,p,q):
        V = (v[0:self.N_buses]+1j*v[self.N_buses:])
        Sinj = y.conj().dot(V.conj())*V
        return hstack(( Sinj[1:].real-p,
                        Sinj[1:].imag-q,
                        [V[0]-self.V_slack],
                        [V[0].imag])).real

    def comp_elec_state(self):
        # Aggregate power injections at buses
        P_devices = array([self.getPModLoad(load) for load in range(self.N_loads)]+[self.getPCurtGen(gen) for gen in range(self.N_gens)])
        Pbus = self.dev2bus.dot(P_devices)/self.baseMVA
        Qbus = self.dev2bus.dot(P_devices*self.qp)/self.baseMVA
        
        # Solve power flow equations
        x = optimize.root(self.pf_eqs, self.V_init, args=(self.Y_bus,Pbus,Qbus), method='lm', options={'xtol' : 1.0e-4}).x
        V_sol = x[0:self.N_buses]+1j*x[self.N_buses:]
        
        # Get the solution
        self.Vmagn = abs(V_sol[1:])
        self.I = array([abs(self.V_to_I[br](V_sol)) for br in range(self.N_branches)])

        # Indicate that the electrical state is now computed.
        self.computed = True

