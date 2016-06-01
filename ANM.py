#!/usr/bin/env python

from numpy import *
from numpy.random import RandomState
from scipy import *
from scipy import optimize
from scipy import sparse
from models import WindSampler, SunSampler
import warnings

class Simulator:
    # ANM.Simulator Simulate the system dynamics of the Active Network
    # Management benchmark available at http://www.montefiore.ulg.ac.be/~anm/ .
    # It takes into account the control actions provided by the user and
    # computes the reward associated to every transition.

    def __init__(self, case, rng=None, wind=WindSampler, sun=SunSampler):
        # __init__ Initialize the data members of the class and 
        # generate the initial state of the system.

        # Load the test case
        self.buses = case["bus"]
        self.branches = case["branch"][case["branch"][:,10]==1.,:]
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
        # P-Q rules of devices
        self.pq_rules = list()
        for load in (self.devices[dev_idx,:] for dev_idx in loads):
            rule = PowerCapabilities("load (fixed pf)")
            rule.qp_ratio = load[Simulator.dcol["Qg"]]/load[Simulator.dcol["Pg"]]
            self.pq_rules.append(rule)
        for gen in (self.devices[dev_idx,:] for dev_idx in gens):
            rule = PowerCapabilities("generator")
            if gen[Simulator.dcol["Pg"]] != 0:
                rule.qp_ratio = gen[Simulator.dcol["Qg"]]/gen[Simulator.dcol["Pg"]]
            if gen[Simulator.dcol["Pmin"]] != gen[Simulator.dcol["Pmax"]]:
                rule.pmax = gen[Simulator.dcol["Pmax"]]
                rule.pmin = gen[Simulator.dcol["Pmin"]]
            if gen[Simulator.dcol["Qmin"]] != gen[Simulator.dcol["Qmax"]]:
                rule.qmax = gen[Simulator.dcol["Qmax"]]
                rule.qmin = gen[Simulator.dcol["Qmin"]]
            if gen[Simulator.dcol["Pc1"]] != gen[Simulator.dcol["Pc2"]]:
                rule.lead_limit = {k: gen[Simulator.dcol[k]] for k in ('Pc1','Pc2','Qc1max','Qc2max')}
                rule.lag_limit = {k: gen[Simulator.dcol[k]] for k in ('Pc1','Pc2','Qc1min','Qc2min')}
            self.pq_rules.append(rule)

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
        self.V_to_I_params  = []
        self.ratings = empty((self.N_branches,))
        self.links = []
        for i in range(self.N_buses):
            # Shunt admittance at buses
            self.Y_bus[i,i] = (self.buses[i,4]+1.0j*self.buses[i,5])
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
            self.V_to_I_params += [(k,m,tap,y_line)]

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
        self.sun = sun()
        self.wind = wind()

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
        self.Q_setpoints = empty(self.N_curt)
        self.Q_setpoints[:] = NaN
        self.Q_insts = empty(self.N_curt)
        self.Q_insts[:] = NaN
        self.q = 1

    def getActiveLosses(self):
        #return sum([self.g_ij[branch_idx]*(self.getI(branch_idx)**2) for branch_idx in range(self.N_branches)])
        if not self.computed:
            self.comp_elec_state()
        Plosses = self.S_slack.real+array([self.getPModLoad(load) for load in range(self.N_loads)]+[self.getPCurtGen(gen) for gen in range(self.N_gens)]).sum()
        return Plosses

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

    def getLossPrice(self,timesteps=0):
        return mean([self.getCurtPrice(gen, timesteps) for gen in self.curtIdInGens])

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

    def getQCurtGen(self, gen):
        # GETQCURTGENS Get the reactive power injection of the generator, in
        # MW.
        try:
            return self.Q_setpoints[self.curtIdInGens.index(gen)]
        except ValueError:
            return self.pq_rules[self.N_loads+gen].qp_ratio*self.getPCurtGen(gen)
    
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

    def getQModLoad(self, load):
        # GETQMODLOAD Get the reactive power injection of loads, in MW,
        # accounting for prospective modulations.
        return self.pq_rules[load].qp_ratio*self.getPModLoad(load)

    def getQuarter(self, timesteps=0):
        # GETQUARTER Get the quarter of an hour in day of the current
        # time period.

        return ((self.q+timesteps-1)%96)+1

    def getCost(self):
        if not self.computed:
            self.comp_elec_state()

        # Compute last transition's operational costs
        curt_cost = sum([self.getCurtPrice(gen) * (self.getPGen(gen)-self.getPCurtGen(gen)) / 4.0 for gen in self.curtIdInGens])
        flex_cost = sum([self.getFlexCost(load)*activated for load, activated in zip(self.flexIdInLoads, self.last_flex_act.tolist())])
        loss_cost = self.getLossPrice()*self.getActiveLosses()/4.0
        return curt_cost + flex_cost + loss_cost

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

    def setQgen(self, gen, q_setpoint):
        try:
            self.Q_insts[self.curtIdInGens.index(gen)] = q_setpoint
        except IndexError:
            raise ValueError("No reactive power service for generator "+str(gen)+".")

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

        # Stochastic transition of the consumption of loads.
        self.Ploads = array([f(self.rng) for f in self.Ploads_fcts])

        # Stochastic transition of the weather (irradiance and wind speed)
        self.ws = self.wind(self.rng)
        self.ir = self.sun(self.rng)

        # Distributed active generation from weather data
        self.Pgens = array([f(self.ir, self.ws) for f in self.Pgens_fcts])

        # Increment the quarter of hour of the day
        self.q = (self.q % 96) + 1

        # Update state of flexible loads.
        self.flex_state = maximum(self.flex_state-1, 0) + self.flex_act*self.flexT

        # Update curtailment state.
        self.prod_limit = self.curt_insts

        # Update reactive setpoints.
        self.Q_setpoints = copy(self.Q_insts)
        for gen in self.curtIdInGens:
            q_setpoint = self.Q_setpoints[self.curtIdInGens.index(gen)]
            if not isnan(q_setpoint):
                self.check_power_setpoint(gen, self.getPCurtGen(gen), q_setpoint)
        
        # Compute value of load modulation signals.
        self.dP_flex = array([signal[T-s] if s > 0 else 0.0 for signal, T, s in zip(self.flexSignals,self.flexT,self.flex_state)])

        # Reset instructions for next period.
        self.last_flex_act = self.flex_act
        self.flex_act = zeros(self.N_flex, dtype=integer)
        self.curt_insts = 1e2*ones(self.N_curt)
        self.Q_insts[:] = NaN

    def check_power_setpoint(self, gen, p_setpoint, q_setpoint):
        # Check that the reactive power setpoint is valid for the generator.
        rule = self.pq_rules[self.N_loads+gen]
        try:
            if (q_setpoint > rule.qmax):
                warnings.warn(
                    "Reactive setpoint %g for generator '%s' has been changed to %g (qmax)."
                    % (q_setpoint, rule.description, rule.qmax)
                )
                self.Q_setpoints[self.curtIdInGens.index(gen)] = rule.qmax
                q_setpoint = rule.qmax # for following checks
        except AttributeError:
            pass
        try:
            if (q_setpoint < rule.qmin):
                warnings.warn(
                    "Reactive setpoint %g for generator '%s' has been changed to %g (qmin)."
                    % (q_setpoint, rule.description, rule.qmin)
                )
                self.Q_setpoints[self.curtIdInGens.index(gen)] = rule.qmin
                q_setpoint = rule.qmin # for following checks
        except AttributeError:
            pass
        try:
            lead_slope, lead_off = rule.lead_limit
            if (q_setpoint > lead_slope*self.getPCurtGen(gen)+lead_off):
                new_prod_limit = (q_setpoint-lead_off)/lead_slope
                if new_prod_limit > self.getPCurtGen(gen):
                    warnings.warn(
                        "Active setpoint %g for generator '%s' has been changed to %g (PCurt)."
                        % (new_prod_limit, rule.description, self.getPCurtGen(gen))
                    )
                    new_prod_limit = self.getPCurtGen(gen)
                self.prod_limit[self.curtIdInGens.index(gen)] = new_prod_limit
        except AttributeError:
            pass
        try:
            lag_slope, lag_off = rule.lag_limit
            if (q_setpoint < lag_slope*self.getPCurtGen(gen)+lag_off):
                new_prod_limit = (q_setpoint-lag_off)/lag_slope
                if new_prod_limit > self.getPCurtGen(gen):
                    warnings.warn(
                        "Active setpoint %g for generator '%s' has been changed to %g (PCurt)."
                        % (new_prod_limit, rule.description, self.getPCurtGen(gen))
                    )
                    new_prod_limit = self.getPCurtGen(gen)
                self.prod_limit[self.curtIdInGens.index(gen)] = new_prod_limit
        except AttributeError:
            pass

    def pf_eqs(self,v,y,p,q):
        V = (v[0:self.N_buses]+1j*v[self.N_buses:])
        Sinj = y.conj().dot(V.conj())*V
        return hstack(( Sinj[1:].real-p,
                        Sinj[1:].imag-q,
                        [V[0]-self.V_slack],
                        [V[0].imag])).real


    def current_from_voltage(self, V_buses,k,m,tap,y_line):
        return ((tap*conj(tap))*V_buses[k]-conj(tap)*V_buses[m])*y_line

    def comp_elec_state(self):
        # Aggregate power injections at buses
        P_devices = array([self.getPModLoad(load) for load in range(self.N_loads)]+[self.getPCurtGen(gen) for gen in range(self.N_gens)])
        Q_devices = array([self.getQModLoad(load) for load in range(self.N_loads)]+[self.getQCurtGen(gen) for gen in range(self.N_gens)])
        Pbus = self.dev2bus.dot(P_devices)/self.baseMVA
        Qbus = self.dev2bus.dot(Q_devices)/self.baseMVA

        # Solve power flow equations
        x = optimize.root(self.pf_eqs, self.V_init, args=(self.Y_bus,Pbus,Qbus), method='hybr', options={'xtol' : 1.0e-4}).x
        V_sol = x[0:self.N_buses]+1j*x[self.N_buses:]
        
        # Get the solution
        self.Vmagn = abs(V_sol[1:])
        self.I = array([abs(self.current_from_voltage(V_sol, *self.V_to_I_params[br]))
                        for br in range(self.N_branches)])
        self.S_slack = (self.Y_bus.conj().dot(V_sol.conj())*V_sol)[0]*self.baseMVA

        # Indicate that the electrical state is now computed.
        self.computed = True

    # Mapping from header names to array columns in the input array of devices.
    dheaders = ["bus","Pg","Qg","Qmax","Qmin","Vg","mBase","status","Pmax","Pmin","Pc1","Pc2",
                "Qc1min","Qc1max","Qc2min","Qc2max","ramp_agc","ramp_10","ramp_30","ramp_q","apf"]
    dcol = dict(zip(dheaders,range(len(dheaders))))

class PowerCapabilities(object):

    def __init__(self, description="Unknown device"):
        self.description = description

    @property
    def qp_ratio(self):
        try:
            return self._qp
        except AttributeError:
            raise AttributeError("'qp_ratio' is not set for '%s'." % self.description)
    @qp_ratio.setter
    def qp_ratio(self, value):
        self._qp = value

    @property
    def pmin(self):
        try:
            return self._pmin
        except AttributeError:
            raise AttributeError("'pmin' is not set for '%s'." % self.description)
    @pmin.setter
    def pmin(self, value):
        self._pmin = value

    @property
    def pmax(self):
        try:
            return self._pmax
        except AttributeError:
            raise AttributeError("'pmax' is not set for '%s'." % self.description)
    @pmax.setter
    def pmax(self, value):
        self._pmax = value

    @property
    def qmin(self):
        try:
            return self._qmin
        except AttributeError:
            raise AttributeError("'qmin' is not set for '%s'." % self.description)
    @qmin.setter
    def qmin(self, value):
        self._qmin = value

    @property
    def qmax(self):
        try:
            return self._qmax
        except AttributeError:
            raise AttributeError("'qmax' is not set for '%s'." % self.description)
    @qmax.setter
    def qmax(self, value):
        self._qmax = value

    @property
    def lag_limit(self):
        try:
            return (self._lag_slope, self._lag_offset)
        except AttributeError:
            raise AttributeError("'lag_limit' is not set for '%s'." % self.description)
    @lag_limit.setter
    def lag_limit(self, params):
        assert type(params) is dict
        dQ = params["Qc2min"] - params["Qc1min"]
        dP = params["Pc2"] - params["Pc1"]
        self._lag_slope = dQ/dP
        self._lag_offset = params["Qc1min"] - dQ / dP * params["Pc1"]

    @property
    def lead_limit(self):
        try:
            return (self._lead_slope, self._lead_offset)
        except AttributeError:
            raise AttributeError("'lead_limit' is not set for '%s'." % self.description)
    @lead_limit.setter
    def lead_limit(self, params):
        assert type(params) is dict
        dQ = params["Qc2max"]-params["Qc1max"]
        dP = params["Pc2"]-params["Pc1"]
        self._lead_slope = dQ/dP
        self._lead_offset = params["Qc1max"] - dQ/dP*params["Pc1"]