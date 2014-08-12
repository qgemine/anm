function [actions, model] = example_policy(sys, model)

    % Clone the object of the system to keep it unchanged.
    obj = sys.copy();    
    % planning horizon
    T = 15;
    % Number of scenarios (i.e. of clusters) to build from trajectories.
    N_scens = 5;
    
    %% Build model of the mathematical program.
    
    if isempty(model)
        display('Building optimization model (done once per simulation)...');
        display(' * this step might last from a few seconds to a few minutes * ');
        
        % Number of nodes in the scenarios tree.
        N_nodes = 1 + T * N_scens; % the constant '1' accounts for the root node.
        % Variables representing the power curtailed at generators.
        delta_curt = cell(N_nodes,1);
        % Variables representing the margin between generation limits and
        % power generation.
        delta_margin = cell(N_nodes,1);
        % Variables representing the activation of flexible loads.
        act = binvar(obj.N_loads,1);
        % Parameters of forecasted power injections.
        P_forecast = cell(N_nodes,1);
        % Parameters of values of modulation signals.
        P_modulation = sdpvar(obj.N_loads,T, 'full');
        % Parameters indicating the flexibility services currently running
        on = sdpvar(obj.N_loads,1);
        % Parameters of the cost of curtailement.
        curt_cost = sdpvar(T,1);
        % Parameters of the probability of occurence of each scenario.
        prob_scens = sdpvar(N_scens,1);
        
        % Number of action vectors (= number of edges)
        N_actions = (T-1)*N_scens+1;
        % Variables of the instructions of maximal power injection of generators.
        P_max = cell(N_actions,1);

        for i = 2:N_nodes
            % Active power injection forecast for node i.
            P_forecast{i} = sdpvar(obj.N_gens+obj.N_loads, 1);
            % Power curtailment variables for node i.
            delta_curt{i} = sdpvar(obj.N_gens, 1);
            % Curtailement margin variables for node i.
            delta_margin{i} = sdpvar(obj.N_gens, 1);
        end

        for i = 1:N_actions
            % Curtailment instruction variables for node i.
            P_max{i} = sdpvar(obj.N_gens, 1);
        end
        
        % Indices of generators in vectors of device variables.
        gens_range = (obj.N_loads+1):(obj.N_gens+obj.N_loads);
        % Indices of loads in vectors of device variables.
        loads_range = 1:obj.N_loads;

        cost = sum((15*abs(obj.dPnom)).*(act-on));
        constraints = [];

        for t=1:T
            for i = 1:N_scens
                % Node ID in the tree.
                node = 1+(t-1)*N_scens+i;
                % Node ID of the closest parent in the tree.
                node_a = max(1,1+(t-2)*N_scens+i);
                % Cost for curtailement instructions.
                cost = cost + prob_scens(i) * ( curt_cost(t)*sum(delta_curt{node}) + 1e-3*sum(delta_curt{node}.^2) - 1e-2*sum(delta_margin{node}) + 1e-3*sum(delta_margin{node}.^2));

                % delta_P >= max(P_forecast-P_curt, 0)
                constraints = [constraints, delta_curt{node} >= 0];
                constraints = [constraints, delta_curt{node} >= P_forecast{node}(gens_range) - P_max{node_a}];
                
                % delta_margin >= max(P_curt-P_forecast, 0)
                constraints = [constraints, delta_margin{node} >= 0];
                constraints = [constraints, delta_margin{node} >= P_max{node_a} - P_forecast{node}(gens_range)];
                
                % Curtailment instructions must be positive and do not need
                % to be higher than 10 MW.
                constraints = [constraints, zeros(obj.N_gens,1) <= P_max{node_a} <= 10*ones(obj.N_gens,1)];
                
                % Flexible loads that are already active must stay active.
                constraints = [constraints, act >= on];
                
                % Crude approximation of operational constraints.
                constraints = [constraints, sum(P_max{node_a}) + sum( P_forecast{node}(loads_range) + act .* P_modulation(:,t) ) <= 18];

            end
        end

        % Build the optimization model, which produces and solves a mathematical
        % program as a function of some parameters (e.g. forecast of power
        % injections).
        model = optimizer(constraints, cost, sdpsettings('verbose', 0, 'solver', 'gurobi'), [P_forecast(2:end); {P_modulation}; {on}; {curt_cost}; {prob_scens}] , [P_max(1); {act}]);
    end
    
    %% Generate scenarios of stochastic power injections.
    
    % Number of trajectories to generate 
    N_samples = 250;

    N_dev = obj.N_gens+obj.N_loads;
    timeseries = zeros(N_samples, T*N_dev);

    display('Sampling trajectories of the system...');
    % Generate trajectories
    for i = 1:N_samples
        inst = obj.copy();
        for t=1:T
            inst.transition();
            timeseries(i, ((t-1)*N_dev+1):(t*N_dev)) = [inst.getPLoads()' inst.getPGens()'];
        end
    end
    
    %% Aggregate scenarios into clusters.

    % Perform k-means clustering
    idx = [];
    restart = 1;
    display('Clustering trajectories...');
    while restart
        try
            idx = kmeans(timeseries,N_scens);
            restart = 0;
        catch
            restart = 1;
        end
    end
    clusters = cell(N_scens,1);
    scens = cell(N_scens,1);
    probs = zeros(N_scens,1);
    for i = 1:N_scens
        clusters{i} = timeseries(idx==i,:);
        
        % Compute mean scenario for the cluster
        if sum(idx==i) == 1 % special case when there is only a single timeserie.
            scens{i} = reshape(clusters{i}, (obj.N_gens+obj.N_loads), T)';
        else
            scens{i} = reshape(mean(clusters{i}), (obj.N_gens+obj.N_loads), T)';
        end
        
        % Probability of occurence of the mean scenario of the cluster,
        % estimated using the samples that have been drawn.
        probs(i) = size(clusters{i},1)/N_samples;
    end
    
    
    
    %%
    % Build parameters of the optimization model as a function of the
    % system's state.
    
    % Parameters of forecasted power injections.
    P_forecast = cell(T*N_scens,1);
    
    % Parameters of modulation signals of loads.
    Tmin = max([max(obj.Tflex) T]);
    rel_time = repmat(1:Tmin, obj.N_loads, 1);
    Tflex_mat = repmat(obj.Tflex, 1, Tmin);
    dPnom_mat = repmat(obj.dPnom, 1, Tmin);
    P_modulation = (rel_time <= Tflex_mat) .* dPnom_mat .* sin( 1.8*pi*(rel_time-0.5*(Tflex_mat+1))./(Tflex_mat-1) );
    is_active = obj.getFlexState() > 0;
    flex_shift = (obj.Tflex - obj.getFlexState() + 1);
    
    for l = 1:obj.N_loads
        if is_active(l)
            P_modulation(l, :) = [P_modulation(l, flex_shift(l):end) zeros(1,flex_shift(l)-1)];
        end
    end
    
    % Forecasts of power injections
    for t=1:T
        for i = 1:N_scens
            node = (t-1)*N_scens+i; % Node ID in the tree.
            P_forecast{node} = scens{i}(t,:)';
        end
    end
    
    % Parameters of curtailment prices.
    curt_cost = obj.getCurtPrice(mod(obj.getQuarter():(obj.getQuarter()+T-1),96)+1)/4;
    
    % Solve the mathematical program using current parameters and get
    % control actions.
    display('Solving mathematical program...');
    actions = model{[P_forecast; {P_modulation(:,1:T)}; {obj.getFlexState() > 0}; {curt_cost}; {probs}]};
    actions{2} = actions{2}.*(obj.getFlexState()==0);
    
end