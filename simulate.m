function rewards = simulate(policy, n_steps, figures, randstream)

% Initialize random generator or use to one provided as parameter.
if nargin >= 4
    RandStream.setGlobalStream(randstream);
else
    RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', mod(round(now*1e6), 2^32)));
end

% Check if a custom time horizon has been specified.
T = 192;
if nargin >= 2
    T = n_steps;
end

% Initialize output vector (vector of observed rewards).
rewards = zeros(T,1);

% Get an instance of the ANM system.
syst = ANM_System();

% Initialize matrices of quantities to plot.
V = NaN(96, 77);
I = NaN(96, 76);
P = zeros(96, 4);
r = zeros(96,1);

% Get current state of the system.
V(end,:) = syst.getV()';
I(end,:) = syst.getI()'./syst.ratings';
P(end,:) = [sum(syst.getPLoads()) sum(syst.getPModLoads()) sum(syst.getPGens()) sum(syst.getPCurtGens())]';

% Initial x-axis values for plots.
x = -94:1;

% Data structure of the policy (optional) is initially empty.
model = [];

% Check if plots have to be displayed.
disp_figs = 0;
if (nargin >= 3) && (figures ~= 0)
    disp_figs = 1;
end
    
if disp_figs
    figure('units','normalized','outerposition',[0 0 1 1]);
end

% Run the simulation
for t=1:T
    
    % Get control actions using the policy.
    [actions, model] = policy(syst, model);
    syst.setCurtailment(max(actions{1},0));
    syst.setFlexAct(round(actions{2}));
    % Trigger a transition of the system.
    syst.transition();
    % Get electrical quantities for the new state and last reward.
    r(1:end) = [r(2:end); syst.getReward()];
    V(:,:) = [V(2:end,:); syst.getV()'];
    I(:,:) = [I(2:end,:); syst.getI()'./syst.ratings'];
    P(:,:) = [P(2:end,:); [sum(syst.getPLoads()) sum(syst.getPModLoads()) sum(syst.getPGens()) sum(syst.getPCurtGens())]];
    
    rewards(t) = r(end);
    display(['Time step ' int2str(t) ' (reward = ' num2str(rewards(t)) ').']);
    
    if disp_figs
        clf;
        x = x+1;

        axes('Position', [0.06 0.52 0.44 0.44]);  hold on;
        title('Evolution of voltage magnitudes (in p.u.)', 'FontSize', 18, 'Interpreter','latex');
        plot(x, 1.05*ones(size(x)), 'r--');
        plot(x, 0.95*ones(size(x)), 'r--');
        plot(x,V(:,2:end), 'k');
        xlim([x(1) x(end)]);
        ylim([0.94 1.06]);
        set(gca, 'XTickLabel', {'     -24h', '-12h', 'Real-time              '},'XTick',[x(1) x(48) x(end)], 'TickDir', 'out', 'FontSize', 13)
        box on;

        axes('Position', [0.52 0.52 0.44 0.44]); hold on;
        title('Evolution of current magnitudes (in \% of link ratings)', 'FontSize', 18, 'Interpreter','latex');
        plot(x, 100*ones(size(x)), 'r--');
        plot(x,100.0*I, 'k');
        xlim([x(1) x(end)]);
        ylim([0 120]);
        set(gca, 'XTickLabel', {'     -24h', '-12h', 'Real-time              '},'XTick',[x(1) x(48) x(end)], 'TickDir', 'out', 'YAxisLocation', 'right', 'FontSize', 13)
        box on;

        axes('Position', [0.06 0.01 0.44 0.44]); hold on;
        title('Evolution of power injections and withdrawals (in MW)', 'FontSize', 18, 'Interpreter','latex');
        X = [x,fliplr(x)];
        curt = [P(:,4)',fliplr(P(:,3)')];
        set(fill(X,curt,'r', 'DisplayName', 'Curt'),'EdgeColor','None');
        modulation = [P(:,2)',fliplr(P(:,1)')];
        set(fill(X,modulation,'y', 'DisplayName', 'Mod'),'EdgeColor','None');
        plot(x, zeros(size(x)), 'k-.', 'HandleVisibility','off');
        plot(x, P(:,2), 'Color', 'b', 'DisplayName', 'Loads');
        plot(x, P(:,1), 'b:', 'HandleVisibility','off');
        plot(x, P(:,4), 'Color', 'k', 'DisplayName', 'DGs');
        plot(x, P(:,3), 'k:', 'HandleVisibility','off');
        xlim([x(1) x(end)]);
        ylim([min(-22.5,min(P(:,2))) max(30,max(P(:,3)))]);
        set(gca, 'XTickLabel', [],'XTick', [], 'TickDir', 'out', 'FontSize', 12)
        legend('Location','NorthWest');
        box on;

        axes('Position', [0.52 0.01 0.44 0.44]); hold on;
        title('Evolution of instantaneous rewards ( $-log(-r_{t}+1)$ )', 'FontSize', 18, 'Interpreter','latex');
        [xx,yy] = stairs(x, -log(-r+1));
        patch([x(1); xx; x(end)], [0; yy; 0], 'k');
        xlim([x(end)-95 x(end)]);
        ylim([-12.5 1]);
        set(gca, 'XTickLabel', [],'XTick', [], 'TickDir', 'out', 'YAxisLocation', 'right', 'FontSize', 12)
        box on;

        drawnow;
    end
end

end