% Initialize random generator.
s = RandStream('mt19937ar', 'seed', mod(round(now*1e6), 2^32));
RandStream.setGlobalStream(s);

obj = ANM_System();
%uncomment next line for YALMIP to use Ipopt (if installed)
%obj = ANM_System('ipopt');

V = NaN(96, 77);
V(end,:) = obj.getV()';

I = NaN(96, 76);
I(end,:) = obj.getI()'./obj.ratings';

P = zeros(96, 4);
P(end,:) = [sum(obj.getPLoads()) sum(obj.getPModLoads()) sum(obj.getPGens()) sum(obj.getPCurtGens())]';

r = zeros(96,1);

x = -94:1;

model = [];

for i=1:1000
    x = x+1;
    
    obj.transition();
    r(1:end) = [r(2:end); obj.getReward()];
    V(:,:) = [V(2:end,:); obj.getV()'];
    I(:,:) = [I(2:end,:); obj.getI()'./obj.ratings'];
    P(:,:) = [P(2:end,:); [sum(obj.getPLoads()) sum(obj.getPModLoads()) sum(obj.getPGens()) sum(obj.getPCurtGens())]];
    
    clf;
    axes('Position', [0.06 0.52 0.44 0.44]); hold on;
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
    mod = [P(:,2)',fliplr(P(:,1)')];
    set(fill(X,mod,'y', 'DisplayName', 'Mod'),'EdgeColor','None');  
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