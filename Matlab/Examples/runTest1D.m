close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusion.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParams(jsonstruct);

G = cartGrid(100);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusion(paramobj);

doplotgraph = false;
if doplotgraph
    cgt = ComputationalGraphTool(model);
    [g, edgelabels] = cgt.getComputationalGraph();
    plot(g, 'edgelabel', edgelabels);
end

% setup schedule
total = 100;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;;
vols = G.cells.volumes;

initcase = 1;
switch initcase
  case 1
    cA      = zeros(nc, 1);
    cA(1)   = sum(vols);
    cB      = zeros(nc, 1);
    cB(end) = sum(vols);
    cC = zeros(nc, 1);
  case 2
    cA = ones(nc, 1);
    cB = ones(nc, 1);
    cC = zeros(nc, 1);
end

initstate.A.c = cA;
initstate.B.c = cB;
initstate.C.c = cC;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);


%%

ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure

for istate = 1 : numel(states)
    state = states{istate};
    cla
    hold on
    plot(model.G.cells.centroids(:, 1), state.A.c, 'displayname', 'cA');
    plot(model.G.cells.centroids(:, 1), state.B.c, 'displayname', 'cB');
    plot(model.G.cells.centroids(:, 1), state.C.c, 'displayname', 'cC');
    legend show
    drawnow
    pause(0.1);
end
