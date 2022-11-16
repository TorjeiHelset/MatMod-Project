close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionGliaCells.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParamsGliaCells(jsonstruct);

G = cartGrid([50, 50]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionGliaCells(paramobj);


% setup schedule
total = 200;
n  = 200;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
vols = G.cells.volumes;

initcase = 1;
switch initcase
  case 1
      cN = zeros(nc,1);
      cN(25 + 49*50) = sum(vols)/2;
      cN(26 + 49*50) = sum(vols)/2;
      cR = zeros(nc, 1);
      cR(1:50) = ones(50, 1);
      cC = zeros(nc, 1);
      cT = zeros(nc , 1);
      edges = [(1:48)*50 + 1, (2:49)*50];
      %edges = [(0:48)*50 + 1, (1:49)*50];
      cT(edges) = ones(length(edges), 1);
      cR(edges) = zeros(length(edges),1);
      cN_inac = zeros(nc,1);
      cTb = zeros(nc,1);
end

initstate.N.c = cN;
initstate.R.c = cR;
initstate.C.c = cC;
initstate.T.c = cT;
initstate.Tb.c = cTb;
initstate.N_inac.c = cN_inac;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);



%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4); figure(5); figure(6);

for istate = 1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    plotCellData(model.G, state.N.c);
    colorbar
    title('N concentration')
    
    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.R.c);
    colorbar
    title('R concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.C.c);
    colorbar
    title('C concentration')

    set(0, 'currentfigure', 4);
    cla
    plotCellData(model.G, state.T.c);
    colorbar
    title('T concentration')

    set(0, 'currentfigure', 5);
    cla
    plotCellData(model.G, state.Tb.c);
    colorbar
    title('Tb concentration')

    set(0, 'currentfigure', 6);
    cla
    plotCellData(model.G, state.N_inac.c);
    colorbar
    title('N_{inac} concentration')
    drawnow
    pause(0.1);
    
end
