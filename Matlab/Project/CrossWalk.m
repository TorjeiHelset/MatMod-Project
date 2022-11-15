close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionCrossWalk.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParamsCrossWalk(jsonstruct);

G = cartGrid([50, 50]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionCrossWalk(paramobj);


% setup schedule
total = 100;
n  = 100;
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
      cN(25 + 25*50) = sum(vols)/4;
      cN(26 + 25*50) = sum(vols)/4;
      cN(25 + 24*50) = sum(vols)/4;
      cN(26 + 24*50) = sum(vols)/4;
      cN2 = zeros(nc,1);
      cR = ones(nc, 1);
      cR2 = ones(nc,1);
      cC = zeros(nc, 1);
      cC2 = zeros(nc,1);
      cT = zeros(nc , 1);
      cT2 = zeros(nc,1);
      edges = [1:50, (1:48)*50 + 1, (2:49)*50, (1:50) + 50*49];
      hole = [1 + 23*50, 1 + 24*50, 1+25*50, 1+26*50]; % Hole in glia cells where neurotransmitters can move cross-walk
      hole2 = [(24*50), (25*50), (26*50), (27*50)]-1;

      %edges = [(0:48)*50 + 1, (1:49)*50];
      cT(edges) = ones(length(edges), 1);
      cT(hole) = zeros(4,1);
      cT2(edges) = ones(length(edges),1);
      cT2(hole2) = zeros(4,1);
      cN_inac = zeros(nc,1);
      CN_inac2 = zeros(nc,1);
      cTb = zeros(nc,1);
      cTb2 = zeros(nc,1);
end

initstate.N.c = cN;
initstate.N2.c = cN2;
initstate.R.c = cR;
initstate.R2.c = cR2;
initstate.C.c = cC;
initstate.C2.c = cC2;
initstate.T.c = cT;
initstate.T2.c = cT2;
initstate.Tb.c = cTb;
initstate.Tb2.c = cTb2;
initstate.N_inac.c = cN_inac;
initstate.N_inac2.c = cN_inac;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);



%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4); figure(5); figure(6);
figure(7); figure(8); figure(9); figure(10); figure(11); figure(12);

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

    set(0, 'currentfigure', 7);
    cla
    plotCellData(model.G, state.N2.c);
    colorbar
    title('N2 concentration')
    
    set(0, 'currentfigure', 8);
    cla
    plotCellData(model.G, state.R2.c);
    colorbar
    title('R2 concentration')

    set(0, 'currentfigure', 9);
    cla
    plotCellData(model.G, state.C2.c);
    colorbar
    title('C2 concentration')

    set(0, 'currentfigure', 10);
    cla
    plotCellData(model.G, state.T2.c);
    colorbar
    title('T2 concentration')

    set(0, 'currentfigure', 11);
    cla
    plotCellData(model.G, state.Tb2.c);
    colorbar
    title('Tb2 concentration')

    set(0, 'currentfigure', 12);
    cla
    plotCellData(model.G, state.N_inac2.c);
    colorbar
    title('N_{inac}2 concentration')
    drawnow
    pause(0.001);
    
end
