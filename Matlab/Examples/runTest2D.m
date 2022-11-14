close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusion2.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParams(jsonstruct);

G = cartGrid([50, 50]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusion(paramobj);


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
    cA      = zeros(nc, 1);
    cA(1)   = sum(vols);
    cB      = zeros(nc, 1);
    cB(end) = sum(vols);
    cC = zeros(nc, 1);
  case 2
    cA = ones(nc, 1);
    cB = ones(nc, 1);
    cC = zeros(nc, 1);
  case 3
      cA = zeros(nc,1);
      cA(25 + 25*50) = sum(vols)/4;
      cA(26 + 25*50) = sum(vols)/4;
      cA(25 + 24*50) = sum(vols)/4;
      cA(26 + 25*50) = sum(vols)/4;
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

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3);

for istate = 1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    plotCellData(model.G, state.A.c);
    colorbar
    title('A concentration')
    
    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.B.c);
    colorbar
    title('B concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.C.c);
    colorbar
    title('C concentration')

    drawnow
    pause(0.1);
    
end
