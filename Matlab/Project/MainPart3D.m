close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionMainPart.json');
jsonstruct = jsondecode(jsonfile);


paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct);

Lxy = 2*220e-9;
Lz = 15e-9;
kappa = 8e10-7;

T = Lxy^2 / kappa;

Avo =6.023e23;
epsilon = 10e-9;
rho = 10e15;
N_init = 5000 / (Lxy^2 * epsilon * Avo);
R_0 = rho / (epsilon * Avo);

G = cartGrid([50,50,10], [1, 1, Lz/Lxy]);

G = computeGeometry(G);
%figure, plotGrid(G), view(10, 45)
paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);


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

initcase = 3;
switch initcase
  case 1
    cN      = zeros(nc, 1);
    cN(1)   = sum(vols);
    cR      = zeros(nc, 1);
    cR(end) = sum(vols);
    cC = zeros(nc, 1);
  case 2
    cN = ones(nc, 1);
    cR = ones(nc, 1);
    cC = zeros(nc, 1);
  case 3
      cN = zeros(nc,1);
      cN(25 + 25*50) = sum(vols)/4;
      cN(26 + 25*50) = sum(vols)/4;
      cN(25 + 24*50) = sum(vols)/4;
      cN(26 + 24*50) = sum(vols)/4;
      cR = ones(nc, 1);
      cC = zeros(nc, 1);
end

initstate.N.c = cN;
initstate.R.c = cR;
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
    plotCellData(model.G, state.N.c);view(30,60);
    colorbar
    title('N concentration')

    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.R.c);view(30,60);
    colorbar
    title('R concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.C.c);view(30,60);
    colorbar
    title('C concentration')

    drawnow
    pause(0.01);
    
end