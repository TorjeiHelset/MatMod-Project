close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionMainPartScaled.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct);
Lxy = 2*220e-9;
Lz = 15e-9;

kappa = paramobj.N.D;
%paramobj.N.D = kappa * (1/Lxy); %Is this necessary?

T = (Lxy^2) / kappa;
Avo = 6.023e23;
% To get concentration per litre we imagine the 2d grid has a height
% epsilon, we let all neurotransmitters reach this level

epsilon = 10e-9;
rho = 10e15;
N0 = 5000/Avo; % Number of moles released in one excitation
R0 = rho * Lxy^2 / Avo; % Number of moles of receptors
% Assuming R0 solved in layer of height epsilon and N0 solved in all of
% the synaptic cleft
nxy = 50;
dxy = Lxy/nxy;

N_max = N0 / (dxy^2); % Concentration if N0 placed in one cell
N_unif = N_max / (nxy^2); % Concentration if N0 uniformly distributed
R_max = R0 / (dxy^2);
R_unif = R_max / (nxy^2); % Concentration if if R0 uniformely distributed


G = cartGrid([nxy, nxy],[Lxy, Lxy]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);


% setup schedule
total = 5e-8;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
vols = G.cells.volumes;

initcase = 4;
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
    case 4
        middle = floor(nxy/2);
        initN = [middle + middle*nxy, middle+1+middle*nxy,...
                 middle+(middle+1)*nxy, middle+1+(middle+1)*nxy];
        cN = zeros(nc,1);
        cN(initN) = N_max/length(initN) .* ones(length(initN),1);
        cR = R_unif .* ones(nc,1);
        cC = zeros(nc,1);
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

    drawnow
    pause(0.01);
    
end
