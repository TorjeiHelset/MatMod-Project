close all

mrstModule add ad-core mrst-gui 

%jsonfile = fileread('diffusionMainPartScaled.json');
%jsonstruct = jsondecode(jsonfile);
%paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct);

paramobj = ReactionDiffusionInputParamsMainPart([]);
paramobj.k1   = 4e6*(mol/litre)*(1/second);
paramobj.k_1   = 5*(1/second);
paramobj.N.D   = 8e-7*(meter^2/second);
paramobj.R.D   = 0*(meter^2/second);
paramobj.C.D = 0*(meter^2/second);


Lxy = 2*220e-9;
Lz = 15e-9;


kappa = paramobj.N.D;
paramobj.N.D = kappa * (1/Lxy^2);

% Using physical quantities
% Gradient operator does not take into account steplenghts so we multiply
% steplenghts into diffusion coefficient

% Constants used to scale equations
T = (Lxy^2) / kappa;

Avo = 6.023e23;
epsilon = 10e-9;
rho = 10e15;
N0 = 5000 /Avo; % Number of moles released in one excitation
R0 = rho * Lxy^2 / Avo; % Number of moles of receptors


nxy = 50;
nz = 10;
dxy = Lxy/nxy;
dz = Lz/nz;

N_max = N0 / (dxy^2 * dz); % Concentration if N0 placed in one cell
N_unif = N_max / (nxy^2);
R_max = R0 / (dxy^2 * dz);
R_unif = R_max / (nxy^2); % Concentration if if R0 uniformely distributed in one layer

G = cartGrid([nxy,nxy,nz],[Lxy, Lxy, Lz]); % z-scale not fully correct

G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);


% setup schedule
total = 5e-2; % Only simulate till timescale we found
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state
bottom = nxy * nxy * (nz-1); % First cell on bottom layer
middle = floor(nxy*nxy*0.5) - 1 + floor(nxy*0.5); % Middle cell of one layer
indexN = [middle-nxy, middle-nxy+1,middle+2*nxy, middle+2*nxy+1,...
       middle-1, middle, middle+1,middle+2,...
       middle+nxy-1, middle+nxy,middle+nxy+1,middle+nxy+2];

% receptors at bottom, transmitters at initN

nc     = G.cells.num;
vols   = G.cells.volumes;
initCR = 1000*((micro*meter)^2)/sum(G.cells.volumes(1:bottom)); %fix this according to dimension
initN = (5000) / sum(G.cells.volumes(indexN));

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
      bottom = nxy * nxy * (nz-1); % First cell on bottom layer
      middle = floor(nxy*nxy*0.5) - 1 + floor(nxy*0.5); % Middle cell of one layer

      %cN(1:nxy^2) =  (N_max/(nxy^2)).*ones(nxy^2,1);
      initN = [middle-nxy, middle-nxy+1,middle+2*nxy, middle+2*nxy+1,...
               middle-1, middle, middle+1,middle+2,...
               middle+nxy-1, middle+nxy,middle+nxy+1,middle+nxy+2];
      cN(initN) = (N_max/length(initN)) .* ones(length(initN),1);
      %cN(middle) = N_max/4;
      %cN(middle + 1) = N_max/4;
      %cN(middle+nxy) = N_max/4;
      %cN(middle+nxy+1) = N_max/4;
      cR = zeros(nc, 1);
      cR(bottom+1:bottom+nxy^2) = R_unif.*ones(nxy^2,1);
      cC = zeros(nc, 1);

  case 4
      cN = zeros(nc,1);
      cN(indexN) = initN;
      cR = zeros(nc,1);
      cR(1:bottom) = initCR;
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
