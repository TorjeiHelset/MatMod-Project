close all

mrstModule add ad-core mrst-gui 

% Constants used to scale equations
Lxy = 2*220e-9;
Lz = 15e-9;
kappa = 8e-7;

T = Lxy^2 / kappa;

Avo =6.023e23;
epsilon = 10e-9;
rho = 10e15;
N0 = 5000 / (Lxy^2 * epsilon * Avo);
R0 = rho / (epsilon * Avo);

jsonfile = fileread('diffusionMainPart.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParamsMainPartScaled(jsonstruct);
% With T as chosen diffusion coefficient becomes zero, but reaction
% coefficients needs to be updated
paramobj.k1 = paramobj.k1 * T * N0;
paramobj.k_1 = paramobj.k_1 * T;
paramobj.k1_N = paramobj.k1_N * T * R0;
paramobj.k_1_N = paramobj.k_1_N * T * R0/N0;

nx = 50;
ny = 50;
nz = 10;
dz = (Lz/Lxy) / (nz-1);
N_start = nx*ny*epsilon/dz;
R_start = epsilon/dz;


G = cartGrid([nx,ny,nz], [1, 1, Lz/Lxy]);

G = computeGeometry(G);
%figure, plotGrid(G), view(10, 45)
paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPartScaled(paramobj);


% setup schedule
total = 1;
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
      bottom = nx * ny * (nz-1);
      middle = floor(nx*ny*0.5) - 1 + floor(nx*0.5);
      %cN(top + middle) = N_start/4;
      %cN(top + middle + 1) = N_start/4;
      %cN(top + middle+n) = N_start/4;
      %cN(top + middle+n+1) = N_start/4;
      cN(middle) = N_start/4;
      cN(middle + 1) = N_start/4;
      cN(middle+nx) = N_start/4;
      cN(middle+nx+1) = N_start/4;
      cR = zeros(nc, 1);
      cR(bottom+1:bottom+nx*ny) = R_start.*ones(nx*ny,1);
      cC = zeros(nc, 1);

  case 4
      cN = zeros(nc,1);
      cR = zeros(nc,1);
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
