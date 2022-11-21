close all

mrstModule add ad-core mrst-gui 


paramobj       = ReactionDiffusionInputParamsMainPart([]);
paramobj.k1   = 4e6*(mol/litre)*(1/second);
paramobj.k_1   = 5*(1/second);
paramobj.N.D   = 8e-7*(meter^2/second);
paramobj.R.D   = 0*(meter^2/second);
paramobj.C.D = 0*(meter^2/second);

Lxy = 2*0.22*micro*meter;
nxy = 100; % Grid points in x and y direction
G = cartGrid([nxy,nxy],[Lxy, Lxy]);
G = computeGeometry(G);

paramobj.G = G;
paramobj = paramobj.validateInputParams();
model = ReactionDiffusionMainPart(paramobj);

% setup schedule
total = 10*nano*second;
n     = 100;
dt    = total/n;
step  = struct('val', dt*(1:n), 'control', ones(n, 1));
%step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

G = model.G;
Rindex = (1:nxy^2)';    
Nindex = [nxy/2+nxy*(nxy/2), nxy/2+nxy*(nxy/2)+1,...
          nxy/2+nxy*(1+nxy/2), nxy/2+nxy*(1+nxy/2)+1]';

% setup initial state

A = 6.02214076e23; % Avogadro constant

nc     = G.cells.num;
vols   = G.cells.volumes;
Vr = sum(vols(Rindex));
rTot = (1000 / (micro*meter)^2) * (Lxy^2); % Total amount of receptors
%initR = (rTot/A*Vr);
initR = rTot/(A*length(Rindex));

Vn = sum(vols(Nindex));
nTot = 5000; % Released in one excitation
%initN = (nTot/(A*Vn));
initN = 5000 /(A*length(Nindex));

initcase = 1;

switch initcase
  case 1
    cR = zeros(nc, 1);
    cR(Rindex) = initR;
    cN = zeros(nc, 1);
    cN(Nindex)  = initN;
    cC = zeros(nc, 1);
  case 2
    cR   = ones(nc, 1);
    cN   = ones(nc, 1);
    cC = zeros(nc, 1);
end

initstate.R.c   = cR;
initstate.N.c   = cN;
initstate.C.c = cC;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);


% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3);

for istate = 1 : numel(states)

    state = states{istate};

    set(0, 'currentfigure', 1);
    cla
    plotCellData(model.G, state.R.c);
    colorbar
    title('R concentration')
    
    set(0, 'currentfigure', 2);
    cla
    plotCellData(model.G, state.N.c);
    colorbar
    title('N concentration')

    set(0, 'currentfigure', 3);
    cla
    plotCellData(model.G, state.C.c);
    colorbar
    title('C concentration')

    drawnow
    pause(0.1);
    
end