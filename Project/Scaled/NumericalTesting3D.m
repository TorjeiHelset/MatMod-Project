close all

mrstModule add ad-core mrst-gui 

Avo = 6.02214e23; % Avogadro constant

paramobj = ReactionDiffusionInputParamsMainPart([]);
paramobj.k1 = 4e6*(mol/litre)*(1/second);
paramobj.k_1 = 5*(1/second);
paramobj.N.D = 8e-7*(meter^2/second);
paramobj.R.D = 0*(meter^2/second);
paramobj.C.D = 0*(meter^2/second);

Lxy = 2*0.22*micro*meter;
Lz = 0.15*micro*meter;
nxy = 50;
nz = 10;
G = cartGrid([nxy,nxy, nz],[Lxy, Lxy, Lz]);
G = computeGeometry(G);


paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);

% setup schedule
total = 10*nano*second;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
%step  = struct('val', dt*(1:n), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
vols = G.cells.volumes;

rIndex = ((nxy^2)*(nz-1)+1:nxy*nxy*nz)';
nIndex = [nxy/2+nxy*(nxy/2), nxy/2+nxy*(nxy/2)+1,...
          nxy/2+nxy*(1+nxy/2), nxy/2+nxy*(1+nxy/2)+1]';


initcase = 1;
switch initcase
  case 1
    vR = sum(vols(rIndex));
    vN = sum(vols(nIndex));
    rTot = (1000 / (micro*meter)^2) * (Lxy^2)/Avo; % Total moles of receptors in system
    nTot = 5000 / Avo; % Total moles of transmitters in system
    r_start = rTot / vR; % Start concentration of receptos
    n_start = nTot / vN; % Start concentration of neurotransmitters 

    cN      = zeros(nc, 1);
    cN(nIndex)   = n_start;
    cR      = zeros(nc, 1);
    cR(rIndex) = r_start;
    cC = zeros(nc, 1);
  case 2
    % Testing how much must be transmitted for a signal to be sent
    vR = sum(vols(rIndex));
    vN = sum(vols(nIndex));
    rTot = (1000 / (micro*meter)^2) * (Lxy^2)/Avo; % Total moles of receptors in system
    nTot = 96 / Avo; % Total moles of transmitters in system
    r_start = rTot / vR; % Start concentration of receptos
    n_start = nTot / vN; % Start concentration of neurotransmitters 

    cN      = zeros(nc, 1);
    cN(nIndex)   = n_start;
    cR      = zeros(nc, 1);
    cR(rIndex) = r_start;
    cC = zeros(nc, 1);
  case 3
      cN = zeros(nc,1);
      cN(nx/2 + ny*nx/2) = sum(vols)/4;
      cN(nx/2+1+ny*nx/2) = sum(vols)/4;
      cN(nx/2 + ny*(nx/2+1)) = sum(vols)/4;
      cN(nx/2 +1+ ny*(nx/2+1)) = sum(vols)/4;
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

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

nbeginning = sum(cN);
rbeginning = sum(cR);
cbeginning = sum(cC);
rsignal = rbeginning / 2;

switch initcase
  case 1
      figure();
      fraction = zeros(numel(states),1);
       for istate = 1 : numel(states)
           state = states{istate};
           fraction(istate) = sum(state.R.c) / rbeginning;
       end

       % Plotting fraction of bounded receptors as function of time
       timesteps = (1:numel(states)) * dt;
       plot(timesteps, fraction)
       hold on;
       plot([0, numel(states)*dt], [0.5, 0.5], 'r--')
       hold off;
       title("Fraction of receptors left in the system as a function of time")
       xlabel("Time [s]")

   case 2
       for istate = 1 : numel(states)
        state = states{istate};
        if (sum(state.R.c) < rsignal)
            X = ['Signal transmitted after ',num2str(istate*dt),' seconds.'];
            disp(X)
            disp(istate*dt)
            break;
        end
      end
end

