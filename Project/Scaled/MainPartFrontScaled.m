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
nx = 100; % Grid points in x and y direction
ny = 100;
nxy = 100;
G = cartGrid([nx,ny],[Lxy, Lz]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);


% setup schedule
total = 10*nano*second;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
% For the dimensions to work out in 2d we multiply volume by Lxy
G.cells.volumes = G.cells.volumes .* Lxy;
vols = G.cells.volumes;
nIndex = [nxy/2 + nxy*(nxy-1), nxy/2+1+nxy*(nxy-1)];
rIndex = (1:nxy);

epsilon = 1*nano*meter;
vR = sum(vols(rIndex));
vN = sum(vols(nIndex));
rTot = (1000 / (micro*meter)^2) * (Lxy^2)/Avo; % Total moles of receptors in system
nTot = 5000 / Avo; % Total moles of transmitters in system
r_start = rTot / vR; % Start concentration of receptos
n_start = nTot / vN; % Start concentration of neurotransmitters 

initcase = 1;
switch initcase
  case 1
      cN = zeros(nc,1);
      cN(nIndex) = n_start;
      cR = zeros(nc, 1);
      cR(rIndex) = r_start;
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
framerate = 5 / numel(states);
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
    frame1 = getframe(1);
    frame2 = getframe(2);
    frame3 = getframe(3);

    im1 = frame2im(frame1);
    im2 = frame2im(frame2);
    im3 = frame2im(frame3);

    [imind1,cm1] = rgb2ind(im1,256);
    [imind2,cm2] = rgb2ind(im2,256);
    [imind3,cm3] = rgb2ind(im3,256);

    if istate == 1
         imwrite(imind1,cm1,'Scaledmain2dFrontN.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind2,cm2,'Scaledmain2dFrontR.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
         imwrite(imind3,cm3,'Scaledmain2dFrontC.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
    else
         imwrite(imind1,cm1,'Scaledmain2dFrontN.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind2,cm2,'Scaledmain2dFrontR.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind3,cm3,'Scaledmain2dFrontC.gif','gif','DelayTime',framerate,'WriteMode','append');
    end
    pause(0.01);
    
end
