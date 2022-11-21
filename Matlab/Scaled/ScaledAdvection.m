close all

mrstModule add ad-core mrst-gui 

Avo = 6.02214e23; % Avogadro constant

paramobj = AdvectionReactionDiffusionInputParams([]);
paramobj.k1 = 4e6*(mol/litre)*(1/second);
paramobj.k_1 = 5*(1/second);
paramobj.k_t1 = 5e5*(mol/litre)*(1/second);
paramobj.k_t_1 = 0.5*(1/second);
paramobj.k_t2 = 2*(1/second);

paramobj.N.D = 8e-7*(meter^2/second);
paramobj.R.D = 0*(meter^2/second);
paramobj.C.D = 0*(meter^2/second);
paramobj.T.D = 0*(meter^2/second);
paramobj.Tb.D = 0*(meter^2/second);

Lxy = 2*0.22*micro*meter;
Lz = 0.15*micro*meter;
nx = 100; % Grid points in x and y direction
ny = 100;
nxy = 100;
G = cartGrid([nx,ny],[Lxy, Lxy]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = AdvectionReactionDiffusion(paramobj);


% setup schedule
total = 10*nano*second;;
n  = 100;
dt = total/n;
step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));

control.none = [];
schedule = struct('control', control, 'step', step);

% setup initial state

nc = G.cells.num;
% For the dimensions to work out in 2d we need to imagine everything to
% have a height epsilon. We set epsilon = 1nm;
epsilon = 1*nano*meter; 
G.cells.volumes = G.cells.volumes .* epsilon;
vols = G.cells.volumes;
rIndex = (1:nc);
nIndex = [nxy/2+nxy*(nxy/2), nxy/2+nxy*(nxy/2)+1,...
          nxy/2+nxy*(1+nxy/2), nxy/2+nxy*(1+nxy/2)+1]';

vR = sum(vols(rIndex));
vN = sum(vols(nIndex)) * (Lz/epsilon);
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
      cT = zeros(nc , 1);
      edges = [1:nxy, (1:(nxy-2))*nxy + 1, (2:(nxy-1))*nxy,...
              (1:nxy) + nxy*(nxy-1)];
      cT(edges) = ones(length(edges), 1);
      hole = [1 + (ny/2-2)*nx, 1 + (ny/2-1)*nx, 1+(ny/2)*nx, 1+(ny/2+1)*nx];
      cT(hole) = zeros(4,1);
      cTb = zeros(nc,1);
end

initstate.N.c = cN;
initstate.R.c = cR;
initstate.C.c = cC;
initstate.T.c = cT;
initstate.Tb.c = cTb;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);

%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4); figure(5);
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

    drawnow
    frame1 = getframe(1);
    frame2 = getframe(2);
    frame3 = getframe(3);
    frame4 = getframe(4);
    frame5 = getframe(5);

    im1 = frame2im(frame1);
    im2 = frame2im(frame2);
    im3 = frame2im(frame3);
    im4 = frame2im(frame4);
    im5 = frame2im(frame5);

    [imind1,cm1] = rgb2ind(im1,256);
    [imind2,cm2] = rgb2ind(im2,256);
    [imind3,cm3] = rgb2ind(im3,256);
    [imind4,cm4] = rgb2ind(im4,256);
    [imind5,cm5] = rgb2ind(im5,256);

    if istate == 1
         imwrite(imind1,cm1,'ScaledAdvectionN.gif','gif', 'DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind2,cm2,'ScaledAdvectionR.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind3,cm3,'ScaledAdvectionC.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
         imwrite(imind4,cm4,'ScaledAdvectionT.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind5,cm5,'ScaledAdvectionTb.gif','gif','DelayTime',framerate, 'Loopcount',inf);
    else
         imwrite(imind1,cm1,'ScaledAdvectionN.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind2,cm2,'ScaledAdvectionR.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind3,cm3,'ScaledAdvectionC.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind4,cm4,'ScaledAdvectionT.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind5,cm5,'ScaledAdvectionTb.gif','gif','DelayTime',framerate,'WriteMode','append');
    end
    pause(0.001);
    
end
