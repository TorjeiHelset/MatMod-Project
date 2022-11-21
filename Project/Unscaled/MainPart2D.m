close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionMainPart.json');
jsonstruct = jsondecode(jsonfile);

paramobj = ReactionDiffusionInputParamsMainPart(jsonstruct);

nx = 100;
ny = 100;
G = cartGrid([nx, ny]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = ReactionDiffusionMainPart(paramobj);

% setup schedule
total = 1;
n  = 100;
dt = total/n;
%step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
step  = struct('val', dt*(1:n), 'control', ones(n, 1));

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
         imwrite(imind1,cm1,'main2dN.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind2,cm2,'main2dR.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
         imwrite(imind3,cm3,'main2dC.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
    else
         imwrite(imind1,cm1,'main2dN.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind2,cm2,'main2dR.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind3,cm3,'main2dC.gif','gif','DelayTime',framerate,'WriteMode','append');
    end

    pause(0.01);
    
end
