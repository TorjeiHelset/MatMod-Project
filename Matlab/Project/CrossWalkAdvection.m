close all

mrstModule add ad-core mrst-gui 

jsonfile = fileread('diffusionGliaCells.json');
jsonstruct = jsondecode(jsonfile);

paramobj = AdvectionReactionDiffusionInputParams(jsonstruct);

G = cartGrid([100, 50]);
G = computeGeometry(G);

paramobj.G = G;

paramobj = paramobj.validateInputParams();

model = AdvectionReactionDiffusion(paramobj);


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
      cN = zeros(nc,1);
      cN(70 + 25*100) = sum(vols)/4;
      cN(76 + 25*100) = sum(vols)/4;
      cN(75 + 24*100) = sum(vols)/4;
      cN(76 + 24*100) = sum(vols)/4;
      cR = ones(nc, 1);
      cC = zeros(nc, 1);
      cT = zeros(nc , 1);
      edges = [1:100, (1:48)*100 + 1,(1:48)*100 + 50, (1:48)*100 + 51, (2:49)*100, (1:100) + 100*49, ];
      %edges = [(0:48)*50 + 1, (1:49)*50];
      hole = [(24:27)*100, (24:27)*(100)+1] + 50;
      cT(edges) = ones(length(edges), 1);
      cT(hole) = zeros(length(hole),1);
      %cR(edges) = zeros(length(edges),1);
      cN_inac = zeros(nc,1);
      cTb = zeros(nc,1);
end

initstate.N.c = cN;
initstate.R.c = cR;
initstate.C.c = cC;
initstate.T.c = cT;
initstate.Tb.c = cTb;
initstate.N_inac.c = cN_inac;

% run simulation

nls = NonLinearSolver();
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'NonLinearSolver', nls);



%%

% Remove empty states (could have been created if solver did not converge)
ind = cellfun(@(state) ~isempty(state), states);
states = states(ind);

figure(1); figure(2); figure(3); figure(4); figure(5); figure(6);
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

    set(0, 'currentfigure', 6);
    cla
    plotCellData(model.G, state.N_inac.c);
    colorbar
    title('N_{inac} concentration')
    drawnow
    frame1 = getframe(1);
    frame2 = getframe(2);
    frame3 = getframe(3);
    frame4 = getframe(4);
    frame5 = getframe(5);
    frame6 = getframe(6);

    im1 = frame2im(frame1);
    im2 = frame2im(frame2);
    im3 = frame2im(frame3);
    im4 = frame2im(frame4);
    im5 = frame2im(frame5);
    im6 = frame2im(frame6);

    [imind1,cm1] = rgb2ind(im1,256);
    [imind2,cm2] = rgb2ind(im2,256);
    [imind3,cm3] = rgb2ind(im3,256);
    [imind4,cm4] = rgb2ind(im4,256);
    [imind5,cm5] = rgb2ind(im5,256);
    [imind6,cm6] = rgb2ind(im6,256);

    if istate == 1
         imwrite(imind1,cm1,'CW_AdvectionN.gif','gif', 'DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind2,cm2,'CW_AdvectionR.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind3,cm3,'CW_AdvectionC.gif','gif', 'DelayTime',framerate,'Loopcount',inf);
         imwrite(imind4,cm4,'CW_AdvectionT.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind5,cm5,'CW_AdvectionTb.gif','gif','DelayTime',framerate, 'Loopcount',inf);
         imwrite(imind6,cm6,'CW_AdvectionN_inac.gif','gif','DelayTime',framerate, 'Loopcount',inf);
    else
         imwrite(imind1,cm1,'CW_AdvectionN.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind2,cm2,'CW_AdvectionR.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind3,cm3,'CW_AdvectionC.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind4,cm4,'CW_AdvectionT.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind5,cm5,'CW_AdvectionTb.gif','gif','DelayTime',framerate,'WriteMode','append');
         imwrite(imind6,cm6,'CW_AdvectionN_inac.gif','gif','DelayTime',framerate,'WriteMode','append');
    end
    pause(0.01);
    
end
