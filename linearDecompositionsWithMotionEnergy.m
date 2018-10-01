function linearDecompositionsWithMotionEnergy
% linearDecompositionsWithMotionEnergy runs apparent motion stimuli over a
% simple motion energy model, and also runs two different linear
% decompositions of these stimuli through the model, to test for
% differences. 
%
% This code reproduces Figure S6B of the manuscript
% Salazar-Gatzimas E, Agrochao M, Fitzgerald JE, Clark DA (2018) Decorrelation of parallel motion pathways explains the neuronal basis of an illusory motion percept. Current Biology.
%
% This code was tested with Matlab 2015a, but could work with previous
% versions and should work with future versions. The ConfAxis utility is
% provided in the utils/ folder and serves to prettify some of these
% figures.


% Initialize apparent motion stimulus linear decompositions--these
% components will have a positive contrast of 1
% This is the first bar
bar1 = [zeros(1,2),ones(1,8),zeros(1,2)];
% This is the second bar
bar2 = [zeros(1,3),ones(1,7),zeros(1,2)];
% This is the piece of the first bar that appears before the second bar
% appears
bar1InitialComponent = [zeros(1,2),ones(1,1),zeros(1,9)];

% Initialize stimulus to mean background
stimInit = zeros(2,12);
stimSpaceDecomposition2 = stimInit;
stimSpaceDecomposition1 = stimInit;

stimTimeDecomposition1 = stimInit;
stimTimeDecomposition2 = stimInit;

% DECOMPOSITION A: BY SPACE
% When summed, these two stimuli vectors define a displacement from row 1
% (where bar 1 is) to row 2 (where bar 2 is)
stimSpaceDecomposition1(1,:) = bar1;
stimSpaceDecomposition2(2,:) = bar2;
% This function will compute and plot how the motion energy model responds
% to the space decomposition in both directions
modelMotionEnergyResponse(stimSpaceDecomposition1,stimSpaceDecomposition2,'Spatial');

% DECOMPOSITION B: BY TIME
% This is the component of the apparent motion stimulus when both bars are
% present (so both parts look like bar 2, and when summed with the initial
% bar 1 component the same apparent motion stimulus displacing from row 1
% to row 2 gets formed
stimTimeDecomposition1(1,:) = bar1InitialComponent;
stimTimeDecomposition2(1,:) = bar2;
stimTimeDecomposition2(2,:) = bar2;
% This function will compute and plot how the motion energy model responds
% to the time decomposition in both directions
modelMotionEnergyResponse(stimTimeDecomposition1,stimTimeDecomposition2,'Temporal');


function modelMotionEnergyResponse(component1,component2, decompName)
% Assumes (for legend purposes) that the components sum create a down
% stimulus
out = MakeFigure;
decompType = [decompName ' ' 'Decomposition'];
set(out, 'Name', decompType);

% Compute the apparent motion stimulus by adding the two components
motionDirDown = component1+component2;

% Plot summed components (the apparent motion) in one direction
pltH = subplot(3,3,1);
imagesc(motionDirDown); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Full Stim $\downarrow$');
xlabel('time');

% Plot individual components
pltH = subplot(3,3,2);
imagesc(component1); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Component 1 $\downarrow$');

pltH = subplot(3,3,3);
imagesc(component2); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Component 2 $\downarrow$');

% Compute the opposite direction of stimulus displacement
component1OppDir = component1(end:-1:1, :);
component2OppDir = component2(end:-1:1, :);
motionDirUp = component1OppDir + component2OppDir;

% Plot summed components (the apparent motion) in opposite direction
pltH = subplot(3,3,4);
imagesc(motionDirUp); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Full Stim $\uparrow$');
xlabel('time');

% Plot individual components
pltH = subplot(3,3,5);
imagesc(component1OppDir); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Component 1 $\uparrow$');

pltH = subplot(3,3,6);
imagesc(component2OppDir); colormap gray;ConfAxis;
set(pltH, 'CLim', [-1 1]);
title('Component 2 $\uparrow$');

% Compute the motion energy response to the components and the sums in both
% directions, and then average the responses to see DS
motDownResponse = mean(motionEnergyModel(motionDirDown));
motUpResponse = mean(motionEnergyModel(motionDirUp));
motComp1Response = mean(motionEnergyModel(component1));
motComp2Response = mean(motionEnergyModel(component2));
motComp1OppResponse = mean(motionEnergyModel(component1OppDir));
motComp2OppResponse = mean(motionEnergyModel(component2OppDir));

% Plot the motion energy responses
rspAx = subplot(3, 1, 3);
hold on;
bar(1, motDownResponse);
bar(2, motUpResponse);
bar([3 4], [motComp1Response motComp2Response; motComp1OppResponse motComp2OppResponse], 0.8, 'stacked');

% Label the responses
rspAx.XTick = 1:4;
rspAx.XTickLabel = {'M($S_{\downarrow}$)', 'M($S_{\uparrow}$)', 'M($s_{\downarrow,A,1}+s_{\downarrow,A,2}$)', 'M($s_{\uparrow,A,1}+s_{\uparrow,A,2}$)'};
ConfAxis

% We need some latex interpreting to make things look nice
subpans = out.findobj('Type', 'Axes');
xlabs = [subpans.XLabel];
ylabs = [subpans.YLabel];
titles = [subpans.Title];

[xlabs.Interpreter] = deal('latex');
[ylabs.Interpreter] = deal('latex');
[titles.Interpreter] = deal('latex');
[subpans.TickLabelInterpreter] = deal('latex');


function out = motionEnergyModel(stim)
% Split stimulus into two locadtions in space
s1 = stim(1,:);
s2 = stim(2,:);

% Two simple derivative arms, one delayed with respect to the other
f1 = [0,1,-1];
f2 = [1,-1];

temp = filter(f1,1,s1) + filter(f2,1,s2); % Temporal filtering

n = 2; % Squaring nonlinearity
out = temp.^n;

% Utility for figure making
function plotH = MakeFigure(varargin)
plotH = figure('Color',[1 1 1],varargin{:});
set(plotH,'Position',[200,500,1000,1000],'WindowStyle','docked');