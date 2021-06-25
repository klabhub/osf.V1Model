%% Adaptation without Plasticity - Psychophysics
% Code to generate figure 1 in
% Quiroga M del M, Morris AP, Krekelberg B. Short-Term Attractive Tilt Aftereffects Predicted by a Recurrent Network Model of Primary Visual Cortex.
% Front. Syst. Neurosci. [Internet] 2019;13:1–14. Available from: https://www.frontiersin.org/article/10.3389/fnsys.2019.00067/full
%
% BK - 2017,2019

%% Setup the model
baseModel  = 'MQ-M';  % Use parameters from Quiroga et al 2017
net = ring(baseModel);
adaptOri        = 20.3906; % Pick one closest to 20, but in the set of simulated neurons to avoid edge effects
testOri         = 0;
adaptDuration   = 200;
testDuration    = 200;
contrast        = 0.5;
%% Simulate
t = 0:adaptDuration+testDuration;
net.stimulusSequence =@(t) adaptationSequence(t,contrast,adaptDuration,testDuration,adaptOri,testOri,net.nrNeurons);
net = simulate(net,'tSpan',[0 adaptDuration+testDuration],'t',t);
% Decode the population using a simple instantaneous vector sum.
[oriPerT] = decode(net,'window',t,'model','vectorsum');
% Simulate adaptation induced rate suppression as it occurs in V1 (~20% with recovery time ~100 ms).
[oriPerTWithSuppression,rateWithSuppression] = decode(net,'window',t,'model','vectorsum','suppression',0.2,'tauRecovery',100,'adaptWindow',[100 200]);
% Simulate adaptation induced rate suppression as it occurs in V1 for longer presentations (i.e. as could
% happen in a typical TAE design (~50% with recovery time ~500 ms).
[oriPerTLongTerm,rateLongTerm] = decode(net,'window',t,'model','vectorsum','suppression',0.5,'tauRecovery',500,'adaptWindow',[100 200]);


%% Show Results
figure(1);
clf
R=2;C=2;
if exist('journalFigure','file')
    [fig,sub]= journalFigure('columns',2,'journal','jov','aspectRatio',1,'panels',{{R,C,[1 2]},{R,C,3},{R,C,4}},'panelLabelDx',-1,'panelLabelDy',0.4,'fontsize',8);
else
    sub(1) = subplot(2,2,[1 2]);
    sub(2) = subplot(2,2,3);
    sub(3) = subplot(2,2,4);
end
nrTimes = 15;
startShow = 5;
stopShow = 300;
lims = [-90  90];
timesToShow = round(linspace(startShow,stopShow,nrTimes));

%% Panel A: time course
axes(sub(1)); hold on;
fromColor = [0.929 0.694 0.125];
toColor  =  [0 0.447 0.741];
for i=1:length(timesToShow)
    color = fromColor + (i-1)*(toColor-fromColor)/length(timesToShow);
    plot(net.preferredOrientation,net.rate(:,timesToShow(i),1),'color',color,'LineWidth',2);
end
xlabel('Preferred Orientation (\circ)')
ylabel('Response (spk/s)')
xlim(lims)
set(gca,'XTick',-90:45:90,'YLim',[ 0 10]);
plot([0 0 ],[0 8],'k--','linewidth',2)

% Time bar
x = 0:1:90;
height = 2;
yOffset  =9;
clrMap = repmat(fromColor',[1 numel(x)]) + repmat(toColor'-fromColor',[1 numel(x)]).*repmat(x./numel(x),[3 1]);
patch([x fliplr(x)],(yOffset-1)+[ones(1,numel(x)) height*ones(1,numel(x))],[x./numel(x) fliplr(x)./numel(x)]);
x0 = adaptDuration/max(timesToShow)*max(x);
h=arrow([x0, yOffset-1.5],[x0 yOffset],40);
h.Color = 'k';
colormap(clrMap')
text(x(1),yOffset,[num2str(startShow) ' ms'],'VerticalAlignment','Top','FontSize',8,'HorizontalAlignment','Right');
text(x(end),yOffset,[num2str(stopShow) ' ms'],'VerticalAlignment','Top','FontSize',8,'HorizontalAlignment','Left')
hold on
text(mean(x),yOffset+0.25*height,'Time since Adapt Onset(ms)','VerticalAlignment','Middle','Fontsize',8,'HorizontalAlignment','Center')

%% Panel B:  Centroid time course
axes(sub(2));hold on;
[~,peakIx] = max(net.rate);
centerOfMass = net.preferredOrientation'*(net.rate./repmat(sum(net.rate),[net.nrNeurons 1]));
for i=1:4:numel(t)
    colorIx = i/stopShow;
    if colorIx>1
        color = [0 0 0];
    else
        color = fromColor + (i./stopShow)*(toColor-fromColor);
    end   
    plot(t(i)-adaptDuration,centerOfMass(i),'.','color',color,'markersize',15);
end
xlabel ('Time since Test Onset (ms)')
ylabel 'Centroid (\circ)'
set(gca,'XLim',[-50 testDuration])

%% Panel C: TAE 
axes(sub(3));hold on
tae  = 2*(oriPerT-sequence(t)); % 2* becuase we compare 20 and -20 adapter
h= plot(t-adaptDuration,tae,'b');
h.LineWidth =1.5;
xlabel ('Time since Test Onset (ms)')
ylabel ('repulsive \leftarrow TAE (\circ) \rightarrow attractive' )

taeWithSuppression = 2*(oriPerTWithSuppression-sequence(t)); % 2* becuase we compare 20 and -20 adapter
plot(t-adaptDuration,taeWithSuppression,'r','linewidth',2);
taeLongTerm = 2*(oriPerTLongTerm-sequence(t)); % 2* becuase we compare 20 and -20 adapter
plot(t-adaptDuration,taeLongTerm,'r:','linewidth',2);
plot(xlim,[0 0],'k:')
set(gca,'XLim',[50 testDuration],'YLIm',[-15 15])
legend('Shift Only','Shift + 20%/100 ms','Shift + 50%/500 ms')

%% Helper function to define a sequence of oriented inputs
function v = adaptationSequence(t,contrast,adaptDuration, testDuration,adaptOri,testOri,nrNeurons)
if t<0
    v = zeros(nrNeurons,1);
elseif t <adaptDuration
    v = ring.singleGratingToInput(adaptOri,contrast,nrNeurons);
elseif t<adaptDuration+testDuration
    v = ring.singleGratingToInput(testOri,contrast,nrNeurons);
else
    v = zeros(nrNeurons,1);
end
end

