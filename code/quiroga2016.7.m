%% Adaptation without Plasticity - Mechanisms
%
%  Code to generate figure 7 in
%  Quiroga M del M, Morris AP, Krekelberg B. Adaptation without plasticity. Cell Rep. 2016;17:58–68.
%  Available from: http://dx.doi.org/10.1016/j.celrep.2016.08.089
%
% BK -May 2016

baseModel = 'MQ-M';
simulationDuration = 175; %
t = 1:simulationDuration;
stimDuration =50;
adapter = -25;
limits = true;      % Set to true to adjust limits that match the example in the paper.
tPoint = 65;         % Show the state at this time
testOrientation = []; % [] means use shifted preferred. Otherwise set a value here.
shrink = [1 1];
adapterContrast = 0.5;
testContrast = 0.5;
%% Simulate the adaptation protocol

a = ring(baseModel);
a.shrinkFactor = shrink;

testOri = -90:1:90;
corticalInputPerOri  = nan(a.nrNeurons, numel(t),numel(testOri));
lgnInputPerOri= nan(a.nrNeurons, numel(t),numel(testOri));
unadaptedCorticalInputPerOri= nan(a.nrNeurons, numel(t),numel(testOri));
[~,zeroNeuronIx] = min(abs(a.preferredOrientation-0));
zeroNeuronPO  = a.preferredOrientation(zeroNeuronIx);
rate= nan(a.nrNeurons, numel(t),numel(testOri));
popPeak = nan(numel(t),numel(testOri));
cntr= 0;
for ori=testOri
    cntr = cntr+1;
    a = ring(baseModel); % Start fresh - adapted
    a.shrinkFactor = shrink;
    a.stimulusSequence = ring.adaptSequence(t,adapter, ori, stimDuration, simulationDuration,0,adapterContrast,testContrast);
    a = simulate(a,'tSpan',[0 max(t)],'t',t);
    
    
    u  =ring(baseModel); % Start fresh - unadapted
    u.shrinkFactor = shrink;
    u.stimulusSequence = ring.adaptSequence(t,adapter, ori, stimDuration, simulationDuration,0,0,testContrast);
    u = simulate(u,'tSpan',[0 max(t)],'t',t);
    
    % Collect
    unadaptedCorticalInputPerOri(:,:,cntr) = u.corticalInput(t);
    corticalInputPerOri(:,:,cntr) = a.corticalInput(t);
    lgnInputPerOri(:,:,cntr) = a.lgnOutput(t);
    [~,ix] = max(a.rate);
    popPeak(:,cntr) = a.preferredOrientation(ix);
    rate(:,:,cntr)   = a.rate(:,t);    
end
preferredIx = find(testOri==zeroNeuronPO);


if isempty(testOrientation)
    % Find the test ori
    [v,shiftedPreferredIx] = max(squeeze(mean(rate(zeroNeuronIx,t > stimDuration & t <=2*stimDuration,:),2)));
    shiftedPreferred =testOri(shiftedPreferredIx);
    testOrientation = shiftedPreferred;
    testOrientationIx  = shiftedPreferredIx;
else
    [~,testOrientationIx] =min(abs(testOrientation-testOri));
end


%% Plot the results
figure(3);
clf;
% Set some graphics properties
R=2;C=5;
green = [0.467 0.675 0.188];
blue = [0 0.447 0.741];
red = [0.851 0.325 0.098];
yellow = [1 0.843 0];
linewidth =1.5;
linestyle = ':';
if ~exist('journalFigure','file')
    [fig,sub]= journalFigure('columns',2,'journal','neuron','fontsize',8,'panels',{{R,C,[1 1.7] },{R,C,3:5},{R,C,[6 6.1]},{R,C,[7.2 10]}});
else
    sub(1) = subplot(R,C,[1 1.7]);
    sub(2) = subplot(R,C,3:5);
    sub(3) = subplot(R,C,[6 6.1]);
    sub(4) = subplot(R,C,[7.2 10]);        
end

%% Panel A.  Show the  input to 0-neuron
axes(sub(1));hold off; cla;hold on
h1= plot(testOri,squeeze(corticalInputPerOri(zeroNeuronIx,t==tPoint,:)),'Color',green,'LineWidth',2);
h1b= plot(testOri,squeeze(unadaptedCorticalInputPerOri(zeroNeuronIx,t==tPoint,:)),'Color',green,'LineWidth',2,'linestyle','-.');
h2= plot(testOri,squeeze(lgnInputPerOri(zeroNeuronIx,t==tPoint,:)),'Color',red,'LineWidth',2);

plot([0 0],ylim,'k:')
plot(xlim,[0 0],'k:')
xlabel 'Test Orientation (\circ)'
ylabel 'Input to 0^\circ neuron (mV)'
box on
if limits
    set(gca,'XTick',-90:45:90)
end
legend([h1 h1b, h2 ],'Recurrent Input (adapted)','Recurrent Input (non-adapted)','Feedforward Input','Location','NorthOutside');

%% Panel B -  Show the state of the network at the start of the test stimulus.
axes(sub(2));hold off; cla;hold on
% Plot rate and cortical input for 0 stimulus
[ax,h1,h2] = plotyy(a.preferredOrientation, rate(:,t==tPoint,testOri==0),...
    a.preferredOrientation, corticalInputPerOri(:,t==tPoint,testOri==0)); %#ok<PLOTYY>
hold on
set(h1,'linestyle','-','linewidth',linewidth,'Color',blue);
set(h2,'linestyle','-','linewidth',linewidth,'Color',green);
% Add rate for shifted stimulus
h6 = plot(a.preferredOrientation,rate(:,t==tPoint,testOrientationIx),'LineStyle',linestyle,'LineWidth',linewidth,'Color',blue);
ylabel ('Response (spk/s)')
xlabel 'Preferred Orientation (\circ)'
% Arrow
plot([adapter adapter],[1 0],'LineWidth',2,'LineStyle','-','Color','k');
plot([adapter-cosd(45)*5 adapter],[sind(45)*0.5 0],'LineWidth',linewidth,'LineStyle','- ','Color','k');
plot([adapter+cosd(45)*5 adapter],[sind(45)*0.5 0],'LineWidth',linewidth,'LineStyle','- ','Color','k');
set(ax,'YColor','k')

if limits
    set(ax(1),'YTick',0:4:8,'Ylim',[0 8])
    ax(2).YLim = [-2 2];
    ax(2).YTick = [-2 0 2];
end
ax(2).FontSize = ax(1).FontSize;
set(gca,'XTick',-90:45:90)

axes(ax(2));
% On the right axes, add feedforward drive for the two stimuli
hold on
h3 = plot(ax(2),a.preferredOrientation,lgnInputPerOri(:,t==tPoint,testOri==0));
set(h3,'LineStyle','-','linewidth',linewidth,'Color',red);
h4 = plot(ax(2),a.preferredOrientation,lgnInputPerOri(:,t==tPoint,testOrientationIx));
set(h4,'LineStyle',linestyle,'linewidth',linewidth,'Color',red);
% Add cortical input for the shifted stimulus.
h5 = plot(a.preferredOrientation, corticalInputPerOri(:,t==tPoint,testOrientationIx),'color',green,'linewidth',linewidth,'LineStyle',linestyle);
plot(xlim,[0 0],':k') % 0 input line
ylabel 'Input (mV)'
legend([h1 h6 h3 h4 h2 h5],'Rate \phi=0','Rate \phi=10','FF \phi=0','FF \phi=10','Rec \phi=0','Rec \phi=10','Location','SE');

plot(xlim,[0 0],':k')
plot([0 0],ylim,'k:')
xlabel 'Preferred Orientation (\circ)'
ylabel 'Input (mV)'

%% Panel D. Show the dynamics of rate and input
axes(sub(4));hold off; cla;hold on

[ax,h1,h2]=plotyy(t,squeeze(rate(zeroNeuronIx,:,[preferredIx testOrientationIx])),t,squeeze(corticalInputPerOri(zeroNeuronIx,:,[preferredIx testOrientationIx]))); %#ok<PLOTYY>
set(h1,'LineStyle','-','linewidth',linewidth,'color',blue);
set(h1(2),'LineStyle',linestyle),
set(h2,'LineStyle','-','linewidth',linewidth,'Color',green);
set(h2(2),'LineStyle',linestyle),

xlabel 'Time (ms)'
ylabel ('Response (spk/s)')


axes(ax(2));
hold on
plot(xlim,[0 0],linestyle,'Color',green)
h = plot(t',squeeze(lgnInputPerOri(zeroNeuronIx,:,[preferredIx testOrientationIx])));
set(h,'LineWidth',linewidth,'Color',red)
set(h(2),'LineStyle',linestyle)
ylabel 'Input (mV)'
if limits
    %     set(ax(2),'YLim',[-0.5 1.7],'YTick',-0.5:0.5:1.7,'XLim',[40 simulationDuration]);
    %     set(ax(1),'YTick',2:2:8,'YLim',[2 8],'XLim',[40 simulationDuration],'XTick',(0:3).*stimDuration);
    set(ax(2),'YLim',[-1.5 1.7],'YTick',-1.5:1:1.7,'XLim',[40 simulationDuration]);
    set(ax(1),'YTick',2:2:8,'YLim',[0 8],'XLim',[40 simulationDuration],'XTick',(0:3).*stimDuration);
    % Get rid of right side tick marks
    set(ax(1),'Box','off');
    plot(ax(1),xlim,[8 8],'Color',0.15*[1 1 1 ],'lineWidth',0.5)
end
ax(2).FontSize = ax(1).FontSize;
%legend([h1 h2(1) h(1)],'Firing Rate','Recurrent Input','Feedforward Input','Location','SouthEast');

legend([h1' h'  h2' ],'Rate \phi=0','Rate \phi=10','FF \phi=0','FF \phi=10','Rec \phi=0','Rec \phi=10','Location','Best');

%% Panel C. Show the population peak location dynamics
axes(sub(3));hold off; cla;hold on
blueish = [0.306  0.396 0.58];
blueish2 = [0.729 0.831 0.957];

y = popPeak(:,[preferredIx testOrientationIx]);
plot(xlim,[0 0],'k:');
ix =find(sum(y,2)==0,1,'first');
plot([t(ix) t(ix)],[y(ix,1) y(ix,2)],'LineWidth',linewidth,'Color','k')
patch([t(1:ix) t(ix:-1:1)]',[y(1:ix,1); y(ix:-1:1,2)],'g','FaceAlpha',0.3,'EdgeColor','w');

patch([t(ix+1:end) t(end:-1:ix+1)]',[y(ix+1:end,1); y(end:-1:ix+1,2)],'r','FaceAlpha',0.3,'EdgeColor','w');
h1= plot(t,y,'Color','k','Linewidth',2);
set(h1(2),'LineStyle',linestyle);

set(gca,'Xlim',[40 simulationDuration],'XTick',(0:3).*stimDuration)
ylabel ' Population Peak (\circ)'
xlabel 'Time (ms)'
box on
legend(h1,'\phi=0','\phi=10','Location','SE')
