%% Adaptation without Plasticity - Parameter Exploration
%
%  Code to generate figure 3 in
%  Quiroga M del M, Morris AP, Krekelberg B. Adaptation without plasticity. Cell Rep. 2016;17:58–68.
%  Available from: http://dx.doi.org/10.1016/j.celrep.2016.08.089
%
% BK -May 2016

baseModel  = 'MQ-C';
net = ring(baseModel);
% Common parameters
pv = {'window',[0 20]','testDuration',20,'adapterDuration',20,'plot',false,'fit',true,'adapterOrientation',-82.5:15:0,'testOrientation',-82.5:15:82.5,'subtractMean',true};
net.shrinkFactor = [1 1 ];
[~,~,shift,r2] = tuningCurve(net,pv{:});
baseShift = max(shift);
%% Cortical strength
jCortex = 0:0.25:3;
jCortexShift = nan(numel(jCortex),1);
for i=1:numel(jCortex)
    net = ring(baseModel);
    net.cortexStrength = jCortex(i);
    [~,~,shift,r2] = tuningCurve(net,pv{:});
     stay = r2>0.9;
    if all(stay)
        jCortexShift(i) = max(abs(shift(stay)));
    end   
end

%% I/E ratio
ieRatio = 0.8:0.1:1.8;
ieRatioShift = nan(numel(ieRatio),1);
for i=1:numel(ieRatio)
    net = ring(baseModel);
    net.ieRatio = ieRatio(i);
    [~,~,shift,r2] = tuningCurve(net,pv{:});
     stay = r2>0.9;
    if all(stay)
        ieRatioShift(i) = max(abs(shift(stay)));
    end   
end

%% E tuning;
shrink = 0.8:0.05:1.2;
eWidthShift = nan(numel(shrink),1);
for i=1:numel(shrink)
    net = ring(baseModel);
    net.shrinkFactor(1) = shrink(i);
    [~,~,shift,r2] = tuningCurve(net,pv{:});
    stay = r2>0.9;
    if all(stay)
        eWidthShift(i) = max(abs(shift(stay)));
    end  
end

%% I tuning
iWidthShift = nan(numel(shrink),1);
for i=1:numel(shrink)
    net = ring(baseModel);  
    net.shrinkFactor(2) = shrink(i);
    [~,~,shift,r2] = tuningCurve(net,pv{:});
  
    stay = r2>0.9;
    if all(stay)
       iWidthShift(i) = max(abs(shift(stay)));
    end 
end


%% Show results in a simple figure
figure(1);clf
R=2;C=2;
if exist('journalFigure','file')
    [fig,sub]= journalFigure('columns',1,'journal','neuron','fontsize',8,'panels',{{R,C,1},{R,C,2},{R,C,3},{R,C,4}},'aspectRatio',1);
    journalFigure('figure',fig,'operation','pinch','opArgument',[0.01 0.01]);
else
    fig=gcf;
    sub(1) = subplot(R,C,1);
    sub(2) = subplot(R,C,2);
    sub(3) = subplot(R,C,3);
    sub(4) = subplot(R,C,4);
end

yl = [0 6];
yt =  0:2:6;

net = ring(baseModel);

axes(sub(1));hold on
plot(jCortex,jCortexShift,'k.')
plot(net.cortexStrength,baseShift,'kx');
xlabel ('J_{cortex}')
ylabel 'Shift (\circ)'


axes(sub(2));hold on
plot(ieRatio,ieRatioShift,'k.')
plot(net.ieRatio,baseShift,'kx');
xlabel ('r_{IE}')

axes(sub(3));hold on
plot(1./shrink,eWidthShift,'k.');
plot(1/net.shrinkFactor(1),baseShift,'kx');
xlabel ('Exc. Stretch Factor')
ylabel 'Shift (\circ)'
xlim([min(shrink) max(shrink)]);

axes(sub(4));hold on
plot(1./shrink,iWidthShift,'k.');
plot(1/net.shrinkFactor(2),baseShift,'kx');
xlim([min(shrink) max(shrink)]);
xlabel ('Inh. Stretch Factor')


set(sub,'YLim',yl,'YTick',yt);

