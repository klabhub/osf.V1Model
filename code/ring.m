% Class to investigate the dynamics of the ring model for V1 orientation
% tuning.
%
% MQ, BK, 2016.
classdef ring
    %% Properties
    properties (Access=public)
        %Network properties
        nonlinearity char   = 'piecewise';  % 'piecewise','halfsquare'
        alpha double        = 1;            % Conversion factor: Spikes/mV
        tau double          = 15;           % membrane time constant (ms)
        lgnStrength double  = 10;           % Strength of the LGN input
        cortexStrength double =1.5;         % Strength of the cortical connections
        
        % Threshold adaptation
        thresholdFactor double     = 0;     % Strenght of threshold adaptation . 0 means no threshold adaptation.
        tauThreshold double        = 100;   % Time course of threshold adaptation (> 0 even for adaptFactor =0)
        
        % Divisive normalization  r= r/(sigma+gainFactor*<r>)  where gain linearly
        % integrates the rate of "nearby" neurons (within normWidth).
        gainFactor double           = 0;    % 0 means no divisive normalization
        sigma double                = 1;    % sigma=1, gainFactor =0 means no normalization
        tauGain double              = 100;  % time course of normalization (>0 even for gainfactor =0)
        normWidth double                   = 2.5;  % Normalization width (K)
        
        % Lateral connectivity
        lgnWidth double          = 0.7;     % Tuning width of the LGN input (von Mises K)
        excitatoryWidth double   = 1.8;     % Width of the excitatory intracortical connections (K)
        inhibitoryWidth double   = 1.5;     % Width of the inhibitory intracortical connections (K)
        inhibitoryCenter double  = 0;       % Offset of the peak of the inhibitory connections (deg)
        ieRatio double           = 1;       % Ratio of Inhibition to Excitation. For two population model: 1 = inhibition of excitatory, 2 = inhibitoion of inhibitory
        shrinkFactor double      =[1 1];    % Becuase the von Mises is scaled to the sum, changing the width also changes the peak.
        % To change only the width, use a shrinkFactor [Exc Inh]
        
        % Stimulus properties
        % For simple gratings specify stimulus orientation and contrast in
        % the call to simulate, for more complex sequences, define a
        % function and assign to stimulus sequence. The function will be called with
        % a time and should return the activity of each input neuron
        % [nrNeurons 1]
        stimulusSequence = [];              % A function returning orientation and activity.
        
        % Model "attention"
        attentionSequence = @(t)(0);         %  A function returning attentional modulation; the additional input to the cortical neurons at time t.
        
        % Transmission delays are 0 by default
        delayRetina=0;        % Delays in ms - time until LGN gets its first input.
        delayAfferent = 0;        % Delay of transmission along the axon from LGN to V1
        delayLateral =0;        % Delay along lateral, intra-V1 connections
        
        % Noise is 0 by default
        noiseVoltage =0;       % Standard deviation of the Gaussian noise added to the voltage each time step.
        noiseRate = 0;          % Standard deviation of the Gaussian noise added to the rate each time step.
    end
    
    properties (GetAccess =public, SetAccess=protected)
        time;           % Simulation time steps
        rate;           % Firing rate for each neuron and time step [nrNeurons nrTimeSteps]
        voltage;        % Membrane potential [nrNeurons nrTimeSteps]
        adaptation;     % Adaptation state [nrNeurons nrTimeSteps]
        gain;           % Normalization state [nrNeurons nrTimeSteps]
        
        % Each networks is derived from some base model:
        baseModel char;                     % MQ-C, MQ-M, MQ-SLOW, TEICH, CARANDINI
        
        preferredOrientation;  % The preferred orientation of the neurons in the network
        
        excitatoryWeights;      % The excitatory connection strengths (translation invariant)
        inhibitoryWeights;      % The excitatory connection strengths
        cortexWeights;          % The combined E+I strengths
        normalizationWeights;   % The weights used to calculate normalization factos
        connectionsScale;       % Used for non-translation invariant connections.
    end % Properties
    
    properties (GetAccess= protected, SetAccess=protected)
        convolutionWeights;         % The circshifted weights for the convolution - for translation invariant connections
        matrixWeights;              % Used for non translation invariant connections. See strengthen. [nrNeurons nrNeurons]. Row 1 is the connectivity going into neuron 1.
        preComputedSequence;        % A sequence of activities of the input units [nrNeurons nrTimeSteps]
    end
    
    properties (Dependent)
        nrNeurons;                  % Total number of neurons in the cortical layer
        nrOrientations;             % Same as the number of neurons...
        isTwoPopulations;           % Boolean indicating whether the model contains separate Excitatory and Inhibitory populations
    end
    
    methods %get/set
        function v =get.nrNeurons(o)
            v = numel(o.preferredOrientation);
        end
        
        function v =get.nrOrientations(o)
            v = numel(o.preferredOrientation);
        end
        
        function o=set.nrNeurons(o,v)
            o.preferredOrientation = linspace(-90,90,v)';
        end
        function v=get.isTwoPopulations(o)
            v = numel(o.excitatoryWidth)>1;
        end
        
        function o = set.stimulusSequence(o,v)
            % Set the sequence of stimuli that will be presented.
            % v = usually a function that takes a time as input and
            % returns the activity at each orientation. [nrNeurons 1]
            % v can also be a [2 T] vector with T the number of
            % milliseconds since time 0. The first row represents the
            % orientation, the second the activity/contrast presented at
            % that time.
            %  [0 90;1 0.5] -> a horizontal grating with contrast 1 at t=0
            %  and a vertical grating with contrast 0.5 at t=1 ms.
            if isa(v,'function_handle')
                o.stimulusSequence =v;
            elseif isnumeric(v)
                o.preComputedSequence = v; %#ok<MCSUP> % Store the values
                % Define a function that will return the approriate input
                % values for each time.
                o.stimulusSequence = @(t) (inputFromPrecompute(o,t));
            end
        end
        
    end
    
    %% Methods
    methods
        %% Constructor function.
        function o = ring(bm,twoPopulations)
            % Specify the base model and whether to use separate Excitatory
            % and Inhibtory populations
            if nargin<2
                twoPopulations = false;
            end
            o.baseModel =bm;
            o= setBaseParameters(o,twoPopulations);
        end
        
        function o=strengthen(o,varargin)
            % Strengthen the connections to a subset of neurons
            p = inputParser;
            p.addParameter('preferredOrientation',0);
            p.addParameter('preferredOrientationWidth',10);
            p.addParameter('scale',1);
            p.parse(varargin{:});
            
            if o.isTwoPopulations
                error('Not implemented yet');
            end
            
            if isscalar(p.Results.scale)
                o.connectionsScale = ones(o.nrNeurons,1);% Start with same weights for all.
                stay = abs(o.preferredOrientation-p.Results.preferredOrientation)<p.Results.preferredOrientationWidth/2;
                o.connectionsScale(stay,:) = p.Results.scale; % Within width/2 from po, all get teh same scale.
            elseif numel(p.Results.scale)==o.nrNeurons
                o.connectionsScale = p.Results.scale(:);% Force col.
            else
                error('The ''scale'' does not match the number of neurons');
            end
            
        end
        function input = attendTo(o,ori,strength,width,disinhibition)
            % Model attention by boosting the activity of a subset of
            % neurons (or reducing their inhibition).
            if nargin <4
                disinhibition = false;
            end
            input = zeros(o.nrNeurons,1);
            input(abs(o.preferredOrientation-ori)/90<0.5*width)=strength;
            if o.isTwoPopulations
                if disinhibition
                    % reduced the inhibitory neurons
                    input = [zeros(o.nrNeurons,1) -abs(input)];
                else
                    % boost to the excitatory neurons
                    input = [abs(input) zeros(o.nrNeurons,1)];
                end
            end
        end
        
        function v = lgnOutput(o,t,po)
            % Returns the output of the LGN layer at time t for all
            % neurons, or for the neuron with preferred orientation po.
            if nargin<3
                nIx = 1:o.nrNeurons;
            else
                nIx = find(ismember(o.preferredOrientation,po));
            end
            v =(nan(numel(nIx),1+o.isTwoPopulations,numel(t)));
            
            for i=1:numel(t)
                % Becuase the model assumes that the LGN processes
                % instantenously, we can combine the retinal delay and the
                % delay between LGN and V1. I.e. this is a pure delay line
                stimulus  = o.stimulusSequence(t(i)-o.delayRetina-o.delayAfferent);
                
                switch upper(o.baseModel)
                    case {'TEICH','CARANDINI'}
                        lgnFilter =exp(-(o.preferredOrientation./(2*o.lgnWidth)).^2);
                    case 'COS2'
                        %TODO: check   o.lgnOutput = o.lgnStrength.*(1-o.lgnWidth +o.lgnWidth.*(cosd(2*(ori-o.preferredOrientation))));
                        o.lgnOutput = repmat(o.lgnOutput,[1 2]);
                        return;
                    case {'MQ-C','MQ-M','MQ-SLOW'}
                        lgnFilter = ring.vonMises(o.preferredOrientation,0,o.lgnWidth,180);
                    otherwise
                        error(['Unknown base model: ' o.baseModel]);
                end
                % Convolve stimulus with LGN filter
                lgnFilter = circshift(lgnFilter,round(o.nrNeurons/2));
                tmp = ifft2(fft2(stimulus).*fft2(lgnFilter));
                
                if o.isTwoPopulations
                    v(:,:,i) = o.lgnStrength.*repmat(tmp(nIx,:),[1 2]);
                else
                    v(:,:,i) = o.lgnStrength.*tmp(nIx);
                end
            end
            v=squeeze(v);
        end
        
        function o= setBaseParameters(o,twoPopulations)
            % Setup the base parameters based on the choice of baseModel
            switch upper(o.baseModel)
                case 'MQ-C'
                    % The "cat" model; best fit to the Felsen et al. data.
                    o.preferredOrientation = linspace(-90,90-180/256,256)';
                    o.alpha             = 10.606627806236400;
                    o.tau               = 10.762315360263232;
                    o.lgnStrength       = 9.569804305270075;
                    o.lgnWidth          = 1.560433795865845;
                    o.cortexStrength    = 1.706513465281997;
                    o.ieRatio           = 1.178813258661855;
                    o.excitatoryWidth   = 1.586832104297276;
                    o.inhibitoryWidth   = 1.158469310525126;
                    o.inhibitoryCenter  = 0;
                    
                case 'MQ-M'
                    % The monkey model; best fit to the Patterson et
                    % al.data.
                    o.preferredOrientation = linspace(-90,90-180/256,256)';
                    o.alpha             = 3.882189013814953;
                    o.tau               = 8;
                    o.lgnStrength       = 11.041389802178394;
                    o.lgnWidth          = 0.473559847094274;
                    o.cortexStrength    = 2.835352731049699;
                    o.ieRatio           = 1.242695980763933;
                    o.excitatoryWidth   = 1.118193314120349;
                    o.inhibitoryWidth   = 0.561309663822524;
                    o.inhibitoryCenter = 0;
                case 'MQ-SLOW'
                    % A model with slow dynamics
                    o.preferredOrientation = linspace(-90,90-180/256,256)';
                    o.alpha             = 4;
                    o.tau               = 15;
                    o.lgnStrength       = 8;
                    o.lgnWidth          = 0.5;
                    o.cortexStrength    = 1.7;
                    o.ieRatio           = 1.14;
                    o.excitatoryWidth   = 2.2;
                    o.inhibitoryWidth   = 1.0;
                    o.inhibitoryCenter  = 0;
                case 'TEICH'
                    % Implementation of the Teich & Qian model
                    o.preferredOrientation = linspace(-90,90-180/128,128)';
                    o.alpha = 10;
                    o.lgnStrength = 10;
                    o.lgnWidth = 45;
                    o.cortexStrength = 1.1;
                    o.ieRatio = 1.2;
                    o.excitatoryWidth = 2.2;
                    o.inhibitoryWidth = 1.4;
                    o.inhibitoryCenter = 0;
                    
                case 'CARANDINI'
                    % Implementation of the Carandini model
                    o.preferredOrientation = linspace(-90,90-180/512,512)';
                    o.alpha = 15;
                    o.lgnStrength = 3.2;
                    o.lgnWidth = 23;
                    o.cortexStrength = 0.115;
                    o.ieRatio = 50/23;
                    o.excitatoryWidth = 7.5;
                    o.inhibitoryWidth = 60;
                    o.inhibitoryCenter = 0;
                    
                case 'COS2'
                    o.preferredOrientation = linspace(-90,90-180/256,256)';
                    o.alpha = 10;
                    o.lgnStrength = 10;
                    o.lgnWidth = 0.1; % epsion
                    o.cortexStrength = 1;
                    o.ieRatio = [1 1];
                    o.excitatoryWidth = [1.5 1.5];
                    o.inhibitoryWidth = [1.8 1.8];
                    o.inhibitoryCenter = 0;
                    
                otherwise
                    error(['Unknown base model: ' o.baseModel]);
            end
            
            if twoPopulations
                % By default we simply duplicate the parameters for each of
                % the two populations.
                o.excitatoryWidth = [o.excitatoryWidth o.excitatoryWidth];
                o.inhibitoryWidth = [o.inhibitoryWidth o.inhibitoryWidth];
                o.ieRatio         = [o.ieRatio o.ieRatio];
                o.lgnStrength     = [o.lgnStrength o.lgnStrength];
            end
            
        end
        
        
        
        function o= initializeWeights(o)
            % Setuo the weights based on the parameter settings.
            switch upper(o.baseModel)
                case {'MQ-C','MQ-M','MQ-SLOW'}
                    if o.isTwoPopulations
                        ee = ring.vonMises(o.preferredOrientation,0,o.excitatoryWidth(1),180); % from e to e
                        ei = ring.vonMises(o.preferredOrientation,o.inhibitoryCenter,o.excitatoryWidth(2),180);
                        ie = ring.vonMises(o.preferredOrientation,0,o.inhibitoryWidth(1),180); % from i to e
                        ii = ring.vonMises(o.preferredOrientation,o.inhibitoryCenter,o.inhibitoryWidth(2),180);
                        ee = ee./sum(ee);ei =ei/sum(ei);ii=ii/sum(ii);ie = ie/sum(ie);
                        o.excitatoryWeights = [ee ei];
                        o.inhibitoryWeights = [ie ii];
                        o.cortexWeights = o.cortexStrength.*[ee -o.ieRatio(1).*ie ei  -o.ieRatio(2).*ii]; % ee, ie, ei, ii;
                    else
                        o.excitatoryWeights = ring.vonMises(o.preferredOrientation,0,o.excitatoryWidth,180);
                        o.inhibitoryWeights = ring.vonMises(o.preferredOrientation,o.inhibitoryCenter,o.inhibitoryWidth,180);
                    end
                case 'TEICH'
                    o.excitatoryWeights = (cosd(2*o.preferredOrientation)+1).^o.excitatoryWidth;
                    o.inhibitoryWeights = (cosd(2*o.preferredOrientation)+1).^o.inhibitoryWidth;
                case 'CARANDINI'
                    o.excitatoryWeights = exp(-o.preferredOrientation.^2./(2*o.excitatoryWidth^2));
                    o.inhibitoryWeights = exp(-o.preferredOrientation.^2./(2*o.inhibitoryWidth^2));
                    o.inhibitoryWeights(abs(o.preferredOrientation)>o.inhibitoryWidth) = 0;
                    
                case 'COS2'
                    ee = ring.warpCosd(2.*o.preferredOrientation,0); % from e to e
                    ei = ring.warpCosd(2.*o.preferredOrientation,0);
                    ii = ring.warpCosd(2.*o.preferredOrientation,0);
                    ie = ring.warpCosd(2.*o.preferredOrientation,0);
                    ee = ee./sum(ee);ei =ei/sum(ei);ii=ii/sum(ii);ie = ie/sum(ie);
                    o.excitatoryWeights = [ee ei];
                    o.inhibitoryWeights = [ie ii];
                    o.cortexWeights = o.cortexStrength.*[ee -o.ieRatio(1).*ie ei  -o.ieRatio(2).*ii]; % ee, ei, ie, ii;
                    
                otherwise
                    error(['Unknown base model: ' o.baseModel]);
            end
            
            
            
            if ~o.isTwoPopulations
                % Normalization
                o.excitatoryWeights = o.excitatoryWeights./sum(o.excitatoryWeights); % Each weight vector is normalized to a sum of 1
                o.inhibitoryWeights = o.inhibitoryWeights./sum(o.inhibitoryWeights);
                o.cortexWeights = o.cortexStrength.*(o.excitatoryWeights-o.ieRatio.*o.inhibitoryWeights);
            end
            
            if any(o.shrinkFactor ~=1)
                o = shrink(o); % recalculare weights using the current shrnk
            end
            
            o.normalizationWeights = circshift(ring.vonMises(o.preferredOrientation,0,o.normWidth,180),round(o.nrNeurons/2));
            if isempty(o.connectionsScale)
                % Translation invariant connectivity, use FFT to do
                % integration in dynamics.
                o.convolutionWeights = circshift(o.cortexWeights ,round(o.nrNeurons/2));
            else
                o.matrixWeights = nan(o.nrNeurons,o.nrNeurons);
                o.matrixWeights(1,:) = circshift(o.cortexWeights ,round(o.nrNeurons/2));
                for i=2:o.nrNeurons
                    o.matrixWeights(i,:) = circshift(o.matrixWeights(i-1,:),1);
                end
                o.matrixWeights =o.matrixWeights.*repmat(o.connectionsScale,[1 o.nrNeurons]);
            end
            
        end
        
        
        function [oriPerT,r] =decode(o,varargin)
            % Simple decoder mapping population activity to "perceived"
            % orientation based on the vector sum, or peak.
            % This decoded can also include some suppression (e.g. due to
            % threshold adaptation).
            p = inputParser;
            p.addParameter('window',[-inf inf]);
            p.addParameter('adaptWindow',[0 200]);
            p.addParameter('model','vectorsum');
            p.addParameter('suppression',0);
            p.addParameter('tauRecovery',200);
            p.parse(varargin{:});
            
            % Implement a simple form of plasticity: take the mean rate in
            % the adaptation window, and subtract a fraction of this from
            % the rate. Let this effect decay exponentially  with
            % tauRecovery.
            if p.Results.suppression >0
                adStay =  o.time >= min(p.Results.adaptWindow) &  o.time <=max(p.Results.adaptWindow);
                adapt = repmat(mean(o.rate(:,adStay),2),[1 size(o.rate,2)]);
                adapt = adapt.*repmat(min(1,exp(-(o.time'-max(p.Results.adaptWindow))/p.Results.tauRecovery)),[o.nrNeurons 1]);
                r = max(o.rate - p.Results.suppression*adapt,0);
            else
                % No rate suppression due to adaptation.
                r = o.rate;
            end
            
            stay  = o.time >= min(p.Results.window) &  o.time <=max(p.Results.window);
            switch upper(p.Results.model)
                case 'VECTORSUM'
                    tmp  = 0.5*circ_mean(2*repmat(o.preferredOrientation*pi/180,[1 sum(stay)]),r(:,stay),1);
                    oriPerT = tmp*180/pi;
                case 'PEAK'
                    [~,ix] = max(r(:,stay));
                    oriPerT= o.preferredOrientation(ix)';
            end
            
            
        end
        
        function [r,contrast]= crf(o,varargin)
            p = inputParser;
            p.addParameter('contrast',0:5:100,@isnumeric);
            p.addParameter('duration',250,@isnumeric);
            p.addParameter('orientation',0,@isnumeric);
            p.parse(varargin{:});
            
            r = nan(size(p.Results.contrast));
            ix = o.preferredOrientation == 0; % Always report PO=0 neuron
            for i=1:numel(p.Results.contrast)
                o = simulate(o,'tSpan',[1 p.Results.duration],'t',1:p.Results.duration,'continue',false,'stimulusContrast',p.Results.contrast(i),'stimulusOrientation',p.Results.orientation);
                r(i)= nanmean(o.rate(ix,:));
            end
            if nargout >1
                contrast =p.Results.contrast;
            end
        end
        
        %% Solves model numerically
        function [o,rate] = simulate(o,varargin)
            p = inputParser;
            p.addParameter('tSpan',[0 250]);
            p.addParameter('t',[]);
            p.addParameter('continue',false);
            p.addParameter('stimulusOrientation',[],@isnumeric);
            p.addParameter('stimulusContrast',1,@isnumeric);
            p.addParameter('options',struct,@isstruct); % ddeset or odeset options
            p.parse(varargin{:});
            [vIx,aIx,gIx] = o.vag;
            if any(p.Results.t > max(p.Results.tSpan))
                tSpan = [p.Results.tSpan max(p.Results.t)];
            else
                tSpan = p.Results.tSpan;
            end
            
            if ~isempty(p.Results.stimulusOrientation)
                o.stimulusSequence = @(~) singleGratingToInput(o,p.Results.stimulusOrientation,p.Results.stimulusContrast);
                %else  simulate the current value/sequence
            end
            o  = initializeWeights(o);
            
            if p.Results.continue && ~isempty(o.voltage)
                initialState = [o.voltage(:,end,:) o.adaptation(:,end,:) o.gain(:,end,:)];
            else
                initialState = zeros(o.nrNeurons,(o.isTwoPopulations+1)*3); %[Voltage Adaptation Gain]
            end
            
            defaultOpts = struct('RelTol',1e-3, 'NormControl','on');
            options = odeset(defaultOpts,p.Results.options);
            if o.delayLateral==0
                % No delays, regular ode
                sol = ode45(@dynamics,tSpan,initialState,options);
            else
                % use delay differential equation solver
                history = reshape(initialState,o.nrNeurons,numel(initialState)/o.nrNeurons);
                sol = dde23(@dynamics,o.delayLateral,history, tSpan,options);
            end
            
            if isempty(p.Results.t)
                t = p.Results.tSpan;
            else
                t = p.Results.t;
            end
            
            state = deval(sol,t);
            
            
            % Reshape
            state = reshape(state,o.nrNeurons,(o.isTwoPopulations+1)*3,length(t));
            state = permute(state,[1 3 2]); % [Neurons, Time, Ve Vi Ae Ai Ge Gi]
            
            r = stateToRate(o,state);
            
            
            if p.Results.continue
                o.time = cat(1,o.time,t(:));
                o.rate = cat(2,o.rate,r(:,:));
                o.voltage = cat(2,o.voltage,state(:,:,vIx));
                o.adaptation   = cat(2,o.adaptation,state(:,:,aIx));
                o.gain  = cat(2,o.gain,state(:,:,gIx));
            else
                o.time = t(:);
                o.rate = r;
                o.voltage = state(:,:,vIx);
                o.adaptation   = state(:,:,aIx);
                o.gain = state(:,:,gIx);
            end
            
            
            
            function r= stateToRate(o,state)
                %% function r= stateToRate(o,state)
                % Convert a state (membrange voltage, threshold, gain) of
                % the network to the firing rate
                if ndims(state)==3
                    volt = state(:,:,vIx);
                    thresh= state(:,:,aIx);
                    gn  = state(:,:,gIx);
                else
                    volt = state(:,vIx);
                    thresh= state(:,aIx);
                    gn  = state(:,gIx);
                end
                
                volt = volt + o.noiseVoltage*randn;
                % voltage to rate conversion
                switch o.nonlinearity
                    case 'piecewise'
                        r = o.alpha*max(volt-thresh,0);
                    case 'halfsquare'
                        r = o.alpha*max(volt-thresh,0).^2;
                    case 'heaviside'
                        r = o.alpha*(volt>thresh);
                    otherwise
                        error(['Unknown nonlinearity: ' o.nonlinearity]);
                end
                
                % Simulate spontaneous rate by adding spiking noise
                r = r + o.noiseRate*randn;
                
                % Divisive normalization
                r = r./( o.sigma + gn);
                
            end
            
            
            %% Model: nested function with access to the ring object
            function dstatedt = dynamics(t,state,delayedState)
                if nargin<3
                    % Ignore delays.
                    delayedState = state;
                end
                
                
                [vIx,aIx,gIx] = vag(o);
                state = reshape(state,o.nrNeurons,numel(state)/o.nrNeurons);
                r= stateToRate(o,state);
                
                lgn = o.lgnOutput(t);
                
                % Lateral interactions
                delayedState = reshape(delayedState,o.nrNeurons,numel(delayedState)/o.nrNeurons);
                delayedRate = stateToRate(o,delayedState); %EI  : with repmat the integral is over EIEI
                if isempty(o.matrixWeights)
                    % Use convolution weights. For 2-pop model: [ee ie ei ii]
                    Vcortex = ifft(fft(o.convolutionWeights).*repmat(fft(delayedRate),[1 1+o.isTwoPopulations]));
                else
                    % Use the full matrix
                    Vcortex = o.matrixWeights*repmat(delayedRate,[1 1+o.isTwoPopulations]);
                end
                if o.isTwoPopulations
                    % Two population model; sum over input to e (1:2) and
                    % sum over input to i (3:4)
                    Vcortex = [sum(Vcortex(:,1:2),2) sum(Vcortex(:,3:4),2)]; %e i
                end
                
                % Attention - modeled as an additional input to a set of
                % neurons.
                att = o.attentionSequence(t);
                
                
                % Compute the derivative
                dvdt = ( -state(:,vIx) + lgn + Vcortex +att )./ o.tau;
                dvdt = dvdt(:);
                
                %Threshold dynamics
                dadt = (-state(:,aIx) +max(0,o.thresholdFactor.*r))/o.tauThreshold;
                
                % Population normalization/ gain control
                normalizationSignal = ifft(fft(o.normalizationWeights).*fft(delayedRate));
                dgdt = (-state(:,gIx) + o.gainFactor*normalizationSignal )/o.tauGain;
                
                dstatedt = [dvdt; dadt(:);dgdt(:)];
                
                
            end % dynamics
            if nargout >1
                rate = o.rate;
            end
        end
        
        function ctx = corticalInput(o,t,po)
            if nargin<3
                po = o.preferredOrientation;
            end
            ctx = nan(o.nrNeurons,numel(t));
            for i=1:numel(t)
                tIx = find(o.time==t(i));
                if isempty(tIx)
                    error(['Simulate with t including ' num2str(t)  ' first']);
                end
                nIx = find(ismember(o.preferredOrientation,po));
                if isempty(nIx)
                    error(['No neuron with PO = ' num2str(po)]);
                end
                if o.isTwoPopulations
                    % Two population model
                    error('todo')
                else
                    ctx(:,i) = ifft(fft(o.convolutionWeights).*fft(o.rate(:,tIx,:)));
                end
            end
            ctx = ctx(nIx,:); % Total cortical input to po at time t
        end
        
        
        %% Generate tuning curve
        function  [rate,parms,shift,r2]= tuningCurve(o,varargin)
            p=inputParser;
            p.addParameter('testOrientation',18,@isnumeric);
            p.addParameter('testDuration',max(o.time),@isnumeric);
            p.addParameter('testContrast',1,@isnumeric); % The contrast of the testers
            
            p.addParameter('adapterDuration',0,@isnumeric);
            p.addParameter('adapterOrientation',[],@isnumeric); % A vector of adapter orientations
            p.addParameter('adapterContrast',1,@isnumeric); % The contrast of the adapters
            
            p.addParameter('blankDuration',0,@isnumeric);
            p.addParameter('window',[50 100]',@isnumeric);
            p.addParameter('preferredOrientation',0,@isnumeric);
            p.addParameter('subtractMean',false,@islogical);
            p.addParameter('subtractAdapter',false,@islogical);
            p.addParameter('plot',false,@(x) (islogical(x) || ishandle(x)));
            p.addParameter('fit',false,@islogical);
            
            p.parse(varargin{:});
            
            
            [~,neuronIX] = min(abs(o.preferredOrientation- p.Results.preferredOrientation));
            
            if isscalar(p.Results.testOrientation)
                testOrientations = linspace(-90,90-180/p.Results.testOrientation,p.Results.testOrientation);
            else
                testOrientations =p.Results.testOrientation;
            end
            nrOris= numel(testOrientations);
            nrWindows = size(p.Results.window,2);
            nrAdapters = max(1,numel(p.Results.adapterOrientation));
            
            if isempty(p.Results.adapterOrientation)
                adapterDuration = 0;
            else
                adapterDuration = p.Results.adapterDuration;
            end
            
            rate = nan(nrOris,1+o.isTwoPopulations,nrAdapters,nrWindows);
            adapterResponse = nan(1,1+o.isTwoPopulations,nrAdapters);
            shift = nan(nrAdapters,1);
            parms = nan(nrAdapters,4);
            
            for a =1:nrAdapters
                if adapterDuration>0 && ~isempty(p.Results.adapterOrientation)
                    % Adapt from scratch
                    o = simulate(o,'tSpan',[0 adapterDuration],'continue',false,'stimulusContrast',p.Results.adapterContrast, 'stimulusOrientation',p.Results.adapterOrientation(a),'t',1:adapterDuration);
                    for w=1:size(p.Results.window,2)
                        stayT = o.time >= p.Results.window(1,w) & o.time <= p.Results.window(2,w);
                        adapterResponse(1,:,a,w) = mean(o.rate(neuronIX,stayT,:),2);
                    end
                end
                if p.Results.blankDuration>0
                    % Simulate blank
                    o = simulate(o,'tSpan',[adapterDuration+1 adapterDuration+p.Results.blankDuration],'continue',true,'stimulusContrast',0,'stimulusOrientation',0);
                end
                % Now determine tuning
                for i = 1:nrOris
                    oTest = simulate(o,'continue',true,'t',(adapterDuration+p.Results.blankDuration:adapterDuration+p.Results.blankDuration+p.Results.testDuration),'tSpan',[adapterDuration+p.Results.blankDuration adapterDuration+p.Results.blankDuration+p.Results.testDuration],'continue',true,'stimulusOrientation',testOrientations(i),'stimulusContrast',p.Results.testContrast); % Get a rate for every ms.
                    for w=1:size(p.Results.window,2)
                        testTime = oTest.time-adapterDuration-p.Results.blankDuration;
                        stayT =testTime >= p.Results.window(1,w) & testTime <= p.Results.window(2,w);
                        rate(i,:,a,w) = mean(oTest.rate(neuronIX,stayT,:),2);
                    end
                end
                
            end
            
            if p.Results.subtractAdapter
                rate = rate - repmat(adapterResponse,[nrOris 1 1 1]);
            end
            if p.Results.subtractMean
                rate = rate- repmat(mean(rate,1),[nrOris 1 1 1]);
            end
            
            
            r2 = nan(nrAdapters,1+o.isTwoPopulations);
            h = nan(nrAdapters,1);
            for a =1:nrAdapters
                for pop =1:1+o.isTwoPopulations
                    for w=1:size(p.Results.window,2)
                        [peak,ix] = max(rate(:,pop,a,w),[],1);
                        po = testOrientations(ix);
                        if p.Results.fit
                            fun = @(parms,x) (parms(4) + parms(3).*ring.vonMises(x,parms(1),parms(2),180));
                            guess = [po 1 2*peak 0];
                            options = optimset('Display','off');
                            parms(a,:,pop,w) = lsqcurvefit(fun,guess,testOrientations',rate(:,pop,a,w),[],[],options);
                            po = parms(a,1,pop,w);
                            cc = corrcoef(fun(parms(a,:,pop,w),testOrientations),rate(:,pop,a,w));
                            r2(a,pop,w) = cc(1,2).^2;
                        end
                        shift(a,pop,w) = po-p.Results.preferredOrientation;
                        styles = '+o'; % populations
                        if p.Results.plot ~=0
                            if ishandle(p.Results.plot)
                                h(a) = p.Results.plot(a);
                            else
                                if a==1 && pop==1 && w==1; clf;end
                                R = ceil(sqrt(nrAdapters));
                                C= max(2,R);
                                h(a) = subplot(R,C,a); 
                            end
                            hold on
                            plot(testOrientations,rate(:,:,a,w),styles(pop))
                            
                            plot(p.Results.preferredOrientation*[1 1],ylim,'k-'); % Label (Preferred orientation) of the neuron
                            
                            
                            xlabel 'Orientation (\circ)'
                            ylabel 'Response (Hz)'
                            if ~isempty(p.Results.adapter)
                                plot(p.Results.adapter(a)*[1 1],ylim,'b-'); % Adapter
                                title (['Tuning Curve. Neuron: ' num2str(p.Results.preferredOrientation) '^\circ Adapter: ' num2str(p.Results.adapter(a)) '^\circ']);
                            else
                                title (['Tuning Curve. Neuron: ' num2str(p.Results.preferredOrientation) '^\circ']);
                            end
                            if p.Results.fit
                                plot(testOrientations,fun(parms(a,:,pop,w),testOrientations),'-'); %Fitted curve
                            end
                            plot([po po],ylim, ':'); % preferred orientation
                            quiver(p.Results.preferredOrientation,0,po-p.Results.preferredOrientation,0,0,'LineWidth',2);
                            xlim([-90 90]);
                        end
                    end
                end
            end
            
            if p.Results.plot~=0 &&  ~isempty(p.Results.adapter)
                if ishandle(p.Results.plot)
                    axis(p.Results.plot(a));
                else
                    subplot(R,C,min(a+1,R*C));
                end
                cla;
                hold on
                for w= 1:size(p.Results.window,2)
                    plot(p.Results.adapter,squeeze(shift(:,:,w)));
                end
                xlabel 'Adapter \circ'
                ylabel 'Shift \circ'
            end
            
            
        end
        
        function [ttRate,ttp,lambda,rates,t] = latency(o,varargin)
            % Present a flash to the network and estimate the latency, using
            % a few commond definitions of latency:
            % ttp = Time to peak
            % ttRate = Time to a fixed rate of response.
            % lambda = Time at which the mean response exceeds the baseline response by a threshold number of standard deviations.
            % rates = rates of the neuron stimulated with its preferred
            % orientation.
            % t = time since stimulus onset.
            %
            % The lambda computation only makes sense when there is noise
            % in the network. The ttRate is a shortcut to compute lambda in
            % the absence of noise. (the fixed rate is then the user's estimate of
            % the response that is needed to get above the noise level)
            %
            p=inputParser;
            p.addParameter('nrTrials',1); % Set this to more than 1 to compute lamda
            p.addParameter('threshold',3); % z-scores above baseline
            p.addParameter('stimDuration',50); % Duration of the flashed stimulus
            p.addParameter('baseline',50); % duration of the pre-stimulus onset baseline
            p.addParameter('contrast',1); % Contrast of the stimulus/
            p.addParameter('plot',false); % Set to true, or pass a handle to an axis.
            p.addParameter('fixedRate',10); % The level to use for the ttRate.
            p.parse(varargin{:});
            
            stimOri = 0;  % Stimulate at this orientation, and use the response of the neuron that prefers this ori.
            AFTERSTIM = 2;  % Time period to simulate after the flash is off
            totalSimulationDuration = p.Results.baseline + p.Results.stimDuration + AFTERSTIM*p.Results.stimDuration;
            contrast = [zeros(1,p.Results.baseline)  p.Results.contrast*ones(1,p.Results.stimDuration) +zeros(1,AFTERSTIM*p.Results.stimDuration) ];  % Contrast zero during baseline then
            ori      = stimOri*ones(1,totalSimulationDuration);
            o.stimulusSequence = [ori;contrast];
            [~,neuronIx] = min(abs(o.preferredOrientation-stimOri));
            t = 0:totalSimulationDuration;
            rates = nan(numel(t),p.Results.nrTrials);
            for i=1:p.Results.nrTrials
                o = simulate(o,'tSpan',[min(t) max(t)],'t',t,'continue',false);
                rates(:,i) = o.rate(neuronIx,:)';
            end
            isBaseline = o.time < p.Results.baseline;
            baselineRate = mean(rates(isBaseline,:),1);
            mSpont       = mean(baselineRate,2); % Average across trials
            sdSpont      = std(baselineRate,0,2); % SD across trials
            z = abs(mean(rates,2)-mSpont)/sdSpont; % Express rates as z of spontaneous
            firstAboveIx = find(z>p.Results.threshold,1,'first'); % Find the first threshold crossing.
            lambda = o.time(firstAboveIx)-p.Results.baseline; % Subtract stimulus onset time
            
            %% Time to half peak
            meanRates = mean(rates,2); % Trial average
            [pkVal,ix] = max(meanRates);
            [~,halfIx] = min(abs(meanRates(1:ix,:)-0.5*pkVal));
            ttp = t(halfIx)-p.Results.baseline;
            %% Time to fixed rate
            [~,levelIx] =  min(abs(meanRates(1:ix,:)-p.Results.fixedRate));
            ttRate = t(levelIx)-p.Results.baseline;
            
            t= t-p.Results.baseline;
            if p.Results.plot ~=0
                if ishandle(p.Results.plot)
                    axes(p.Results.plot);
                end
                plot(t,rates);
                xlabel 'Time since stim onset (ms)'
                ylabel 'Response (spk/s)';
                hold on;
            end
        end
        
        
        
        
        %% Make plots like those in Teich and Qian, 2003
        function plotTeichQian(o,t)
            
            %% Figure 1
            f = findobj('name',['Teich&Qian: ' o.baseModel]);
            if isempty(f)
                f=figure;
                set(f,'name',['Teich&Qian: ' o.baseModel]);
            end
            figure(f);
            clf;
            if nargin<2
                t= max(rnn.time);
            end
            subplot(2,2,1);
            plot(o.preferredOrientation,o.lgnStrength.*o.lgnOutput(t),'linewidth',2)
            xlabel('Stimulus orientation (\circ)')
            ylabel('Synaptic potential (mV)')
            xlim([-90 90]);set(gca,'XTick',-90:45:90);
            
            a= subplot(2,2,2);
            % Careful! Multiplies by N/180...
            a.LineStyleOrder = {'-',':'};
            plot(o.preferredOrientation,o.nrOrientations/180*o.excitatoryWeights,'linewidth',2);
            hold on;
            plot(o.preferredOrientation,o.nrOrientations/180*o.inhibitoryWeights);
            legend('exc','inh')
            xlabel('\Delta Preferred Orientation (?)')
            xlim([-90 90]);set(gca,'XTick',-90:45:90);
            ylabel('Probability of connection')
            
            a=subplot(2,2,3);
            a.LineStyleOrder = {'-',':'};
            plot(o.preferredOrientation,o.cortexWeights,'linewidth',2);
            hold on;
            plot(xlim, [0 0],'k--');
            ylabel('Connection profile');
            xlabel('\Delta Preferred orientation (?)');
            xlim([-90 90]);set(gca,'XTick',-90:45:90,'YTick',0);
            if o.isTwoPopulations
                legend('ee','ei','ie','ii');
            end
            
            a= subplot(2,2,4);
            a.LineStyleOrder = {'-',':'};
            [~,tIx] = min(abs(o.time-t));
            plot(o.preferredOrientation,squeeze(o.rate(:,tIx,:)),'linewidth',2)
            xlabel('Preferred Orientation (?)')
            ylabel('Firing rate (spikes/s)')
            xlim([-90 90]);set(gca,'XTick',-90:45:90);
            title('Population Response')
            
            if o.isTwoPopulations
                legend('E','I');
            end
            
            
        end % teichQianPlots
        
        %% Make a video of the timecourse of the population voltage or firing rates
        function timecourseVideo(o,varargin)
            p = inputParser;
            p.addParameter('output','voltage',@(x)(iscell(x) || ismember(x,{'voltage','rate','adaptation','gain'})));
            p.addParameter('method','peak',@(x)(ismember(x,{'peak','centerMass'})));
            p.addParameter('filename','',@ischar); % filename for video  (e.g. .avi)
            p.addParameter('framerate',15,@isnumeric); % Frames per second of video
            p.addParameter('certainIxs',[],@isnumeric);
            p.addParameter('ylim',[],@isnumeric); % If want to match up y axis with other videos
            p.addParameter('showTime',true,@islogical);
            p.addParameter('showLines',true,@islogical);
            p.addParameter('legendLocation','northwest',@ischar); % empty removes the legend
            p.addParameter('timeStep',1,@isnumeric); % Step through o.time with these steps.
            p.addParameter('threeD',false,@islogical);
            p.parse(varargin{:});
            
            if ischar(p.Results.output)
                output = {p.Results.output};
            else
                output = p.Results.output;
            end
            
            var2graph= [];
            for i=1:numel(output)
                switch upper(output{i})
                    case 'VOLTAGE'
                        var2graph = cat(3,var2graph,o.voltage);
                    case 'RATE'
                        var2graph = cat(3,var2graph,o.rate);
                    case 'ADAPTATION'
                        var2graph = cat(3,var2graph,o.adaptation);
                    case 'GAIN'
                        var2graph = cat(3,var2graph,o.gain);
                end
            end
            
            if ~isempty(p.Results.filename)
                vidObj = VideoWriter(p.Results.filename);
                vidObj.FrameRate = p.Results.framerate;
                open(vidObj);
            end
            
            f = findobj('name',['TimeCourseVid ' o.baseModel]);
            if isempty(f)
                f=figure;
                set(f,'name',['TimeCourseVid ' o.baseModel]);
            end
            figure(f);
            clf;
            colors = 'rgbcmy';
            
            if ~isempty(p.Results.ylim)
                forcedYlim = p.Results.ylim;
                forcedZlim = p.Results.zlim; %#ok<NASGU>
                scale =1;
            elseif size(var2graph,3)==1
                maxVal = max(var2graph(:));
                minVal = min(0,min(var2graph(:)));
                if maxVal==minVal
                    minVal = minVal-1;
                    maxVal = maxVal+1;
                end
                forcedYlim = [minVal maxVal+maxVal/10];
                forcedZlim = [minVal maxVal+maxVal/10]; %#ok<NASGU>
                scale = 1;
            else
                maxVal=1;
                minVal =-1;
                forcedYlim = [minVal maxVal+maxVal/10];
                forcedZlim = [minVal maxVal+maxVal/10]; %#ok<NASGU>
                scale = squeeze(max(max(var2graph)));
            end
            h = nan(1,size(var2graph,3));
            offset = min(var2graph(:)); % USed to make sure all scales are positive for threeD
            for i=1:p.Results.timeStep:length(o.time)
                stimOri = o.stimulusSequence(o.time(i));
                stimOri = o.preferredOrientation(stimOri>0); % Map input back to orientation
                for j=1:size(var2graph,3)
                    if p.Results.threeD
                        minSize = 36;
                        maxSize = 720;
                        thisScale = minSize+ (maxSize-minSize)*(var2graph(:,i,j)-offset)./scale(j) ;
                        h(j) =scatter(j*cosd(2*o.preferredOrientation),j*sind(2*o.preferredOrientation),thisScale,['.' colors(j)]);
                        hold on
                        % Current stimulus
                        compass(cosd(2*stimOri),sind(2*stimOri),'k')                       
                    else
                        h(j) = plot(o.preferredOrientation,var2graph(:,i,j)./scale(j),colors(j),'LineWidth',6);hold on
                        plot(o.preferredOrientation(p.Results.certainIxs),var2graph(p.Results.certainIxs,i,j)./scale(j),[colors(j) '.'],'markersize',50);
                        % Current stimulus
                        for thisOri = stimOri
                            if ~isempty(thisOri)
                                plot([ thisOri thisOri],forcedYlim,'k--','lineWidth',4);
                            end
                        end
                    end
                end                
               
                if isempty(stimOri)
                    titleStr = 'Population Dynamics Stimulus OFF';
                else
                     titleStr  = 'Population Dynamics';
                    titleStr = titleStr + string([' Stimulus = ' num2str(round(mod(stimOri(:)',180))) '\circ']);
                end
                if p.Results.showTime
                    titleStr= titleStr + string([' t = ' num2str(round(o.time(i))) ' ms']);
                end
                
                title(titleStr)
                
                
                xlabel('Preferred orientation (\circ)');
                
                if p.Results.threeD
                    ylabel('Preferred orientation (\circ)');
                    ylim([-1 1]*size(var2graph,3));xlim([-1 1]*size(var2graph,3));
                    
                    for stimOri =0:45:179
                        text(size(var2graph,3)*cosd(2*stimOri),size(var2graph,3)*sind(2*stimOri),[num2str(stimOri) '\circ'])
                    end
                    zlabel(p.Results.output)
                    zlim(forcedYlim);
                    set(gca,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','');
                    axis equal; axis square; axis off
                else
                    xlim([-90 90]);
                    set(gca,'XTick',-90:45:90);
                    plot(xlim,[0 0],'k--','LineWidth',2);
                    ylabel(p.Results.output)
                    ylim(forcedYlim);
                end
                if ~isempty(p.Results.legendLocation)
                    legend(h,output,'location',p.Results.legendLocation)
                end
                
                drawnow;               
                if ~isempty(p.Results.filename)
                    currFrame = getframe;
                    writeVideo(vidObj,currFrame);
                end
                hold off
            end
 
            if ~isempty(p.Results.filename)
                close(vidObj);
            end
        end % timecourseVid
        
        
    end
    
    methods (Access= protected)
        
        function [input] = inputFromPrecompute(o,t)
            % precomputed stores t==0 as the first entry, before t ==0
            % input is asumed to be 0.
            % First row is orientation, second row the contrast (defaults
            % to 1).
            phi         = zeros(size(t));
            contrast    = ones(size(t));
            [nrVars,nrTimes] = size(o.preComputedSequence);
            ix = floor(t)+1; % First entry in .preComputed corresponds to t==0, second to t==1
            ix(ix>nrTimes) = nrTimes; % The last stimulus stays "on"
            % Before t==0 zero contrast
            phi(t<0) = 0;
            contrast(t<0) = 0;
            phi(t>=0) = o.preComputedSequence(1,ix(t>=0));
            if nrVars>1
                % Contrast has been precomputed
                contrast(t>0) = o.preComputedSequence(2,ix(t>0));
            end
            % Now map to single grating input
            input = singleGratingToInput(o,phi,contrast);
        end
        
        function o = shrink(o)
            % Shring the lateral interactions (without affecting the peak
            % strength)
            o.excitatoryWeights =shrnk(o.preferredOrientation,o.excitatoryWeights,o.shrinkFactor(1));
            o.inhibitoryWeights =shrnk(o.preferredOrientation,o.inhibitoryWeights,o.shrinkFactor(2));
            if o.isTwoPopulations
                o.cortexWeights = o.cortexStrength.*[o.excitatoryWeights(:,1) -o.ieRatio(1).*o.inhibitoryWeights(:,1) o.excitatoryWeights(:,2) -o.ieRatio(2).*o.inhibitoryWeights(:,2)]; % ee, ei, ie, ii;
            else
                o.cortexWeights = o.cortexStrength.*(o.excitatoryWeights-o.ieRatio.*o.inhibitoryWeights);
            end
            % Nested function to do the actual shrinking iwht interpolation
            function trg = shrnk(x,y,factor)
                trgX = min(x):factor:max(x);
                trg = interp1(x,y,trgX);
                trg = interp1(trgX./factor,trg,x);
                for j=1:size(trg,2)
                    trg(isnan(trg(:,j)),j) =min(trg(:,j));
                end
            end
            
        end
        
        function [v,a,g] = vag(o)
            % Voltage Adaptation Gain state variable
            if o.isTwoPopulations
                v = [1 2];
                a = [3 4];
                g = [5 6];
            else
                v =1;a =2; g =3;
            end
        end
        
        function stimulus = singleGratingToInput(o,ori,contrast)
            % Convenience function to map a ori and contrast value to the
            % networks input pattern. This can also be used to present a
            % single grating ot the network with s.stimulusSequece
            nrOris = numel(ori);
            stimulus = zeros(o.nrNeurons,nrOris);
            for i=nrOris
                if ~isnan(ori(i))
                    [~,ix] = min(abs(mod(o.preferredOrientation,180)-mod(ori(i),180)));
                    stimulus(ix,i) = contrast(i);
                end
            end
        end
        
        
    end
    
    
    
    methods (Static)
        
        
        function v= warpCosd(x,factor)
            % Warped cosine.
            % factor= 0 is standard cosine
            % Positive factors make the cosine more narrow.
            % factor ~4  = cosine with a width of 90
            % negative factors broaden the cosine
            if factor>=0
                f=factor+1;
                v=(cosd(x)+1).^f/(2^(f-1))-1;
            else
                f=abs(factor)+1;
                v=-(cosd(x+180)+1).^f/(2^(f-1))+1;
            end
        end
        
        function seq = adaptSequence(t, adapterOri, testOri, adaptDuration, testDuration,blankDuration,adapterContrast,testContrast)
            %  A standard adaptation sequence
            if nargin <8
                testContrast =1;
                if nargin < 7 
                    adapterContrast =1;
                    if nargin <6
                    blankDuration = 0;
                    end
                end
            end
            v = nan(1,numel(t));
            c = ones(1,numel(t));
            c(t<0) = 0;
            isAdaptPhase= t<=adaptDuration;
            v(isAdaptPhase) = adapterOri;
            c(isAdaptPhase) = adapterContrast;
            
            c(t>adaptDuration & v <= adaptDuration+blankDuration) = 0; % blank            
            
            isTestPhase  =(t>adaptDuration+blankDuration & t<= adaptDuration+blankDuration + testDuration);
            v(isTestPhase) = testOri;
            c(isTestPhase) = testContrast;
            
            c(t>adaptDuration+blankDuration + testDuration) = 0;
            
            seq = [v;c];
        end
        
        
        function y = vonMises(x,mu,k,period,lambda)
            % The vonMises function 
            if nargin==3
                period = 360;
                lambda = 0;
            end
            if nargin==4
                lambda = 0;
            end
            y = (exp(k*cos(2*pi*(x-mu)/period))./(2*pi*besseli(0,k))).*(1+lambda*sin(2*pi*(x-mu)/period));
        end
        
        
        
        function [rate,ph] = retinalResponse(freq,contrast,varargin)
            % Implementation of the Victor (1978) model based on Bernadette& Kaplan (1990) primate
            % retinal ganglion cell temporal response measurements.
            % Default values are the median values for M-cell responses (ignoring ON/OFF differences).
            % freq = frequency in hertz.
            % contrast = contrast. 1 means 100% contrast
            %
            % This can be used to generate more realistic retinal input to
            % the model.         
            p=inputParser;
            p.addParameter('gain',500,@isnumric); % gain per unit contrast  (A)
            p.addParameter('delay',2.33,@isnumeric); % initial delay in ms (D)
            p.addParameter('highPass',1,@isnumeric); % Hs- strength of high pass filter.  (Hs)
            p.addParameter('tauLow',1.68,@isnumeric); % Time-constant of low pass filter.    (tl)
            p.addParameter('tauHighZero',45,@isnumeric);%time constant of high pass stage at zero contrast (T0)
            p.addParameter('c50',0.048,@isnumeric); % Contrast at which tauLow is half its initial value. (C1/2)
            p.addParameter('nLow',25.5,@isnumeric); % Number of low-pass filters  (Nl)
            p.addParameter('plot',false);
            p.parse(varargin{:});
            
            tauHigh = p.Results.tauHighZero./(1+(contrast/p.Results.c50).^2);
            omega  = 2*pi*freq/1000; % mHz
            i = sqrt(-1);
            R  = p.Results.gain.*exp(-i.*omega.*p.Results.delay).*...
                (1-p.Results.highPass./(1+i.*omega.*tauHigh)).*...
                (1./(1+i.*omega.*p.Results.tauLow)).^p.Results.nLow;
            rate = abs(R);
            ph   = phase(R); % Phase advance for this "median" magno cell is much larger than the examples in Bernadette data.
            
            
            if p.Results.plot~=0
                if ishandle(p.Results.plot)
                    axes(p.Results.plot);
                end
                yyaxis left
                plot(freq,rate)
                ylabel('Response (Hz)');
                hold on
                yyaxis right
                plot(freq,ph);
                xlabel 'Frequency (Hz)'
                ylabel('Phase');
                set(gca,'XScale','log')
                hold off
            end
        end
    end
end %