classdef aswProc < Processor
%ASWPROC Apparent Source Width (ASW) prediction processor.
%   This processor predicts apparent source width perception.
%
%   ASWPROC properties:
%       binCue   - Select binaural cue ('itd' or 'ild')
%       wdMethod - Method for estimating the distribution width ('prct' or 'std')
%       percent  - Choose percent values for prctile calculation
%       avMethod - Method for averaging across frequency channels
%       transform - 'Choose transformation of bin. cue: 'none', 'norm' (normalize)  or 'latcomp' (lateral compression)'
%       freqWeighting - Select frequency-dependent weighting of the binaural cue
%
%   See also: Processor, itdProc, ildProc
%

%% Properties

    properties (Dependent = true)
        binCue   % Select binaural cue ('itd' or 'ild')
        wdMethod % Method for estimating the distribution width ('prct' or 'std')
        percent  % Choose percent values for prctile calculation
        avMethod % Method for averaging across frequency channels
        transform % 'Choose transformation of bin. cue: ''none'', ''norm'' (normalize)  or ''latcomp'' (lateral compression)'
        freqWeighting % Select frequency-dependent weighting of the binaural cue
    end
    
    properties (GetAccess = private)
        ymax %Defines the bin edges for the edges vector in the histogram estimation
        
    end
        
    % You can add more properties with different attributes if needed
    
%% Methods
    methods
        function pObj = aswProc(fs,parObj)
        %aswProc   Construct an apparent source width (ASW) prediction processor
        %
        % USAGE:
        %   pObj = aswProc(fs, parObj)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        %
        % OUTPUT ARGUMENTS:
        %   pObj : Processor instance
        %
        % NOTE: Parameter object instance, parObj, can be generated using genParStruct.m
        % User-controllable parameters for this processor and their default values can be
        % found by browsing the script parameterHelper.m
        %
        % See also: genParStruct, parameterHelper, Processor
            
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call super-constructor
            pObj = pObj@Processor(fs,fs,'aswProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                % Add here any additional initialization steps which need only be called
                % once at instantiation. Additional initialization steps that should be
                % performed again after receiving feedback are implemented below
            end
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Requests the processing for a new chunk of
            %                   signal
            %
            %USAGE:
            %    out = processChunk(in)
            %
            %INPUT ARGUMENTS:
            %   pObj : Processor instance
            %     in : Input chunk
            %
            %OUTPUT ARGUMENT:
            %    out : Processor output for that chunk
            
            %% Transformation of bin. cue
            switch pObj.transform
                case 'none' %no transformation
                    %do nothing
                    
                case 'norm' %normalization
                    in = in/pObj.ymax;
                    
                case 'latcomp' %lateral compression
                    in = latCompression(in,pObj.ymax);
                    
                otherwise
                    error(['The transformation method ' pObj.transform ' is not defined!'])
          
            end
            
%             pObj.LowerDependencies;
%             removeUpperDependency(pObj,{'asw'})
%             removeHandleInLowerDependencies(pObj)
            
            %% width
            switch pObj.avMethod
                case 'wdthFirst'
                    %A) calculate width per channel first, then average across channels
                    widthChan = calcDistrWidth(in,pObj.wdMethod,pObj.percent); %calculate width of binCue directly for each frequency channel
                    
                    switch pObj.freqWeighting
                        case 'none'
                            width = mean(widthChan,2); %total width as averaging across all channels
                        otherwise
                            error('This frequency weighting is not defined!')
                    end
                    
                case 'avFirst'
                    %B) average channels first, then calculate width
                    switch pObj.freqWeighting
                        case 'none'
                            meanChan = mean(in,2); %total width as averaging across all channels
                        otherwise
                            error('This frequency weighting is not defined!')
                    end
                    width = calcDistrWidth(meanChan,pObj.wdMethod,pObj.percent);
                    
                otherwise
                    error('Averaging method is not defined!')
            end
            
            out = width;
            
        end
        
        function reset(pObj)
            %reset          Order the processor to reset its internal
            %               states, e.g., when some critical parameters in
            %               the processing have been changed
            %USAGE
            %       pObj.reset()
            %       reset(pObj)
            %
            %INPUT ARGUMENT
            %       pObj : Processor object
            
            % Reset the internal states of your processor here, if any. 
            % Otherwise, leave empty.
            
        end
        
    end
    
    %% Hidden Methods
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            switch pObj.binCue
                case 'itd'
                    pObj.ymax = 1e-3; %1 ms
                    
                case 'ild'
                    pObj.ymax = 20; %20 dB
                
                case 'ic'
                    pObj.ymax = 1; %correlation
                    
                otherwise
                    error(['The function is not specified for signal ' pObj.method '!'])
                
            end
                
        end
        
    end
    
    %% "Getter" methods
    methods
        function binCue = get.binCue(pObj)
            binCue = pObj.parameters.map('asw_binCue');
        end
        function wdMethod = get.wdMethod(pObj)
            wdMethod = pObj.parameters.map('asw_wdMethod');
        end
        function percent = get.percent(pObj)
            percent = pObj.parameters.map('asw_percent');
        end
        function avMethod = get.avMethod(pObj)
            avMethod = pObj.parameters.map('asw_avMethod');
        end
        function transform = get.transform(pObj)
            transform = pObj.parameters.map('asw_transform');
        end
        function freqWeighting = get.freqWeighting(pObj)
            freqWeighting = pObj.parameters.map('asw_freqWeighting');
        end
        
    end
    
    %% Static methods
    % These methods store all hard-coded information regarding your processor, remember to
    % update them with actual information!
    methods (Static)
        
        function dep = getDependency()
            dep = 'itd';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ...
            %                           gammatoneProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'asw_binCue',...
                     'asw_wdMethod',...
                     'asw_percent',...
                     'asw_avMethod',...
                     'asw_transform',...
                     'asw_freqWeighting'
                    };
            
            descriptions = {'Select binaural cue (''itd'' or ''ild'')',...
                            'Method for estimating the distribution width (''prct'' or ''std'')',...
                            'Choose percent values for prctile calculation',...
                            'Method for averaging the distribution across channels (''wdthFirst'' or ''avFirst'')',...
                            'Choose transformation of bin. cue: ''none'', ''norm'' (normalize)  or ''latcomp'' (lateral compression)',...
                            'Select frequency-dependent weighting of the binaural cue',...
                            };
            
            defaultValues = {'itd',...
                             'prct',...
                             [10 90],...
                             'wdthFirst',...
                             'none',...
                             'none',...
                            };
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'ASW Predictor';
            pInfo.label = 'Apparent Source Width Predictor';
            pInfo.requestName = 'asw';
            pInfo.requestLabel = 'Apparent Source Width';
            pInfo.outputType = 'TimeDomainSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
        
end