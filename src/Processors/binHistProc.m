classdef binHistProc < Processor
%BINHISTPROC Binaural Histogram estimation processor.
%   This processor calculates a histogram [container x frequency] based on a binaural cue ('itd' or 'ild') matrix [time x frequency].
%
%   BINHISTPROC properties:
%       method   - Select binaural cue ('itd' or 'ild')
%       nbins    - #equally spaced containers in histogram calculation (default: nbins = 41)
%       brel - Boolean for displaying relative occurance (normalise data to number of frames)',...
%       blatcomp - Boolen for lateral compression
%       chThreshold - Select thresholding method for frequency channels ('none','hearing','energy') 
%       freqWeighting - Select frequency-dependent weighting of the binaural cue
%
%   See also: Processor, itdProc, ildProc
%

%% Properties

    properties (Dependent = true)
        method   % Select binaural cue ('itd' or 'ild')
        nbins    % #equally spaced containers in histogram calculation (default: nbins = 41)
        brel % Boolean for displaying relative occurance (normalise data to number of frames)',...
        blatcomp % Boolen for lateral compression
        chThreshold % Select thresholding method for frequency channels ('none','hearing','energy') 
        freqWeighting % Select frequency-dependent weighting of the binaural cue
    end
    
    properties (GetAccess = private)
        ymax %Defines the bin edges for the edges vector in the histogram estimation
        
    end
        
    % You can add more properties with different attributes if needed
    
%% Methods
    methods
        function pObj = binHistProc(fs,parObj)
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
            pObj = pObj@Processor(fs,fs,'binHistProc',parObj);
            
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
            
            %get active bands
            actbands = getActBands(in);
            
            %ymax
%             pObjPar = getCurrentParameters(pObj,0);
%             dependentProcs = {'itd'};
%             addLowerDependencies(pObj,dependentProcs)
            
            %lateral compression
            if pObj.blatcomp
                in = latCompression(in,pObj.ymax);
            end
            
            %Nhist
            Nhist = binCueHistc(in,'nbins',pObj.nbins,'brel',pObj.brel,...
                    'actbands',actbands,'ymax',pObj.ymax);

            out = Nhist;
            
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
            
            switch pObj.method
                case 'itd'
                    pObj.ymax = 1e-3; %1 ms
                    
                case 'ild'
                    pObj.ymax = 20; %10 dB
                
                case 'ic'
                    pObj.ymax = 1; %correlation
                    
                otherwise
                    error(['The function is not specified for signal ' pObj.method '!'])
                
            end
                
        end
        
    end
    
    %% "Getter" methods
    methods
        function method = get.method(pObj)
            method = pObj.parameters.map('bHist_method');
        end
        function nbins = get.nbins(pObj)
            nbins = pObj.parameters.map('bHist_nbins');
        end
        function brel = get.brel(pObj)
            brel = pObj.parameters.map('bHist_brel');
        end
        function blatcomp = get.blatcomp(pObj)
            blatcomp = pObj.parameters.map('bHist_blatcomp');
        end
        function chThreshold = get.chThreshold(pObj)
            chThreshold = pObj.parameters.map('bHist_chThreshold');
        end
        function freqWeighting = get.freqWeighting(pObj)
            freqWeighting = pObj.parameters.map('bHist_freqWeighting');
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
            
            
            names = {'bHist_method',...
                     'bHist_nbins',...
                     'bHist_brel',...
                     'bHist_blatcomp',...
                     'bHist_chThreshold',...
                     'bHist_freqWeighting',...
                    };
            
            descriptions = {'Select binaural cue (''itd'' or ''ild'')',...
                            '#equally spaced containers in histogram calculation',...
                            'Boolean for displaying relative occurance (normalise data to number of frames)',...
                            'Boolen for lateral compression',...
                            'Select thresholding method for frequency channels (''none'',''hearing'' or ''energy'')',... 
                            'Select frequency-dependent weighting of the binaural cue',...
                            };
            
            defaultValues = {'itd',...
                             41,...
                             true,...
                             false,...
                             'none',...
                             'none',...
                            };
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Binaural Histogram';
            pInfo.label = 'Binaural Histogram Estimation';
            pInfo.requestName = 'binHist';
            pInfo.requestLabel = 'Binaural Histogram';
            pInfo.outputType = 'BinHistogramSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
        
end