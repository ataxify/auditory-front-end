classdef hlProc < Processor
%HLPROC Hearing Level processor.
%   This processor estimates the Hearing Level (HL) either based on 
%   the absolute threshold of hearing  ('ath') or on the energy in each filter ('en').
%   The thresholds are computed for each filter of the auditory filterbank.
%
%   HLPROC properties:
%       Method    - Method used to estimate the thresholds 
%                     absolute threshold of hearing ('ath') or energy-based ('en')
%       enTh      - Threshold for energy-based method
%       wname       - Window shape descriptor
%       wSizeSec    - Window duration in seconds
%       hSizeSec    - Step size between windows in seconds
%
%   See also: Processor, ihcProc
%
    
    properties (Dependent = true)
        Method    % Method used to estimate the thresholds 
                    % absolute threshold of hearing ('ath') or energy-based ('en')
        enTh      % Threshold for energy-based method
        wname       % Window shape descriptor (see window.m)
        wSizeSec    % Window duration in seconds
        hSizeSec    % Step size between windows in seconds
    end
    
    properties (GetAccess = private)
        wSize       % Window duration in samples
        hSize       % Step size between windows in samples
        win         % Window vector
        buffer_l    % Buffered input signals (left ear)
        buffer_r    % Buffered input signals (right ear)
        cfHz        % Center frequencies of filterbank
    end
    
    methods
        function pObj = hlProc(fs,parObj)
		%ildProc   Construct an hearing level extractor processor
        %
        % USAGE:
        %   pObj = hlProc(fs, parObj)
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
            
            % Checking input parameter
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,[],'hlProc',parObj);
            
            if nargin>0
                % Initializa the buffers
                pObj.buffer_l = [];
                pObj.buffer_r = [];
            end
            
        end
        
        function [out_l, out_r] = processChunk(pObj,in_l,in_r)
            %processChunk       Apply the processor to a new chunk of input signal
            %
            %USAGE
            %   [out_l out_r] = pObj.processChunk(pObj,in_l,in_r)
            %
            %INPUT ARGUMENT
            %    pObj : Processor instance
            %    in_l : New chunk of input data from left channel
            %    in_r : New chunk of input data from right channel
            %
            %OUTPUT ARGUMENT
            %   out_l : Corresponding output (left channel)
            %   out_r : Corresponding output (right channel)
            
            % Append provided input to the buffer
            if ~isempty(pObj.buffer_l)
                in_l = [pObj.buffer_l;in_l];
                in_r = [pObj.buffer_r;in_r];
            end
            
            % Quick control of dimensionality
            if max(size(in_l)~=size(in_r))
                error('Buffered inputs should be of same dimension for both ears')
            end
            
            [nSamples,nChannels] = size(in_l);
            
            % How many frames are in the buffered input?
            nFrames = floor((nSamples-(pObj.wSize-pObj.hSize))/pObj.hSize);
            
            % Compute HL:
            
            %init SPL calculations
            frame_l_SPL = zeros(1,nChannels);
            frame_r_SPL = zeros(1,nChannels);
            
            % Pre-allocate outputs
            out_l = NaN(nSamples,nChannels);
            out_r = NaN(nSamples,nChannels);
            
            % Window processing
            winRep = repmat(pObj.win,1,nChannels);
            
            % Get threshold
            switch pObj.Method
                case 'none'
                    actbands = (1:length(pObj.cfHz)); %column-vector with #filters
                    
                case 'ath'
                    ath = calcATH(pObj.cfHz*1e-3); %absolute threshold of hearing
                    xref = 20e-5; %ref for dB SPL definition: xrms = 1 -> Lp = 74dB SPL
                    %(xrms = xref = 20e-5 -> Lp = 0dB SPL)
                    % xref = 1e-5; %ref for digital value definition: xrms = 1 -> Lp = 100dB
                    %              %using this reference will reproduce totEnergy = Lp
                    %              (Lp being the specified level of the input signal)
                    
                    % Loop on the time frame
                    for ii = 1:nFrames
                        % Get start and end indexes for the current frame
                        n_start = (ii-1)*pObj.hSize+1;
                        n_end = (ii-1)*pObj.hSize+pObj.wSize;
                        
                        % Energy in the windowed frame for left and right input
                        frame_l_rms = sqrt(mean(power(winRep.*in_l(n_start:n_end,:),2)));
                        frame_r_rms = sqrt(mean(power(winRep.*in_r(n_start:n_end,:),2)));
                        
                        %update SPL and devide by number of frames for averaging
                        %between two successive frames
                        frame_l_SPL = 20*log10(0.5*(frame_l_rms/xref + 10.^(frame_l_SPL/20))); %acc. SPL left channel
                        frame_r_SPL = 20*log10(0.5*(frame_r_rms/xref + 10.^(frame_r_SPL/20))); %acc. SPL right channel
                        
                    end
                    
                    %select maximum values from left and right channel
                    max_lr_SPL = max(frame_l_SPL,frame_r_SPL);
                    
                    %all filters that are above absolute threshold of hearing
                    bpassth = max_lr_SPL>ath; %bool 
                    
                    %active bands
                    actbands = (1:length(pObj.cfHz)); %column-vector with #filters
                    actbands = actbands(bpassth); %select activated channels from the vector
                    
                case 'en'
                    warning('Energy-based method not implemented yet!')
                    
                otherwise
                    error(['The method ' pObj.Method ' is not defined!'])
            end
            
            % Set all frames to 0 that are below the specified threshold
            out_l(:,actbands) = in_l(:,actbands);
            out_r(:,actbands) = in_r(:,actbands);
            
            %frame_SPL contains total SPL after frame-based processing
            %(check if it corresponds to the specified level)
            totSPL = 10*log10(sum(10.^(max_lr_SPL/10))); %tot. level from all filters left channel
            
            % Update the buffer: the input that was not extracted as a
            % frame should be stored
            pObj.buffer_l = in_l(nFrames*pObj.hSize+1:end,:);
            pObj.buffer_r = in_r(nFrames*pObj.hSize+1:end,:);
            
            
        end
        
        function reset(pObj)
             %reset     Resets the internal states of the ILD extractor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : ILD extractor processor instance
             
             % Only thing needed to reset is to empty the buffer
             pObj.buffer_l = [];
             pObj.buffer_r = [];
             
        end
        
    end
    
    methods (Hidden = true)
        
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for the pre-processor
            %
            %NB: This method is overloading the parent method to deal with multiple
            %outputs
            
            if dObj.isStereo
                sig_l = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            'left');
                sig_r = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            'right');
                dObj.addSignal(sig_l);
                dObj.addSignal(sig_r);
            
                output = {sig_l, sig_r};
                
            else
                sig = feval(pObj.getProcessorInfo.outputType, ...
                            pObj, ...
                            dObj.bufferSize_s, ...
                            pObj.Channel);
            
                dObj.addSignal(sig);
            
                output = {sig};
                
            end
            
        end
        
        function initiateProcessing(pObj)
            %INITIATEPROCESSING    Wrapper calling the processChunk method and routing I/O
            % Because the pre-processor can have two outputs (for stereo signals), it is
            % necessary to overload the parent method here.
            
            
            if size(pObj.Input,2)>1
                [out_l, out_r] = pObj.processChunk( pObj.Input{1,1}.Data('new'),...
                    pObj.Input{1,2}.Data('new'));
            else
                [out_l, out_r] = pObj.processChunk( pObj.Input{1,1}.Data('new'),...
                    []);
            end
            
            pObj.Output{1}.appendChunk(out_l);
            
            if ~isempty(out_r)
                pObj.Output{2}.appendChunk(out_r);
            end
            
        end
        
        function prepareForProcessing(pObj)
            
            % Compute internal parameters
            pObj.wSize = 2*round(pObj.parameters.map('hl_wSizeSec')*pObj.FsHzIn/2);
            pObj.hSize = round(pObj.parameters.map('hl_hSizeSec')*pObj.FsHzIn);
            pObj.win = window(pObj.parameters.map('hl_wname'),pObj.wSize);
            % Output sampling frequency
%             pObj.FsHzOut = 1/(pObj.hSizeSec);
            pObj.FsHzOut = pObj.FsHzIn;
            
            % Access center frequencies for computing acitve bands
            pObj.cfHz = pObj.getDependentParameter('fb_cfHz');
                
        end
        
    end
    
    % "Getter" methods
    methods
        function Method = get.Method(pObj)
            Method = pObj.parameters.map('hl_Method');
        end
        
        function enTh = get.enTh(pObj)
            enTh = pObj.parameters.map('hl_enTh');
        end
        
        function wSizeSec = get.wSizeSec(pObj)
            wSizeSec = pObj.parameters.map('hl_wSizeSec');
        end
        
        function hSizeSec = get.hSizeSec(pObj)
            hSizeSec = pObj.parameters.map('hl_hSizeSec');
        end
        
        function wname = get.wname(pObj)
            wname = pObj.parameters.map('hl_wname');
        end
        
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'innerhaircell';
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
            
            
            names = {'hl_Method',...
                    'hl_enTh',...
                    'hl_wname',...
                    'hl_wSizeSec',...
                    'hl_hSizeSec'};
            
            descriptions = {'Hearing Level Method',...
                    'Threshold for energy-estimation',...
                    'Window name',...
                    'Window duration (s)',...
                    'Window step size (s)'};
            
            defaultValues = {'ath',...
                            0.5,...
                            'hann',...
                            20E-3,...
                            10E-3};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'HL Extractor';
            pInfo.label = 'Hearing Level Extractor';
            pInfo.requestName = 'hl';
            pInfo.requestLabel = 'Hearing Level';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = 2;
            pInfo.hasTwoOutputs = 1;
            
        end
        
    end
        
    
end