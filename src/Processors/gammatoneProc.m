classdef gammatoneProc < Processor
    
    properties
        cfHz            % Filters center frequencies
        Filters         % Array of filter objects
        nERBs           % Distance between neighboring filters in ERBs
        n_gamma         % Gammatone order of the filters
        bwERBs          % Bandwidth of the filters in ERBs
        f_low           % Lowest center frequency used at instantiation
        f_high          % Highest center frequency used at instantiation
        fb_decimation   % Decimation ratio of the filterbank
    end
        
    methods
        function pObj = gammatoneProc(fs,flow,fhigh,irType,nERBs,bAlign,...
                                        n,bw,durSec)
            %gammatoneProc      Construct a gammatone filterbank inheriting
            %                   the "processor" class
            %
            %USAGE
            %   pObj = gammatoneProc(fs,flow,fhigh)
            %   pObj = gammatoneProc(fs,flow,fhigh,irType,nERBs,bAlign)
            %   pObj = gammatoneProc(fs,flow,fhigh,irType,nERBs,bAlign,n,bw,dur)
            %
            %INPUT ARGUMENTS
            %     fs : Sampling frequency (Hz)
            %   flow : Lowest center frequency for the filterbank (in Hz)
            %  fhigh : Highest center frequency for the filterbank (in Hz)
            % irType : 'FIR' to generate finite impulse response Gammatone
            %          filters or 'IIR' for infinite (default: 'FIR')
            %  nERBs : Distance in ERBS between neighboring center
            %          frequencies (default: nERBS = 1)
            % bAlign : Set to true for phase correction and time alignment
            %          between channels (default: bAlign = false)
            %      n : Filter order (default: n = 4)
            %     bw : Bandwidth of the filters in ERBS 
            %          (default: bw = 1.08 ERBS)
            %    dur : Duration of the impulse response in seconds 
            %          (default: dur = 0.128)
            %
            %OUTPUT ARGUMENTS
            %   pObj : Processor object
            %
            % TO DO: Implement solution to allow for different impulse
            % response durations for different filters (if necessary)
            
            if nargin>0  % Failsafe for constructor calls without arguments
            
            % Checking input arguments
            if nargin < 3 || nargin > 9
                help(mfilename);
                error('Wrong number of input arguments!')
            end
            
            % Set default parameter
            if nargin < 9 || isempty(durSec); durSec = 0.128; end
            if nargin < 8 || isempty(bw); bw = 1.08; end
            if nargin < 7 || isempty(n); n = 4; end
            if nargin < 6 || isempty(bAlign); bAlign = false; end
            if nargin < 5 || isempty(nERBs);  nERBs  = 1;     end
            if nargin < 4 || isempty(irType); irType = 'FIR'; end
            
            % ERBs vector, with a spacing of nERBs
            ERBS = freq2erb(flow):double(nERBs):freq2erb(fhigh);
            
            % Conversion from ERB to Hz
            cfHz = erb2freq(ERBS);
            
            % Number of gammatone filters
            nFilter = numel(cfHz); 
            
            % Instantiating the filters
            pObj.Filters = pObj.populateFilters(cfHz,fs,irType,n,bw,...
                            bAlign,durSec);
            
            % Setting up additional properties
            % 1- Global properties
            populateProperties(pObj,'Type','Gammatone filterbank',...
                'Dependencies',getDependencies('gammatone'),...
                'FsHzIn',fs,'FsHzOut',fs);
            % 2- Specific properties
            pObj.cfHz = cfHz;
            pObj.nERBs = nERBs;
            pObj.n_gamma = n;
            pObj.bwERBs = bw;
            pObj.f_low = flow;
            pObj.f_high = fhigh;
            pObj.fb_decimation = 1;
            
            end
        end
        
        function out = processChunk(pObj,in)
            %processChunk       Passes an input signal through the
            %                   Gammatone filterbank
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : Gammatone filterbank object
            %        in : One-dimensional array containing the input signal
            %
            %OUTPUT ARGUMENTS
            %       out : Multi-dimensional array containing the filterbank
            %             outputs
            %
            %SEE ALSO:
            %       gammatoneProc.m
            
            % TO DO: Indicate that this function is not buit to deal with
            % multiple channels. Multiple channels should be treated with
            % multiple instances of the filterbank.
            
            % Check inputs
            if min(size(in))>1
                error('The input should be a one-dimensional array')
            end
            
            % Turn input into column vector
            in = in(:);
            
            % Get number of channels
            nFilter = size(pObj.Filters,2);
            
            % Pre-allocate memory
%             out = zeros(nFilter,floor(size(in,2)*pObj.FsHzOut/pObj.FsHzIn));
            
            % Loop on the filters
            for ii = 1:nFilter
                out(:,ii) = pObj.Filters(ii).filter(in);
            end
            
            % TO DO : IMPLEMENT ALIGNMENT CORRECTION
            
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
            
            nFilter = size(pObj.Filters,2);
            
            % Resetting the internal states of the filters
            for ii = 1:nFilter
                pObj.Filters(ii).reset();
            end
            
        end
        
        function hp = hasParameters(pObj,p)
            %hasParameters  This method compares the parameters of the
            %               processor with the parameters given as input
            %
            %USAGE
            %    hp = pObj.hasParameters(p)
            %
            %INPUT ARGUMENTS
            %  pObj : Processor instance
            %     p : Structure containing parameters to test
            
            %NB: Could be moved to private
            
            % The gammatone processor has the following parameters to be
            % checked: f_low, f_high, nERBs, n, bwERBs, FsHzIn, FsHzOut
            
            p_list = {'f_low','f_high','nERBs','n_gamma','bwERBs'};%,...
                %'fb_decimation'};
            
            % Initialization of a parameters difference vector
            delta = zeros(size(p_list,2),1);
            
            % Loop on the list of parameters
            for ii = 1:size(p_list,2)
                try
                    delta(ii) = abs(pObj.(p_list{ii}) - p.(p_list{ii}));
                    
                catch err
                    % Warning: something is missing
                    warning('Parameter %s is missing in input p.',p_list{ii})
                    delta(ii) = 1;
                end
            end
            
            % Check if delta is a vector of zeros
            if max(delta)>0
                hp = false;
            else
                hp = true;
            end
            
        end  
        
    end
    
    methods (Access = private)
        function obj = populateFilters(pObj,cfHz,fs,irType,n,bw,bAlign,durSec)
            % This function is a workaround to assign an array of objects
            % as one of the processor's property, should remain private

            nFilter = numel(cfHz);
            
            % Preallocate memory by instantiating last filter
            obj(1,nFilter) = gammatoneFilter(cfHz(nFilter),fs,irType,n,...
                                        bw,bAlign,durSec);
            % Instantiating remaining filters
            for ii = 1:nFilter-1
                obj(1,ii) = gammatoneFilter(cfHz(ii),fs,irType,n,...
                                        bw,bAlign,durSec);
            end                        
            
        end
    end
        
end