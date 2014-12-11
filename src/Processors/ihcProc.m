classdef ihcProc < Processor
    
     properties (Dependent = true)
         method        % Label for the IHC model used
     end
     
     properties (GetAccess = private)
         IHCFilter     % Filter involved in the extraction, if any
     end
     
     methods
         function pObj = ihcProc(fs,parObj)
             %ihcProc   Construct a inner haircell (IHC) envelope
             %                  extractor
             %
             %USAGE
             %   pObj = ihcProc(fs,method)
             %
             %INPUT ARGUMENTS
             %     fs : Sampling frequency (Hz)
             % method : Envelope extraction method, among 'halfwave',
             %          'fullwave', 'square', 'hilbert', 'joergensen',
             %          'dau', 'breebart', 'berstein'
             %
             %N.B: The constructor does not instantiate the lowpass filters
             %needed for some of the methods.
             
             % TO DO: Detail the help file more 
             
             if nargin>0
             
             
             
             % Populate the object's properties
             % 1- Global properties
             populateProperties(pObj,'Type','IHC envelope extractor',...
                 'Dependencies',getDependencies('innerhaircell'),...
                 'FsHzIn',fs,'FsHzOut',fs);
             % 2- Specific properties
             pObj.method = parObj.map('IHCMethod');
             
             % Store specific parameters
             pObj.parameters = parObj.getProcessorParameters('ihcProc');
             
             % Instantiate a low-pass filter if needed
             switch pObj.method
                 
                 case 'joergensen'
                     % First order butterworth filter @ 150Hz
                     pObj.IHCFilter = bwFilter(pObj.FsHzIn,1,150);

                 case 'dau'
                     % Second order butterworth filter @ 1000Hz
                     pObj.IHCFilter = bwFilter(pObj.FsHzIn,2,1000);

                 case 'breebart'
                     % First order butterworth filter @ 2000Hz cascaded 5
                     % times
                     pObj.IHCFilter = bwFilter(pObj.FsHzIn,1,2000,[],5);

                 case 'bernstein'
                     % Second order butterworth filter @ 425Hz
                     pObj.IHCFilter = bwFilter(pObj.FsHzIn,2,425);

                 otherwise
                     pObj.IHCFilter = [];
             end
             
             end
            
         end
         
         function out = processChunk(pObj,in)
                        
            % Carry out the processing for the chosen IHC method
            switch pObj.method
                case 'none'
                    out = in;

                case 'halfwave'
                    % Half-wave rectification
                    out = max(in,0);

                case 'fullwave'
                    % Full-wave rectification
                    out = abs(in);

                case 'square'
                    out = abs(in).^2;

                case 'hilbert'
                    out = abs(hilbert(in));

                case 'joergensen'
                    out = pObj.IHCFilter.filter(abs(hilbert(in)));

                case 'dau'
                    out = pObj.IHCFilter.filter(max(in,0));

                case 'breebart'
                    out = pObj.IHCFilter.filter(max(in,0));

                case 'bernstein'
                    env = max(abs(hilbert(in)).^(-.77).*in,0).^2;
                    out = pObj.IHCFilter.filter(env);

                otherwise
                    error('%s: Method ''%s'' is not supported!',upper(mfilename),pObj.IHCMethod)
            end
            
         end
         
         function reset(pObj)
             %reset     Resets the internal states of the IHC envelope
             %          extractor
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : Inner haircell envelope extractor processor instance
             
             % A reset is needed only if the extractor involves filters
             if ~isempty(pObj.IHCFilter)
                 pObj.IHCFilter.reset
             end
         end
         
         function verifyParameters(~,parObj)
             % List of valid methods
             validMeth = {'none',...
                         'halfwave',...
                         'fullwave',...
                         'square',...
                         'hilbert',...
                         'joergensen',...
                         'dau',...
                         'breebart',...
                         'bernstein'};
             
             % Check method name
             if ~ismember(parObj.map('IHCMethod'),validMeth)
                 [~,defaultMethod] = ihcProc.getParameterInfo;
                 warning(['''%s'' is an invalid name for envelope extraction method. '...
                          'Setting it to the default value, ''%s'''],...
                          parObj.map('IHCMethod'),defaultMethod)
                 parObj.map('IHCMethod') = defaultMethod;
             end
             
         end
         
     end
         
     methods (Static)
        function dep = getDependency()
            dep = 'gammatone';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ihcProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'IHCMethod'};
            
            descriptions = {['Inner hair-cell envelope extraction method (''none'', ' ...
                            '''halfwave'', ''fullwave'', ''square'', ''hilbert'', '...
                            '''joergensen'', ''dau'', ''breebart'', ''berstein'')']};
            
            defaultValues = {'dau'};
                
        end
        
        function [name,description] = getProcessorInfo
            
            %Returns a very short name and a short description of the processor function
            name = 'IHC envelope';
            description = 'Inner hair-cell envelope extraction';
            
        end
        
     end
    
     % "Getter" methods
     methods
         function method = get.method(pObj)
             method = pObj.parameters.map('IHCMethod')
         end
     end
     
     
end