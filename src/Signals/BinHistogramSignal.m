classdef BinHistogramSignal < Signal
%BINHISTOGRAM Signal class for two-dimensional, container-frequency representations (frequency-dependent histogram).
%   This children signal class regroups all signal that are some sort of container frequency 
%   representation.
%
%   BINHISTOGRAMSIGNAL properties:
%       cfHz - Center frequencies of audio frequency channels (Hz)
%
% See also Signal, gammatoneProc, drnlProc, ihcProc, icProc, ildProc, itdProc, onsetProc, 
% offsetProc
    
    properties (SetAccess=protected)
        cfHz        % Center frequencies of the frequency channels
        binMethod        % Type of binaural cue
    end
       
    properties (GetAccess = protected)
%         scaling
    end
    
    methods 
        function sObj = BinHistogramSignal(procHandle,bufferSize,channel,data)
            %BinHistogramSignal   Class constructor 
            %
            %USAGE
            %     sObj = BinHistogramSignal(procHandle)
            %     sObj = BinHistogramSignal(procHandle,data)
            %
            %INPUT ARGUMENTS
            % procHandle : Handle to the processor generating this signal as output
            % bufferSize : Size of the ring buffer in s (default: bufferSize = 10)
            %    channel : Flag indicating 'left', 'right', or 'mono' (default: 
            %              channel = 'mono')
            %       data : Array of amplitudes to construct an object from existing data
            %
            %OUTPUT ARGUMENT
            %  sObj : Container-frequency signal object inheriting the signal class
             
            if nargin<4; data = []; end
            if nargin<3||isempty(channel); channel = 'mono'; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1||isempty(procHandle); procHandle = emptyProc; end
            
            sObj = sObj@Signal( procHandle, bufferSize, ...
                                length(procHandle.getDependentParameter('fb_cfHz')));
            
            if nargin>0     % Safeguard for Matlab empty calls
            
            sObj.Dimensions = 'nContainers x nFilters';
            sObj.cfHz = procHandle.getDependentParameter('fb_cfHz');
            sObj.setData( data );
            sObj.Channel = channel;
            
%             % TODO: do something non-specific for the scaling?
%             if isa(procHandle,'ratemapProc')
%                 sObj.scaling = procHandle.scaling;
%             else
%                 sObj.scaling = 'magnitude';
%             end
            
            end
            
        end
        
        function h = plot(sObj,h0,p,varargin)
            %plot       This method plots the data from a container-frequency
            %           domain signal object
            %
            %USAGE
            %       sObj.plot
            %       sObj.plot(h_prev,p,...)
            %       h = sObj.plot(...)
            %
            %INPUT ARGUMENT
            %  h_prev : Handle to an already existing figure or subplot
            %           where the new plot should be placed
            %       p : Structure of non-default plot parameters (generated
            %           from genParStruct.m)
            %
            %OUTPUT ARGUMENT
            %       h : Handle to the newly created figure
            %
            %OPTIONAL ARGUMENTS
            % 'rangeSec' - Vector of time limits for the plot
            
            if ~isempty(sObj.Data)
            
                % Manage plotting parameters
                if nargin < 3 || isempty(p) 
                    % Get default plotting parameters
                    p = Parameters.getPlottingParameters('BinHistogramSignal');
                else
                    defaultPar = Parameters.getPlottingParameters('BinHistogramSignal');
                    defaultPar.replaceParameters(p);
                    p = defaultPar;
                end
                
                % Manage optional arguments
                % Read in optional arguments
                if nargin>3 && ~isempty(varargin)
                    opt = struct;
                    for ii = 1:2:size(varargin,2)
                        opt.(varargin{ii}) = varargin{ii+1};
                    end
                else
                    opt = [];
                end
                
                % Nhist
                Nhist = sObj.Data(:);
                
                % Fixed Parameters
                gravplwidth = 0.1; %width of plot for gravity
                shift = 0.3; %shift in width for imagesc and pos(1) for gravity plot
                
                %number of containers 'nbins' and frequency channels 'Nchan'
                nbins = size(Nhist,1);
                Nchan = size(Nhist,2);
                
                %consider only activated bands
                [~, Nchanact] = getActBands(Nhist);
                
                %ymax 
                switch p.binMethod
                    case 'itd'
                        ymax = 1e-3; %1 ms
                        labels = {'Frequency [Hz]';'ITD'};
                        
                    case 'ild'
                        ymax = 10; %10 dB
                        labels = {'Frequency [Hz]';'ILD'};
                        
                    case 'ic'
                        ymax = 1; %correlation
                        labels = {'Frequency [Hz]';'IC'};
                        
                    otherwise
                        error(['The function is not specified for signal ' p.binMethod '!'])
                        
                end
                
                % Optional Parameters
%                 if ~isempty(opt) && isfield(opt,'ymax')
%                     ymax = opt.ymax;
%                 else
%                     ymax = max(max(abs(Nhist))); %ylim according to data
%                 end
%                 if ~isempty(opt) && isfield(opt,'labels')
%                     labels = opt.labels;
%                 else
%                     labels = {'Frequency [Hz]';'Binaural Cue'}; %labels
%                 end
                if length(labels)<3
                    if(p.brel)
                        labels{3} = 'Rel. occurance'; %(normalised data)
                    else
                        labels{3} = 'Abs. occurance';
                    end
                end

                % Manage handles
                if nargin < 2 || isempty(h0)
                        h = figure;             % Generate a new figure
                    elseif get(h0,'parent')~=0
                        % Then it's a subplot
                        figure(get(h0,'parent')),subplot(h0)
                        h = h0;
                    else
                        figure(h0)
                        h = h0;
                end
                
                % Plot the figure
                %x and y coordinate
                x = 1:Nchan;
                y = linspace(-ymax,ymax,nbins);
                
                if p.bgravity
                    subplot(1,2,1)
                end
                switch p.plMethod
                    case 'SURF'
                        [Xmesh, Ymesh] = meshgrid(x,y); %meshgrid
                        surf(Xmesh,Ymesh,Nhist);
                        zlabel(labels{3},'fontsize',p.map('fsize_title'),'fontname',p.map('ftype'))
                        plotone = gca;
                    
                    case 'IMAGE'
                        imagesc(x,y,Nhist,p.clim);
                        % Set the color map
                        try
                            colormap(p.map('colormap'))
                        catch
                            warning('No colormap %s is available, using ''bone''.',p.map('colormap'))
                            colormap('bone')
                        end
                        colormap(flipud(colormap))
                        if p.bColorbar
                            colorbar
                        end
                        title(labels{3},'fontsize',p.map('fsize_title'),'fontname',p.map('ftype'))
                        plotone = gca;
                end
                
                %labels
                xlabel(labels{1},'fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                ylabel(labels{2},'fontsize',p.map('fsize_label'),'fontname',p.map('ftype'))
                
                %axes
                xlim([1 Nchan])
                ylim([-ymax ymax])
                if p.brel
                    zlim(p.clim)
                end
                
                % Managing frequency axis ticks for auditory filterbank
                %
                % Find position of y-axis ticks
                M = size(sObj.cfHz,2);  % Number of channels
                n_points = 500;         % Number of points in the interpolation
                interpolate_ticks = spline(1:M,sObj.cfHz,...
                    linspace(0.5,M+0.5,n_points));
                %
                % Restrain ticks to signal range (+/- a half channel)
                aud_ticks = p.map('aud_ticks');
                aud_ticks=aud_ticks(aud_ticks<=interpolate_ticks(end));
                aud_ticks=aud_ticks(aud_ticks>=interpolate_ticks(1));
                n_ticks = size(aud_ticks,2);        % Number of ticks
                ticks_pos = zeros(size(aud_ticks)); % Tick position
                %
                % Find index for each tick
                for ii = 1:n_ticks
                    jj = find(interpolate_ticks>=aud_ticks(ii),1);
                    ticks_pos(ii) = jj*M/n_points;
                end
                
                % Set up y-axis
                set(plotone,'XTick',ticks_pos,...
                    'XTickLabel',aud_ticks,'fontsize',p.map('fsize_axes'),...
                    'fontname',p.map('ftype'))
                
                %gravityFct
                if p.bgravity
                    %sum distributions of all frequency channels and normalise to all active channels
                    gravityFctAvFirst = nansum(Nhist,2)/Nchanact;
                    
                    %plot
                    subplot(1,2,2)
                    plot(y,gravityFctAvFirst,'LineWidth',p.map('linewidth_m'),'Color','k')
                    view(90,90)
                    plottwo = gca;
                    pos1 = get(plotone,'Position');
                    pos2 = get(plottwo,'Position');
                    hpytick = get(plotone,'YTick');
                    set(plotone,'Position',[pos1(1) pos1(2) pos1(3)+shift pos1(4)])
                    set(plottwo,'Position',[pos2(1)+shift-0.08 pos2(2) gravplwidth pos2(4)])
                    set(gca,'XTick',hpytick)
                    set(gca,'XTickLabel',[])
                    set(gca,'YTickLabel',{'0';'';num2str(p.gravlim)})
                    xlim([-ymax ymax])
                    ylim([0 p.gravlim])
                    set(gca,'fontsize',p.map('fsize_axes'),...
                    'fontname',p.map('ftype'))
                end
                
            else
                warning('This is an empty signal, cannot be plotted')
            end
                
            
        end
    end
    
    methods (Static)
       
        function sObj = construct(fs,bufferSize,name,label,cfHz,channel,data)
            %construct
            %
            %
            
            if nargin<7; data = []; end
            if nargin<6; channel = []; end
            if nargin<5; cfHz = []; end
            if nargin<4; label = []; end
            if nargin<3; name = []; end
            if nargin<2||isempty(bufferSize); bufferSize = 10; end
            if nargin<1; fs = []; end
            
            % Create a dummy structure with that information to emulate a processor and
            % correctly call the class constructor
            dummyStruct = struct;
            dummyStruct.FsHzOut = fs;
            dummyStruct.getProcessorInfo.requestName = name;
            dummyStruct.getProcessorInfo.requestLabel = label;
            dummyStruct.getDependentParameter = containers.Map;
            dummyStruct.getDependentParameter('fb_cfHz') = cfHz;
            
            % Instantiate the signal
            sObj = BinHistogramSignal(dummyStruct,bufferSize,channel,data);
            
        end
        
        function [names, defaultValues, descriptions] = getPlottingParameterInfo()
            %GETPLOTTINGPARAMETERINFO   Stores plot parameters that are common to all
            %signals.
            
            
            names = {'plMethod',...
                    'binMethod',...
                    'colormap',...
                    'bColorbar',...
                    'clim',...
                    'aud_ticks',...
                    'brel',...
                    'bgravity',...
                    'gravlim',...
                    };
                 
            descriptions = {'Plotting Method (Choose: ''IMAGE'' or ''SURF'')',...
                    'Choose which binaural cue to plot ''itd'', ''ild'' or ''ic''',...
                    'Colormap for binaural histogram plots',...
                    'Boolean for displaying colorbar in binaural histograms plots',...
                    'Dynamic range for time-frequency plots (dB)',...
                    'Auditory ticks for ERB-based representations',...
                    'Boolean for displaying relative occurance (normalise data to number of frames)',...
                    'Boolean for displaying gravity/centroid function next to histogram',...
                    'xlim for gravity plot (ylim before rotation on view)',...
                    };
                
            defaultValues = {'IMAGE',...
                    'itd',...
                    'bone',...
                    0,...
                    [0 0.3],... %limits for colorbar in imagesc for FA2014 presentation %otherwise: [0 1];
                    [100 250 500 1000 2000 4000 8000 16000 32000],...
                    true,...
                    true,...
                    0.2,...
                    };
            
        end
        
        
    end
end
