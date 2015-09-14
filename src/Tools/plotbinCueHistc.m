%plotbinCueHistc.m
%
% plotbinCueHistc(Nhist,gravityFunction,varargin)
%
%This function plots a histogram 'Nhist' calculated by function
%binCueHistc.m 
%Either a surf or image plot can be plotted specified by 'method'.
%
%inputs:
% Nhist - The histogram with dimensions [nbins x Nchan] calculated by binCueHistc.m 
%        The plot will be along dimensions [x y Ndist]
% optional:
% 'method' - 'IAMGE' uses imagesc, 'SURF' uses surf plot
% 'gravity' - bool: plot gravity of histogram(1)
% 'clim' - color scale limits for IMAGE plot (see imagesc.m)
% 'gravlim' - ylim for gravity function plot
% 'brel' - bool: relative occurance (normalise data) (1)
% 'ymax' - maximum value for y-vector
% 'labels' - xlabel and ylabel
% 'S' - plot settings: FontSize
%
%-----------------------------------------------
%Revision:
% 17-10-2013 by JK: including flexible inputs via varargin 
% 24-03-2014 by JK: - changed 'thresh' to 'actbands' (in order to omit a second estimation of the activated bands inside this function)
%                   - changed histogram-calculation from 'hist' to 'histc' 
%-----------------------------------------------
%
%by Johannes K?sbach (JK), DTU, CAHR, 08-October-2013
%
% see also: binCueHistc.m, surf.m, imagesc.m

function plotbinCueHistc(Nhist,gravityFunction,varargin)

%% Definitions (fixed Parameters)
gravplwidth = 0.1; %width of plot for gravity
shift = 0.3; %shift in width for imagesc and pos(1) for gravity plot
xtick = 5; %every 'xtick' tick on xaxis
% fsax = 14; %FontSize Axis -> substituted by S.fsax

%% check inputs
if nargin<2||isempty(gravityFunction) 
    gravity_flag = 0;
else
    gravity_flag = 1; %bool: plot gravityFunction(1)
end

r = struct(varargin{:});
try 
    clim = r.clim;
catch
    % clim = [0 1]; %limits for colorbar in imagesc
    clim = [0 0.3]; %limits for colorbar in imagesc for FA2014 presentation
end
try
    gravlim = r.gravlim;
catch
    gravlim = 0.025; %limit of x-axis (y-axis before rotation on view) of plot for gravity
end
try
    color_flag = r.color_flag;
catch
    color_flag = 0; %bool: plot colorbar(1) 
end
try
    method  = r.method;
catch
    method = 'IMAGE'; %imagesc plot
end
try
    brel  = r.brel;
catch
    brel = 1; %relative occurance (normalise data)
end

try
    ymax = r.ymax;
catch
    ymax = max(max(abs(Nhist))); %ylim according to data
end
try
    labels = r.labels;
catch
    labels = {'Filterbank Channel';'Binaural Cue'}; %labels
end
try
    xticklabel = r.f0;
catch
    %do nothing
end
try
    markers = r.markers;
catch
    markers = 0;
end
try
    S = r.S;
catch
    load S %common plot settings
end
% end

if length(labels)<3
    if(brel)
        labels{3} = 'Rel. occurance'; %(normalised data)
    else
        labels{3} = 'Abs. occurance';
    end
end

%% Plots

%number of containers 'nbins' and frequency channels 'Nchan'
nbins = size(Nhist,1);
Nchan = size(Nhist,2);

%x and y coordinate
x = 1:Nchan;
y = linspace(-ymax,ymax,nbins);

figure
if gravity_flag
    subplot(1,2,1)
end
switch method
    case 'SURF'
        [Xmesh, Ymesh] = meshgrid(x,y); %meshgrid
        surf(Xmesh,Ymesh,Nhist);
        zlabel(labels{3},'FontSize',S.fstxt)
        plotone = gca;
    case 'IMAGE'
        imagesc(x,y,Nhist,clim);
        colormap('bone')
        colormap(flipud(colormap))
        if color_flag
            colorbar
        end
        if markers
           hold on
           plot(markers(1),y,'Color','k','Linewidth',S.lw)
           plot(markers(2),y,'Color','k','Linewidth',S.lw)
           hold off
        end
        title(labels{3},'FontSize',S.fstxt)
        plotone = gca;
end

%labels
xlabel(labels{1},'FontSize',S.fstxt)
ylabel(labels{2},'FontSize',S.fstxt)
if ~exist('xticklabel','var')
    xticklabel = x;
end
set(plotone,'XTick',xtick:xtick:length(xticklabel))
set(plotone,'XTickLabel',xticklabel(xtick:xtick:end))

%axes
xlim([1 Nchan])
ylim([-ymax ymax])
if brel
    zlim(clim)
end

set(gca,'FontSize',S.fsax)

%gravity
if gravity_flag
    subplot(1,2,2)
    plot(y,gravityFunction,'LineWidth',S.lw,'Color','k')
    view(90,90)
%     title('Centre of gravity','FontSize',S.fstxt)
    plottwo = gca;
    pos1 = get(plotone,'Position');
    pos2 = get(plottwo,'Position');
    hpytick = get(plotone,'YTick');
    set(plotone,'Position',[pos1(1) pos1(2) pos1(3)+shift pos1(4)])
    set(plottwo,'Position',[pos2(1)+shift-0.1 pos2(2) gravplwidth pos2(4)])
    set(gca,'XTick',hpytick)
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',{'0';'';num2str(gravlim)})
    xlim([-ymax ymax])
    ylim([0 gravlim])
    set(gca,'FontSize',S.fsax)
end
        
end
%EOF