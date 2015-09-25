%binCueHistc.m
%
% [Nhist, gravityFctAvFirst, gravityAvFirst, gravityFctWdthFirst, gravityWdthFirst] = binCueHistc(binCue,varargin)
%
%This function calculates for a binaural cue 'binCue' (itd or ild) with dimensions [frames x channels] 
%the histogram 'Nhist' in respect to 'nbins' containers using the function histc.m.
%The histogram can be plotted using function plotBinCueHistc.m
%The center of gravity of the histogram 'gravityFct' is the accumulated distribution
%across frequency channels
%The variable 'gravity' contains the width of the 'gravityFct' 
%calculated using function calcDistrWidth.m
%The width of the original distribution 'binCue' is the mean of the width
%per channel (also calculated using function calcDistrWidth.m)
%
%
% inputs:
% 'binCue': binaural cue matrix with dimensions [frames x channels]
% 
% optional inputs using 'varargin':
% 'nbins' - %#equally spaced containers in histogram calculation (default: nbins = 41)
% 'wdMethod' - 'prct' using the difference of percentiles (default)
%              'std' using the standard deviation with normalisation N-1
% 'percent' - specify the percent values P used in percentile calculation based on prctile.m 
% 'brel' - bool: relative occurance (normalise histogram to #frames) (1)
% 'actbands' - vector that specifies the activated bands/channels 
%              (bands under a certain threshold are omitted)
% 'ymax' - Threshold for vector 'y' specifying containers (EDGES) in histc.m
%
% output:
% Nhist: histogram of the distribution 'binCue'
% gravityFct: center of gravity, i.e. 
%             - the accumulated distribution across frequency channels 
%               (gravityFctAvFirst)
%             - width of distribution per channel (gravityFctWdthFirst)
%               (uses calcDistrWidth.m)
% gravity:    - width for the respective gravityFct
%               calculated using calcDistrWidth.m (gravityAvFirst)
%             - width calculating mean of width per channel
%               (gravityWdthFirst)
% ----------------------------------------------------------------------
% History:
% by Johannes Käsbach (JK), DTU, CAHR, 22-05-2014, adapted from plotHist3D.m
% updated 26. August 2015, by JK
% ----------------------------------------------------------------------
%
% see also: histc.m, prctile.m, std.m, plotBinCueHistc.m

function [Nhist, gravityFctAvFirst, gravityAvFirst, gravityFctWdthFirst, gravityWdthFirst] = binCueHistc(binCue,varargin)

%% check inputs
opt = struct(varargin{:});
try
    nbins = opt.nbins; %#bins in histogram calculation
catch
    nbins = 41;
end
try
    wdMethod = opt.wdMethod; %choose 'method' for calculating width in calcDistrWidth
catch
    wdMethod = 'prct';
end
try
    percent = opt.percent; %choose percent values for prctile calculation
catch
    percent = [10 90];
end
try
    brel  = opt.brel;
catch
    brel = 1; %relative occurance (normalise data)
end
try
    actbands = opt.actbands;
catch
    actbands = 1:size(binCue,2); %considering only activated bands
end
try
    ymax = opt.ymax;
    if isempty(ymax)
        ymax = max(max(abs(binCue)));
    end
catch
    ymax = max(max(abs(binCue))); %ylim according to data
end

%Implementation
%% Nhist
%Dimensions of the binCue: Nframes x Nchan
Nframes = size(binCue,1);
Nchan = size(binCue,2);

%vector for histc.m
y = linspace(-ymax,ymax,nbins);

%absolute maxium of the distribution
binCueMax = max(max(abs(binCue)));
disp('- binCueHistc.m -')
disp(['Maximum value of distribution: ' num2str(binCueMax)])
disp(['Specified allowed maximum value for histogram: ' num2str(ymax)])
if ymax < binCueMax
    disp(['Histogram is truncated to ''ymax'' = ' num2str(ymax)])
end
    
%init 
Nhist = NaN(nbins,Nchan);

%consider only activated bands
Nchanact = size(actbands,2); %number of activated channels, used to normalise gravity plot
binCue = binCue(:,actbands); %consider just activated bands

%Use histc to put distribution into containers specified in y
Nhist(:,actbands) = histc(binCue,y);

if(brel) %relative occurance: normalise to Nframes (each channel contains counts of all frames)
    Nhist = Nhist/Nframes;
end


%% gravity
gravityFctAvFirst = sum(Nhist,2)/Nchanact; %sum distributions of all frequency channels and normalise to all active channels
gravityAvFirst = calcDistrWidth(gravityFctAvFirst,wdMethod,percent); %calculate width of gravity function

gravityFctWdthFirst = calcDistrWidth(Nhist,wdMethod,percent); %calculate width of each channel of histogram
gravityWdthFirst = mean(gravityFctWdthFirst,2); %total width as averaging across all channels

end

%EOF