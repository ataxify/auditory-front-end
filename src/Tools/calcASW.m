%calcASW.m
%
% asw = (dObj,varargin)
%
%This function predicts apparent source width (ASW) perception from ITD and ILDs
%!!!In future this function should be replaced by a the corresponding
%processor aswProc.m
%
% Inputs:
% dObj - data object containing ITDs and ILDs
% 
% optional inputs (varargin):
% wdMethod - Method for estimating the distribution width ('prct' or 'std')
% percent  - Choose percent values for prctile calculation
% avMethod - Method for averaging across frequency channels
% transform - 'Choose transformation of bin. cue: 'none', 'norm' (normalize)  or 'latcomp' (lateral compression)'
% freqWeighting - Select frequency-dependent weighting of the binaural cue
%
% Outputs:
% asw - Prediction of ASW perception.
%
% by Johanens Käsbach, CAHR, DTU, johk@elektro.dtu.dk, 10. September 2015
%
%---------History-----------
% - no further editing done
%
%---------------------------

function asw = calcASW(dObj,varargin)

%% Inputs
%check inputs
if nargin<1||isempty(dObj)
    error('Please provide an data object!')
end
if isempty(findprop(dObj,'itd'))
    error('Please provide an object containing property ''itd''!')
end
if isempty(findprop(dObj,'ild'))
    error('Please provide an object containing property ''ild''!')
end

%optional parameters
opt = struct(varargin{:});
                             
if isfield(opt,'wdMethod')
    wdMethod = opt.wdMethod;
else
    wdMethod = 'prct';
end
if isfield(opt,'percent')
    percent = opt.percent;
else
    percent = [10 90];
end
if isfield(opt,'transMethod')
    transMethod = opt.transMethod;
else
    transMethod = 'none';
end
if isfield(opt,'freqWeighting')
    freqWeighting = opt.freqWeighting;
else
    freqWeighting = 'none';
end
if isfield(opt,'combMethod')
    combMethod = opt.combMethod;
else
    combMethod = 'itd';
end

%% Definitions
range.itd = [-1e-3 1e-3];
range.ild = [-20 20];

%% Processing
%get itd and ild data
itd = dObj.itd{1}.Data(:);
ild = dObj.ild{1}.Data(:);
cfHz = dObj.itd{1}.cfHz;
Nchan = size(dObj.itd{1}.Data(:),2);

%A) calculate width per channel first, then average across channels
%Width per channel
itdWidthChan = calcDistrWidth(itd,wdMethod,percent); %calculate width of binCue directly for each frequency channel
ildWidthChan = calcDistrWidth(ild,wdMethod,percent); %calculate width of binCue directly for each frequency channel

%Transformation of bin. cue
switch transMethod
    case 'none' %no transformation
        %do nothing
        
    case 'norm' %normalization
        itdWidthChan = itdWidthChan/range.itd(2);
        ildWidthChan = ildWidthChan/range.ild(2);
        
    case 'latcomp' %lateral compression
        itdWidthChan = latCompression(itdWidthChan,range.itd(2));
        ildWidthChan = latCompression(ildWidthChan,range.ild(2));
        
    otherwise
        error(['The transformation method ' transMethod ' is not defined!'])
        
end

%Frequency weighting
switch freqWeighting
    case 'none'
        
    otherwise
        error(['The frequency weighting ' freqWeighting ' is not defined!'])
end

%Comination of both cues
switch combMethod
    case 'itd' %use only itd
        asw = mean(itdWidthChan,2); %average all channels
        
    case 'ild' %use only ild
        asw = mean(ildWidthChan,2); %average all channels
        
    case 'duplex' %combine itd and ild according to duplex theory
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be comined!')
        end
        
        %cross-over frequency in Hz
        fcrossHz = 1500; 
        
        %init
        aswWidthChan = zeros(size(itdWidthChan));
        
        %choose channels below fcrossHz from itds and above from ilds
        bchanitd = (cfHz <= fcrossHz);
        aswWidthChan(:,bchanitd) = itdWidthChan(:,bchanitd);
        aswWidthChan(:,~bchanitd) = itdWidthChan(:,~bchanitd);
        
        %average all channels
        asw = mean(aswWidthChan,2);
        
    case 'dominant' %choose the dominant cue in each channel, i.e. the one with higher variance
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be comined!')
        end
        
        bitdDominance = (itdWidthChan>ildWidthChan);
        
        %init
        aswWidthChan = zeros(size(itdWidthChan));
        
        %choose channels for either dominance
        aswWidthChan(:,bitdDominance) = itdWidthChan(:,bitdDominance);
        aswWidthChan(:,~bitdDominance) = itdWidthChan(:,~bitdDominance);  
        
        %average all channels
        asw = mean(aswWidthChan,2);
        
    otherwise
        error(['The transformation method ' transMethod ' is not defined!'])

end

%EOF