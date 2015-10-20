%calcASW.m
%
% asw = (dObj,varargin)
%
% This function predicts apparent source width (ASW) perception from ITD and ILDs
% !!!In future this function should be replaced by a the corresponding
% processor aswProc.m
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

function [asw, aswReprChan, itdReprChan, ildReprChan] = calcASW(dObj,varargin)

%% Inputs
% check inputs
if nargin<1||isempty(dObj)
    error('Please provide a data object!')
end
if isempty(findprop(dObj,'itd'))||isempty(findprop(dObj,'ild'))
    error('Please provide an object containing property ''itd'' and ''ild''!')
end

% optional parameters
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
if isfield(opt,'itdmax')
    itdmax = opt.itdmax;
else
    itdmax = [-1e-3 1e-3];
end
if isfield(opt,'ildmax')
    ildmax = opt.ildmax;
else
    ildmax = [-20 20];
end
if isfield(opt,'chanRepresent')
    chanRepresent = opt.chanRepresent;
else
    chanRepresent = 'width';
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

%% Processing
% get itd and ild data
itd = dObj.itd{1}.Data(:);
ild = dObj.ild{1}.Data(:);
cfHz = dObj.itd{1}.cfHz;
Nchan = size(dObj.itd{1}.Data(:),2);

% Width per channel
[itdWidthChan, itdLR_boarderChan, itdPrctChan] = calcDistrWidth(itd,wdMethod,percent); %calculate width of binCue directly for each frequency channel
[ildWidthChan, ildLR_boarderChan, ildPrctChan] = calcDistrWidth(ild,wdMethod,percent); %calculate width of binCue directly for each frequency channel

% choose representation of ASW per channel
switch chanRepresent
    case 'width' %width
        itdReprChan = itdWidthChan;
        ildReprChan = ildWidthChan;
        
    case 'lr' %left and right boundary seperately
        itdReprChan = itdLR_boarderChan;
        ildReprChan = ildLR_boarderChan;
        
    case 'prct' %use all percentiles
        if strcmp(chanRepresent,wdMethod)
            itdReprChan = itdPrctChan;
            ildReprChan = ildPrctChan;
        else
            warning(['This representation method is only valid for a '],...
                  ['percentile-based calculation! Using ''lr'' method instead!'])
            itdReprChan = itdLR_boarderChan;
            ildReprChan = ildLR_boarderChan;             
        end  
end

% Transformation of bin. cue
switch transMethod
    case 'none' %no transformation
        %do nothing
        
    case 'norm' %normalization
        itdReprChan = itdReprChan/itdmax(2);
        ildReprChan = ildReprChan/ildmax(2);
        
    case 'latcomp' %lateral compression
        itdReprChan = latCompression(itdReprChan,itdmax(2));
        ildReprChan = latCompression(ildReprChan,ildmax(2));
        
    case 'mapping' %Map boarders (percentiles/std) of ITDs and ILDs to azimuthal angle
        % init
        itd_mapped = zeros(2,Nchan);
        ild_mapped = zeros(2,Nchan);
        angleoffset = -91; %offset to find the correct angle [degree]
        itdload = load('ITD2Azimuth_Subband.mat');
        itd2azim = itdload.mapping.itd2azim;
        ildload = load('ILD2Azimuth_Subband.mat');
        ild2azim = ildload.mapping.ild2azim;
        
        for ii = 1:Nchan %all channels
            for jj = 1:size(itdReprChan,1) %all percentiles/stds
                % itd mapping
                if isempty(find(itdReprChan(jj,ii)<=itd2azim(:,ii))) %exceeds max/min values?
                    if itdReprChan(jj,ii)<=0
                        itd_mapped(jj,ii) = -90;
                    else
                        itd_mapped(jj,ii) = 90;
                    end
                else
                    itd_mapped(jj,ii) = find(itdReprChan(jj,ii)<=itd2azim(:,ii),1) + angleoffset;
                end
                
                % ild mapping
                if isempty(find(ildReprChan(jj,ii)<=ild2azim(:,ii))) %exceeds max/min values?
                    if ildReprChan(jj,ii)<=0
                        ild_mapped(jj,ii) = -90;
                    else
                        ild_mapped(jj,ii) = 90;
                    end
                else
                    ild_mapped(jj,ii) = find(ildReprChan(jj,ii)<=ild2azim(:,ii),1) + angleoffset;
                end
            end
        end
        
        % Overwrite channel representation per cue with results from
        % mapping
        itdReprChan = itd_mapped;
        ildReprChan = ild_mapped;

otherwise
    error(['The transformation method ' transMethod ' is not defined!'])
    
end

% Frequency weighting
switch freqWeighting
    case 'none' %no weighting, i.e. use all bands 
        itdBandSelect = 1:Nchan;
        ildBandSelect = 1:Nchan;
        
    case 'itdlow' %lowpass itds, i.e. use itds only at low freqeuncies
        itdBandSelect = 1:20; %lowpass filter cut-off [#frequency band]
        ildBandSelect = 1:Nchan;
        
    case 'itdE3' %lowpass itds, i.e. use itds only at low freqeuncies
        itdBandSelect = 7:22; %#frequency bands that correspond to the averaging 
                        %performed in IACC_E3, 
                        %i.e. octave bands at    0.5,   1   and 2 kHz
                        %corresponding to bands 7:12, 11:16 and 17:22,
                        %respectively
        ildBandSelect = 1:Nchan;
        
    otherwise
        error(['The frequency weighting ' freqWeighting ' is not defined!'])
end

% Combination of both cues
switch combMethod
    case 'itd' %use only itd
        aswReprChan = itdReprChan;
        asw = nanmean(itdReprChan(:,itdBandSelect),2); %average all channels
        
    case 'ild' %use only ild
        aswReprChan = ildReprChan;
        asw = nanmean(ildReprChan(:,ildBandSelect),2); %average all channels
        
    case 'duplex' %combine itd and ild according to duplex theory
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        % cross-over frequency in Hz
        fcrossHz = 1500; 
        
        % init
        aswReprChan = zeros(size(itdReprChan));
        
        % choose channels below fcrossHz from itds and above from ilds
        bchanitd = (cfHz <= fcrossHz);
        aswReprChan(:,bchanitd) = itdReprChan(:,bchanitd);
        aswReprChan(:,~bchanitd) = itdReprChan(:,~bchanitd);
        
        % average all channels
        asw = nanmean(aswReprChan,2);
        
    case 'dominant' %choose the dominant cue in each channel, i.e. the one with higher variance
        if strcmp(transMethod,'none')
            warning('Bin. cues need to be transformed before they can be combined!')
        end
        
        bitdDominance = (itdReprChan>ildReprChan);
        
        % init
        aswReprChan = zeros(size(itdReprChan));
        
        % choose channels for either dominance
        aswReprChan(:,bitdDominance) = itdReprChan(:,bitdDominance);
        aswReprChan(:,~bitdDominance) = itdReprChan(:,~bitdDominance);  
        
        % average all channels
        asw = nanmean(aswReprChan,2);
        
    otherwise
        error(['The combination method ' combMethod ' is not defined!'])

end

%EOF