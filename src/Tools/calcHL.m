%calcHL.m
%
% ------------------------ Modelling ASW -------------------------------- %
%
% This framework concerns modelling apparent source width (ASW) perception
% based on data measured and presented in Käsbach et al. (2013-2015)
%
% by Johannes Käsbach (JK), Centre for Applied Hearing Research (CAHR),
% Technical University of Denmark (DTU), johk@elektro.dtu.dk
%
% ----------------------------------------------------------------------- %
%
% This processor estimates the Hearing Level (HL) either based on 
% the absolute threshold of hearing  ('ath') or on the energy in each filter ('en').
% The thresholds are computed for each filter of the auditory filterbank.
% Implementation as in hlProc.m!
%
% inputs:
% in        - Left and right input signal [time x frequency]
% method    - Method used to estimate the thresholds 
%             absolute threshold of hearing ('ath') or energy-based ('en')
% enTh      - Threshold for energy-based method
%
% output:
% out       - Left and right input signal with applied threshold
% totEn     - SPL in each frequency-band ('ath')
%             Total energy in each frequency-band ('en')
% actbands  - Vector containing the frequency-bands that are above
%             threshold
%
% --------History-----------
% - 30. September 2015: Modelling ASW framework
% - 23. October 2015
% --------------------------

function [out_l, out_r, totEn ,actbands, weights] = calcHL(in_l,in_r,method,cfHz,FsHzIn,enTh,hl_wname,hl_wSizeSec,hl_hSizeSec)

%% Check Inputs
if nargin<2
    error('Please provide left and right input!')
end
if nargin<3 || isempty(method)
    method = 'none';
end
if nargin<4 || isempty(FsHzIn)
    error('Please provide a center-frequency vector!')
end
if nargin<5 || isempty(FsHzIn)
    FsHzIn = 44.1e3; %Sampling Frequency [Hz]
end
if nargin<6 || isempty(enTh)
    enTh = 0.5; %Threshold for energy-estimation
end
if nargin<7 || isempty(hl_wname)
    hl_wname = 'hann'; %Window name
end
if nargin<8|| isempty(hl_wSizeSec)
    hl_wSizeSec = 20E-3; %Window duration (s)
end
if nargin<9|| isempty(hl_hSizeSec)
    hl_hSizeSec = 10E-3; %Window step size (s)
end

% frame-based estimation?
bframebased = 1;

% get dimensions of input
[nSamples,nChannels] = size(in_l);

% Window
wSize = 2*round(hl_wSizeSec*FsHzIn/2);
hSize = round(hl_hSizeSec*FsHzIn);
win = window(hl_wname,wSize);

% How many frames are in the buffered input?
nFrames = floor((nSamples-(wSize-hSize))/hSize);

% Compute HL:

% init SPL calculations
frame_l_SPL = zeros(1,nChannels);
frame_r_SPL = zeros(1,nChannels);

% Pre-allocate outputs
out_l = NaN(nSamples,nChannels);
out_r = NaN(nSamples,nChannels);

% Window processing
winRep = repmat(win,1,nChannels);

% Estimate SPL
xref = 20e-6; %ref for dB SPL definition: xrms = 1 -> Lp = 94dB SPL
%(xrms = xref = 20e-6 -> Lp = 0dB SPL)
% xref = 1e-5; %ref for digital value definition: xrms = 1 -> Lp = 100dB
%              %using this reference will reproduce totEnergy = Lp
%              (Lp being the specified level of the input signal)

if bframebased
    % Loop on the time frame
    for ii = 1:nFrames
        % Get start and end indexes for the current frame
        n_start = (ii-1)*hSize+1;
        n_end = (ii-1)*hSize+wSize;
        
        % Energy in the windowed frame for left and right input
        frame_l_rms = sqrt(mean(power(winRep.*in_l(n_start:n_end,:),2)));
        frame_r_rms = sqrt(mean(power(winRep.*in_r(n_start:n_end,:),2)));
        
        %update SPL and devide by number of frames for averaging
        %between two successive frames
        frame_l_SPL = 20*log10(0.5*(frame_l_rms/xref + 10.^(frame_l_SPL/20))); %acc. SPL left channel
        frame_r_SPL = 20*log10(0.5*(frame_r_rms/xref + 10.^(frame_r_SPL/20))); %acc. SPL right channel
        
    end
    
else
    % Energy in the windowed frame for left and right input
    frame_l_rms = sqrt(mean(power(in_l,2)));
    frame_r_rms = sqrt(mean(power(in_r,2)));
    
    frame_l_SPL = 20*log10(frame_l_rms/xref); %acc. SPL left channel
    frame_r_SPL = 20*log10(frame_r_rms/xref); %acc. SPL right channel
end

% select maximum values from left and right channel
max_lr_SPL = max(frame_l_SPL,frame_r_SPL);

% frame_SPL contains total SPL after frame-based processing
% (check if it corresponds to the specified level)
totEn = 20*log10(sum(10.^(max_lr_SPL/20))); %tot. level from all filters

% Weights normalized to frequency band with highest SPL
weights = max_lr_SPL/max(max_lr_SPL);

% calculate active frequency bands
switch method
    case 'none'
        actbands = (1:length(cfHz)); %column-vector with #filters
        
    case 'ath'
        ath = calcATH(cfHz*1e-3); %absolute threshold of hearing

        %all filters that are above absolute threshold of hearing
        bpassth = max_lr_SPL>ath; %bool
        
        %active bands
        actbands = (1:length(cfHz)); %column-vector with #filters
        actbands = actbands(bpassth); %select activated channels from the vector
        
    case 'en'
        warning('Energy-based method not implemented yet!')
        
    otherwise
        error(['The method ' pObj.Method ' is not defined!'])
end

% Set all frames to 0 that are below the specified threshold
out_l(:,actbands) = in_l(:,actbands);
out_r(:,actbands) = in_r(:,actbands);

%EOF