function [SIGNALS,FEATURES,STATES] = process_WP2(earSignals,STATES)
%process_WP2   Perform WP2 processing
%
%USAGE
%      [SIGNALS,FEATURES,STATES] = process_WP2(binaural,STATES)
%
%INPUT PARAMETERS
%     binaural : binaural signal [nSamples x 2]
%       STATES : settings
% 
%OUTPUT PARAMETERS
%      SIGNALS : extracted signals (e.g. output of gammatone filterbank)
%     FEATURES : extracted features (e.g. estimated ITD / azimuth)
%       STATES : settings

%   Developed with Matlab 8.2.0.701 (R2013b). Please send bug reports to:
%   
%   Author  :  Tobias May, � 2013
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
% 
%   History :  
%   v.0.1   2014/01/31
%   ***********************************************************************

%% CHECK INPUT ARGUMENTS 
% 
% 
% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Short cut
fsHz = STATES.signal.fsHz;



%% CREATE PERIPHERAL AUDITORY SIGNAL
% 
% 
% Compute peripheral auditory signal
auditorySignals = auditoryPeriphery(earSignals,fsHz,STATES.periphery);


%% EXTRACT FEATURES
% 
% Do not assume any prior knowledge, extract all potential positions
nSources = inf;

% Estimate sound source azimuth
[azim,salience] = estimate_Azimuth(auditorySignals,fsHz,STATES.binaural.winSizeSec,nSources);
        


%% CREATE OUTPUT 
% 
% 
% Signal struct
SIGNALS.fsHz = STATES.signal.fsHz;
SIGNALS.data = auditorySignals;

% Feature struct
FEATURES.azimuth.direction = azim;
FEATURES.azimuth.salience  = salience;

% STATES struct
% STATES = STATES;


