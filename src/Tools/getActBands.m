%getActBands.m
%
% actbands = getActBands(timeFreq)
%
% This function outputs a column-vector of active frequency bands in a
% time-frequency (time x frequency) signal. It consideres all channels to
% be active that have non-zero components.
%
% Inputs:
% timeFreq - Time-frequency signal (works also for container-frequency signals)
% 
% Outputs:
% actbands - Column-vector of active frequency bands
% Nchanact - #active channels
%
% by johannes Käsbach, DTU, CAHR, johk@elektro.dtu.dk, 08. September 2015
%
%---------History-----------
% - no further editing done
%
%---------------------------

function [actbands, Nchanact] = getActBands(timeFreq)

%check inputs
if nargin<1||isempty(timeFreq)
    error('Please provide a time-frequency signal!')
end

%Number of frequency channels 'Nchan'
Nchan = size(timeFreq,2);

%init actbands
actbandsInit = 1:Nchan;

%find non-zero channels
bpassTh = sum(~isnan(timeFreq),1)>0;

%active bands 'actbands'
actbands = actbandsInit(bpassTh);

%#active channels 'Nchanact'
Nchanact = size(actbands,2);

end

%EOF

