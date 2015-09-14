%calcDistrWidth.m
%
% width = calcDistrWidth(distr,wdMethod,percent)
%
%This function calculated the width of a distribution. For a matrix the
%width will be calculated for each row
%
% inputs:
% distr: distribution to be analysed [time x channels]
% wdMethod: 'prct' using the difference of percentiles
%         'std' using the standard deviation with normalisation N-1 (default)
% percent: choose percent values for prctile calculation
%
% output:
% width: the calculated width of the distribution
% 
% ------------------History-------------------
% by Johannes Käsbach, DTU, CAHR, 02-10-2013
% updated 28. August 2015
% --------------------------------------------

function width = calcDistrWidth(distr,wdMethod,percent)

%% Check inputs
if nargin<3||isempty(percent)
    percent = [10 90]; %choose percent values for prctile calculation
end

if nargin<2||isempty(wdMethod)
    wdMethod = 'prct'; %choose 'method' for calculating width
end

%% Implementation
DIM = 1; %calculation along dimension DIM
switch wdMethod
    case 'prct'
        prct = prctile(distr,percent,DIM); %percentiles
        width = abs(prct(end,:) - prct(1,:)); %absolute difference between 90th and 10th percentile
        
    case 'std'
        width = std(distr,0,DIM);
end

end
