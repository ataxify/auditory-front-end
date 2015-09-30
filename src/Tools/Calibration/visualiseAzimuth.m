function h = visualiseAzimuth(dObj,realAz,h0)
%VISUALISEAZIMUTH   This functions displays a plot of the cross-correlation together with
%a distribution of azimuth, and optionally the estimated and actual source positions
%
% USAGE:
%   h = visualiseAzimuth(dObj)
%   h = visualiseAzimuth(dObj,estAz,realAz,h0)
%
% INPUT ARGUMENTS:
%    dObj : AFE data object
% peakVal : Vector of peak cross-correlation values at the estimated source azimuths
%  realAz : Vector of actual source azimuths in degrees
%      h0 : Handle to an existing figure
%
% OUTPUT ARGUMENT:
%       h : Handle to the new figure

% Checking that the provided data actually contains azimuth distribution
if ~isprop(dObj,'azimuthDistribution')
    error('No azimuth distributions were computed for that data object')
end

% Manage figure handle
if nargin < 3 || isempty(h0)
    h = figure;             % Generate a new figure
elseif get(h0,'parent')~=0
    % Then it's a subplot
    figure(get(h0,'parent')),subplot(h0)
    h = h0;
    clf(h)
else
    figure(h0)
    h = h0;
    clf(h)
end

if nargin < 2; realAz = []; end

% TODO: Generate multiple figures if there are several cross-correlation instances in the
% data object

% Plot the cross-correlation summary
h1 = subplot(3,1,1:2);
dObj.azimuthDistribution{1}.plot(h1);

h1Pos = get(h1,'position');

CCF_Warped_Sum = squeeze(mean(mean(dObj.azimuthDistribution{1}.Data(:),1),2));

% Plot the azimuth distribution
h2 = subplot(3,1,3); hold on;
h21 = plot(dObj.azimuthDistribution{1}.azimuth,CCF_Warped_Sum,'-');
set(h21,'color',[0.45 0.45 0.45],'linewidth',1.75)
h2Pos = get(h2,'position');
h2Pos(3) = h1Pos(3);
set(h2,'position',h2Pos)

xlim([-90 90])
ylim([0.95*min(CCF_Warped_Sum) 1.25*max(CCF_Warped_Sum)])
xlabel('Azimuth (deg)')

if isprop(dObj,'localization')
    % Get cross-correlation value at the peaks
    peakVal = interp1(dObj.azimuthDistribution{1}.azimuth,...
                      CCF_Warped_Sum,...
                      dObj.localization{1}.Data(end));

    % Plot estimated source position
    h22 = plot(dObj.localization{1}.Data(end),peakVal,'kx');
    set(h22,'MarkerSize',12,'LineWidth',2.5);
  
end    
   
if ~isempty(realAz)
    % Plot real source location
    realAz = realAz(:).';
    h23 = plot([realAz; realAz],[-1 1],':');
    set(h23,'color',[0.65 0.65 0.65],'linewidth',4)
    
end

hold off

end