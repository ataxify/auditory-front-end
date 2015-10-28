%latCompression.m
%Laterality compression according to Goupell et al. (2007) eq.(7) and (8)
%"Interaural fluctuations and the detection of interaural incoherence.
%III. Narrowband experiments and binaural models" 
%
%"A small static interaural difference leads to a small dis-
%placement in the lateral position of the auditory image from a 
%centered position. A greater interaural difference leads to a 
%greater displacement, but increasing interaural differences produce 
%diminishing returns because the laterality is a compressive function of 
%interaural differences. A perceptual model for fluctuations can easily adopt 
%this effect from static experiments. The compression functions, to be 
%called “laterality compression,” used in the present analysis were 
%exponential fits to the data from Yost's 1981 experiments."
%
%inputs:
% x - dynamic binaural cue
% xrange - maximal dynamic range of x
%
%output:
% xcomp - binaural cue after laterality compression
%
%implemented by Johannes Käsbach, CAHR, DTU, 14-07-2014

function xcomp = latCompression(x,xrange)

if nargin < 2
    error('Please provide the two inputs ''x'' and its dynamic range ''xrange''!')
end

xcomp = 10.*sign(x).*(1-exp(-abs(x)./xrange));

end