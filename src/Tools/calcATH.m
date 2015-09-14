%calcATH.m
%This function calculates the Absolute Threshold of Hearing (ATH) according to equation (1) in
%Terhardt, E. (1979). “Calculating virtual pitch,” Hear. Res. 1, 155–182.
%
%inputs:
% f: frequency [kHz]
%outputs:
% ATH: Absolute Threshold of Hearing [dB SPL]
%
%by Johannes Käsbach, DTU, CAHR, 02-10-2013

function ATH = calcATH(f)

ATH = 3.64*f.^(-0.8)-6.5.*exp(-0.6*(f-3.3).^2)+10^(-3)*f.^4;

end