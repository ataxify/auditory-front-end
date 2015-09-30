function binaural = spatializeAudio(audio,fs,azimuth,room)

% Check for proper input arguments
if nargin ~= 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Number of audio files
[nSamples,nSources] = size(audio);

% Allocate memory
binaural = zeros(nSamples,2);

% Loop over number of audio files
for ii = 1 : nSources
    % Normalize RMS
    sig = audio(:,ii)/sqrt(mean(audio(:,ii).^2));
    
    % Spatialize signal using HRTF processing
    binaural = binaural + filterHRTF(sig,fs,azimuth(ii),room);
end