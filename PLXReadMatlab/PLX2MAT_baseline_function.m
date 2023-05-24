function PLX2MAT_baseline_function(plxfiledir, plxfilename)
%% read in data
exportdir = plxfiledir;
plxdata = readPLXFileC([plxfiledir plxfilename '.plx'], 'all');
%% prepare data set
filename = plxfilename;
wavefreq = plxdata.ContinuousChannels(1).ADFrequency;
eventfreq = plxdata.WaveformFreq;

values = zeros(length(plxdata.ContinuousChannels(1).Values), length(plxdata.ContinuousChannels), 'int16');
for chidx = 1:length(plxdata.ContinuousChannels)
    values(:, chidx) = plxdata.ContinuousChannels(chidx).Values;
end

%marker_ts = [];
%marker_name = [];
%for chidx = 1:plxdata.NumEventChannels
%    marker_ts = [marker_ts, plxdata.EventChannels(chidx).Timestamps'];
%    if length(plxdata.EventChannels(chidx).Timestamps) > 0
%        marker_name = [marker_name, plxdata.EventChannels(chidx).Channel];
%    end
%end
%marker = [marker_name; marker_ts];
%marker = [marker_ts];

unit = plxdata.ContMaxMagnitudeMV / (2 ^ (plxdata.BitsPerContSample - 1) * ...
    plxdata.ContinuousChannels(1).ADGain * plxdata.ContinuousChannels(1).PreAmpGain);

%% export
%save([exportdir filename '.mat'], 'wavefreq', 'eventfreq', 'values', 'marker', 'unit', '-v7.3');
save([exportdir filename '.mat'], 'wavefreq', 'eventfreq', 'values', 'unit', '-v7.3');
disp('exported!')
