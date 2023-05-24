
% Gu lab, IOBS, Fudan University.
% Written by Zhang Yimu in July, 2022.
%% 
clear;
plxfiledir = 'Z:\zym\IS\correlation\random_total\2_chr2_3\';
exportdir = 'Z:\zym\IS\correlation\random_total\2_chr2_3\';

subdirpath = fullfile(plxfiledir, '*.plx');

dats = dir(subdirpath);
for ii = 1:length(dats)
    plxfilename = dats(ii).name;
    plxdata = readPLXFileC(fullfile(plxfiledir,plxfilename),'all');
    disp(fullfile(plxfiledir,plxfilename));
    plxfilename = plxfilename(1:end-4);
    if ~exist([plxfiledir, plxfilename], 'dir')
        mkdir([plxfiledir, plxfilename])
    
    end
    %% 
    ch_num = [1,1,1,1,5,5,5,5,9,9,9,9,13,13,13,13];
    unit_name = 'abcde';
    channel = [];
    ch_idx = [];
    for i = 1:16
        Nunit = getfield(plxdata.SpikeChannels,{i,1},'NUnits');
        channel_idx = getfield(plxdata.SpikeChannels,{i,1},'Channel');
        if channel_idx > 16
            channel_idx = mod(channel_idx,16);
            if channel_idx == 0
                channel_idx = 16;
            end
        end
        if Nunit > 0
            for unit_num = 1:Nunit
                chname = ['SPKC',num2str(ch_num(i),'%02d'),unit_name(unit_num)];
                channel = [channel,{chname}];
                ch_idx = [ch_idx,{channel_idx}];
            end
        end
    end
    channel_unique = unique(channel);
    
    for i = 1:16
        wavelist{i,1} = getfield(plxdata.SpikeChannels,{i,1},'Waves');
        unitlist{i,1} = getfield(plxdata.SpikeChannels,{i,1},'Units'); 
        tslist{i,1} = getfield(plxdata.SpikeChannels,{i,1},'Timestamps');
    end
    wave = struct('waves',wavelist, 'units', unitlist,'timestamps',tslist);
    
    index_cal = [];
    for ch = 1:length(channel)
        channel_name = char(channel(ch));
        chidx = ch_idx{ch};
        chunit = strfind(unit_name, channel_name(7));
        unit_num = double(getfield(wave,{chidx,1},'units'));
        unit_num1 =  repmat(unit_num', 64, 1);
        waveform = double(getfield(wave,{chidx,1},'waves'));
        spk_ts = double(getfield(wave,{chidx,1},'timestamps'));
        waveform1 = waveform.*(unit_num1==chunit);
        mean_waveform = mean(waveform1,2);
        x = 1:1:64;
        x1 = 1:0.01:64;
        y1 = interp1(x,mean_waveform,x1,'spline'); 
        [trough,trough_index] = min(mean_waveform);
        [peak,peak_index] = max(mean_waveform(trough_index:end));
        peak_index = peak_index + trough_index - 1;
        half_peak = 0.5* peak;
        [hppl,hpl]=min(abs(y1(1:peak_index)-half_peak)); 
        hpl_x = x1(hpl);
        [hppr,hpr]=min(abs(y1(peak_index:end)-half_peak));
        hpr_x = x1(hpr+peak_index-1); 
        half_peak_width = abs(hpr_x - hpl_x)*(1.6/64); 
        
        peak_trough_distance = (peak_index-trough_index)*(1.6/64);
        
        spk_ts1 = spk_ts.*(unit_num==chunit)/40000;  %ms
        spk_ts1(all(spk_ts1==0,2),:)=[];
        
        output(ch).unit = channel_name;
        output(ch).channel = cell2mat(ch_idx(ch));
        output(ch).peak_trough_distance = peak_trough_distance;
        output(ch).half_peak_width = half_peak_width;
        output(ch).amplitude = peak;
        output(ch).timestamp = spk_ts1;
        output(ch).mean_waveform = mean_waveform;
    end
    %% 
    max_idx = [];
    for ch_u = 1:length(channel_unique)
        amplitude = [];
        index = find(strcmp(channel,channel_unique(ch_u)));
        for i = 1:length(index)
            amplitude = [amplitude, output(index(i)).amplitude];
        end
        [amp, idx] = max(amplitude);
        max_idx = [max_idx, index(idx)];
    end
    
    figure;
    timestamp = [];
    for i = 1:length(max_idx)
        name = char(channel(max_idx(i)));
        output1(i).unit = name;
        output1(i).channel = cell2mat(ch_idx(max_idx(i)));
        output1(i).peak_trough_distance = output(max_idx(i)).peak_trough_distance;
        output1(i).half_peak_width = output(max_idx(i)).half_peak_width;
        output1(i).amplitude = output(max_idx(i)).amplitude;
        for j = 1:length(max_idx)
            timestamp = setfield(timestamp, name, output(max_idx(i)).timestamp);
        end
        mean_waveform_p = output(max_idx(i)).mean_waveform;
        
        other_unit_idx = find(strcmp(channel,output(max_idx(i)).unit));
        total_mean_waveform = [];
        for other = other_unit_idx
            other_waveform = output(other).mean_waveform;
            total_mean_waveform = [total_mean_waveform, other_waveform];
            x = xlswrite(fullfile([plxfiledir,plxfilename],[output(other).unit, '_total_waveform_data.xlsx']), total_mean_waveform);
        end

        plot(mean_waveform_p);
        title(output(max_idx(i)).unit);
        ylabel('10*Voltage(μV)');
        xticks(0:8:64);
        xticklabels(0:0.2:1.6);
        text(40,min(mean_waveform_p)+10, ...
            ['trough to peak=' num2str(output(max_idx(i)).peak_trough_distance) 'ms' newline ...
            '2nd peak FWHM=' num2str(output(max_idx(i)).half_peak_width) 'ms']);
        saveas(gcf,fullfile([plxfiledir, plxfilename], [output(max_idx(i)).unit, '_waveform.png']));
        s = xlswrite(fullfile([plxfiledir,plxfilename],[output(max_idx(i)).unit, '_max_waveform_data.xlsx']), mean_waveform_p);
        close;
    end
    
    start_kbd = double((plxdata.EventChannels(1).Timestamps)/40000);
    stop_kbd = double((plxdata.EventChannels(2).Timestamps)/40000);
    
    timestamp = setfield(timestamp, 'start', start_kbd);
    timestamp = setfield(timestamp, 'stop', stop_kbd);
    if ~isempty(plxdata.EventChannels(11).Timestamps)
        kbd1 = double(plxdata.EventChannels(11).Timestamps)/40000;
        timestamp = setfield(timestamp, 'kbd1', kbd1);
    end
    wave_table= struct2table(output1);
    writetable(wave_table, fullfile([plxfiledir,plxfilename,'waveform.xlsx']));
    save([plxfiledir, plxfilename,'timestamp.mat'],'timestamp','-v7.3');
    disp('export successfully!');
    %clearvars output output1 timestamp
end