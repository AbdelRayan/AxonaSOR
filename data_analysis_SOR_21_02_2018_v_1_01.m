% The aim of this script is to read the LFP data provided by Andrew Lehr
% during object recognition task 
% Abdelrahman Rayan, abdelrahman.rayan@rub.de 
%% Loading the trial files 
clear all, close all, clc
% Choose the folder to be analyzed 
data_path = uigetdir;
cd(data_path)
% Choose the intertrial files
files            = dir( fullfile(data_path,'*.egf') );
TXTfiles         = dir( fullfile(data_path,'*.txt') );
TXTfilesName     = {TXTfiles.name}'; % extract the nae of the txt files 
flnmroot         = {files.name}'; % extract the names of the files from the folder
[flnmroot,INDEX] = sort_nat(flnmroot);
indexDay         = strfind(flnmroot{1},'_');
AnimalCode       = flnmroot{1}(1:indexDay(1)+1);

%% Reading the trial LFP data into one big data structure
LFP_data = [];
eegNum  = 1;
eegNum_2 = 2;

for iFile = 1:numel(flnmroot)
    Struct_name             = flnmroot{iFile}(1:end-4);                        % file name
    % Find '-' and replace them with '_'
    Struct_name             = strrep(Struct_name,'-','_');
    flnmrootNew{iFile,1}      = Struct_name;
    LFP_data.(Struct_name)  = [];                                          % File name as a part of a structure
    set_file                = [data_path '\' flnmroot{iFile}(1:end-4) '.set'];            % set file name 
    egf_file                = [data_path '\' flnmroot{iFile}(1:end-4) '.egf'];            % egf file name
    egf2_file               = [data_path '\' flnmroot{iFile}(1:end-4) '.egf2'];          % efg2 file name
    pos_file                = [data_path '\' flnmroot{iFile}(1:end-4) '.pos']; 
    setfile_header          = read_binary_data(set_file,'set');
    %     Reading the first eeg electrode data
    egfData.Animal_Code     = AnimalCode;
    egfData.FilePath        = flnmroot{iFile}(1:end-4);
    egfData.trial_date      = char(setfile_header.KeyValue('trial_date'));
    egfData.trial_time      = char(setfile_header.KeyValue('trial_time'));
    egfData.trial_duration  = setfile_header.KeyValue('duration','num');
    [egfData.eeg_and_time, egfData.header, egfData.downsampledFs] = read_egf_file_abdel(egf_file);                                        % Redaing the egf data without converting it into volts
    conversion_to_egf_file 
    % This script convert the raw values into voltage           
    egfData.Coversion_Factor = factor;
    egfData.eeg_and_time(:,1)= egfData.eeg_and_time(:,1) * 1000000;
     egfData.channel_num     = setfile_header.KeyValueExact(['EEG_ch_' num2str(eegNum)],'num');                         % Extracting the channel number
    egfData.channel_gain     = setfile_header.KeyValueExact(['gain_ch_', num2str(egfData.channel_num-1)],'num');
    egfData.sample_rate      = egfData.header.KeyValue('sample_rate',  'num');
    egfData.trial_start_time = egfData.eeg_and_time(1,2) - (1 / egfData.downsampledFs);                              % trial start time
    egfData.trial_end_time   = egfData.eeg_and_time(end,2);
    clear mean_egf std_egf ADC_fullscale_mv bytesBerSamp chn convFact eegChanGain factor fracAboveRange maxIntVal
    %     Reading the second EEG electrode data
    [egfData.eeg_and_time2, egfData.header2, egfData.downsampledFs] = read_egf_file_abdel(egf2_file);                                        % Redaing the egf data without converting it into volts
    conversion_to_volts_egf2                                                                                       % This script convert the raw values into voltage           
    egfData.Coversion_Factor2 = factor;
    egfData.eeg_and_time2(:,1)= egfData.eeg_and_time2(:,1) * 1000000;
     egfData.channel_num_2    = setfile_header.KeyValueExact(['EEG_ch_' num2str(eegNum_2)],'num');                                         % Extracting the channel number
    egfData.channel_gain_2    = setfile_header.KeyValueExact(['gain_ch_', num2str(egfData.channel_num_2-1)],'num');
    egfData.sample_rate2      = egfData.header2.KeyValue('sample_rate',  'num');                                          % Sampling frequency of the data                                                                % The range of EEg data
    clear mean_egf2 std_egf2 ADC_fullscale_mv bytesBerSamp chn convFact eegChanGain factor fracAboveRange maxIntVal
    LFP_data.(Struct_name)    = egfData;
    % Reading the position data 
    [posData.led_pos, posData.led_pix, posData.header] = read_pos_file(pos_file);
    [posData.xy, posData.dir, posData.speed, posData.times, posData.jumpyPercent, posData.num_LEDs, posData.head_dir] = ...
    postprocess_pos_data(posData, 1, 0.4, setfile_header);
    posData.pixels_per_metre       = posData.header.KeyValue('pixels_per_metre',  'num');
    posData.sample_rate            = posData.header.KeyValue('sample_rate',  'num');
    posData.acceleration           = makeAcceleration(posData.speed,posData.sample_rate,4);
    LFP_data.(Struct_name).PosData = posData;
end
clearvars -except LFP_data flnmroot animal AnimalCode DayStr data_path setfile_header flnmrootNew TXTfiles TXTfilesName
%% Removing local drifts of the data
for iFilt = 1 :numel(flnmroot)
   LFP_data.(flnmrootNew{iFilt,1}).eeg_filtered  = locdetrend(LFP_data.(flnmrootNew{iFilt,1}).eeg_and_time(:,1),LFP_data.(flnmrootNew{iFilt,1}).downsampledFs,[.1 .05]); 
   LFP_data.(flnmrootNew{iFilt,1}).eeg_filtered2 = locdetrend(LFP_data.(flnmrootNew{iFilt,1}).eeg_and_time2(:,1),LFP_data.(flnmrootNew{iFilt,1}).downsampledFs,[.1 .05]);
end
%% Plotting all trials
for iStruc = 1: numel(flnmrootNew)

    % Figure out how many subplots do we need
    figure(iStruc)
    
    % Plotting the data
    plot(LFP_data.(flnmrootNew{iStruc}).eeg_and_time(:,2),LFP_data.(flnmrootNew{iStruc}).eeg_filtered,'k','LineWidth',2)
    hold on
    plot(LFP_data.(flnmrootNew{iStruc}).eeg_and_time(:,2),LFP_data.(flnmrootNew{iStruc}).eeg_filtered2+500,'r','LineWidth',2)
    xlabel('Time [seconds]')
    ylabel('Amplitude [\muV]')
    lgd = legend('CA1','DG');
    set(lgd,'Box','off','Location','best','FontWeight','bold')
    suptitle(flnmroot{iStruc}(1:end-4))
    set(gca,'FontSize',20)
    set(gcf, 'Color','w','Position', [200, 200, 1049, 895])
end
clear iStruc
%% Plotting rat trajectories 
for iPos = 1: numel(flnmroot)
    Struct_name       = flnmrootNew{iPos}; 
    x                 = LFP_data.(Struct_name).PosData.xy(:,2);
    y                 = LFP_data.(Struct_name).PosData.xy(:,1);
    T                 = LFP_data.(Struct_name).PosData.times;
    figure
    plot(x,y,'.','Color','k');
    set(gcf,'Color','w')
    axis off
end
clearvars -except LFP_data flnmroot animal AnimalCode DayStr data_path setfile_header flnmrootNew TXTfiles TXTfilesName
%% Extract the epochs 
for iTxt = 1:numel(flnmrootNew)
    StrtoLook = 'trial';
    FindStr   = strfind(flnmrootNew{iTxt},StrtoLook);
    if ~isempty(FindStr)
        TrialNumI  = FindStr+5;
        TrialNum   = flnmrootNew{iTxt}(TrialNumI);
        Num        = str2num(TrialNum);
        ObjectData = load(char(TXTfilesName{Num}));
        LFP_data.(flnmrootNew{iTxt}).ObjectData = ObjectData;
    end
end