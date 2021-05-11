function [OUT] = importMAT_data(filename)

% Function that imports MAT data, in terms of first contact (FC) and last
% contact (LC) for both feet. Events are reported in samples. Furthermore,
% events are divided in LAPS.
%
% INPUT:
%   filename = .xlsx file containing MAT data
%
% OUTPUT:
%   OUT = structure with all relevant mAT data nested in.
%
% ________________________________________________________________________
%% Import data
[MAT,CLL] = xlsread(filename);

MAT_smpls = MAT(find(MAT(:,1) == 1):end,3:end); 
MAT_param = MAT(4:find(MAT(:,1) == 1)-2,3:end);

% ________________________________________________________________________
%% Define laps and right/left sides properly.
RL   = CLL(22:end,2)';
LAPS = str2double(extractBefore(RL,':'));

RL(contains(RL,'Right')) = {'Right'};
RL(contains(RL,'Left'))  = {'Left'};

% ________________________________________________________________________
%% Features strings 
% Create a string of features that I want to save as sampled (for each foot
% contact) and as final statistics (mean, std, etc). 

FEAT = CLL(6,3:end);
% first_contact         First contact [sec]
% last_contact          Last contact [sec]
% int_pressure          Integ. Pressure [p * sec]
% ftt_length            Foot length [cm]
% ftt_length_perc       Foot length [%]
% ftt_width             Foot width [cm]
% ftt_area              Foot area [cm * cm]
% ftt_angle             Foot angle [degrees]
% IO_angle              Toe in/out angle [degrees]
% toe_x                 Foot toe x location [cm]
% toe_y                 Foot toe y location [cm]
% heel_x                Foot heel x location [cm]
% heel_y                Foot heel y location [cm]
% ftt_x                 Foot center x location [cm]
% ftt_y                 Foot center y location [cm]
% SPL                   Step length [cm]
% abs_SPL               Absolute step length [cm]
% SL                    Stride length [cm]
% SW                    Stride width [cm]
% SPT                   Step time [sec]
% ST                    Stride time [sec]
% SV                    Stride velocity [cm/sec]
% DOP                   DOP [degrees]
% cycle_time            Gait cycle time [sec]
% stn_time              Stance time [sec]
% stn_perc              Stance [%]
% swn_time              Swing time [sec]
% swn_perc              Swing [%]
% SS                    Single support [sec]
% SS_perc               Single support [%]
% ini_DS                Init. D. support [sec]
% ini_DS_perc           Init. D. support [%]
% end_DS                Terminal D. support [sec]
% end_DS_perc           Terminal D. support [%]
% tot_DS                Total D. support [sec]
% tot_DS_perc           Total D. support [%]
% CISP                  CISP Time [sec]
% CISP_AP               CISP AP [%]
% CISP_ML               CISP ML [%]
% stn_COP_dist          Stance COP dist. [cm]
% SS_COP_dist           SS COP dist. [cm]
% DS_COP_dist           DP COP dist. [cm]     
% stn_COP_dist_perc     Stance COP dist. [%]
% SS_COP_dist_perc      SS COP dist. [%]
% stn_COP_pathEff       Stance COP path eff [%]
% SS_COP_pathEff        SS COP path eff [%]
% DS_COP_pathEff        DS COP path eff [%]

smpl = {'first_contact','last_contact','int_pressure','ftt_length', ...
        'ftt_length_perc','ftt_width','ftt_area','ftt_angle',...
        'IO_angle', 'toe_x','toe_y','heel_x','heel_y','ftt_x','ftt_y',...
        'SPL','abs_SPL','SL','SW','SPT','ST','SV','DOP','cycle_time', ...
        'stn_time','stn_perc','swn_time','swn_perc','SS','SS_perc',...
        'ini_DS','ini_DS_perc','end_DS','end_DS_perc','tot_DS',...
        'tot_DS_perc','CISP','CISP_AP','CISP_ML','stn_COP_dist',...
        'SS_COP_dist','DS_COP_dist','stn_COP_dist_perc',...
        'SS_COP_dist_perc','stn_COP_pathEff','SS_COP_pathEFF',...
        'DS_COP_pathEff'};
    
smpp = {'First Contact (sec.)','Last Contact (sec.)', ...
        'Integ. Pressure (p x sec.)','Foot Length (cm.)', ...
        'Foot Length %','Foot Width (cm.)','Foot Area (cm. x cm.)', ...
        'Foot Angle (degrees)','Toe In/Out Angle (degrees)', ...
        'Foot Toe X Location (cm.)','Foe Toe Y Location (cm.)', ...
        'Foot Heel X Location (cm.)','Foot Heel Y Location (cm.)', ...
        'Foot Center X Location (cm.)','Foot Center Y Location (cm.)', ...
        'Step Length (cm.)','Absolute Step Length (cm.)', ...
        'Stride Length (cm.)','Stride Width (cm.)','Step Time (sec.)', ...
        'Stride Time (sec.)','Stride Velocity (cm./sec.)', ...
        'DOP (degrees)','Gait Cycle Time (sec.)','Stance Time (sec.)',...
        'Stance %','Swing Time (sec.)','Swing %', ...
        'Single Support (sec.)','Single Support %', ...
        'Initial D. Support (sec.)','Initial D. Support %', ...
        'Terminal D. Support (sec.)','Terminal D. Support %', ...
        'Total D. Support (sec.)','Total D. Support %', ...
        'CISP Time (sec.)','CISP AP (%)','CISP ML (%)', ...
        'Stance COP Dist. (cm.)','SS COP Dist. (cm.)', ...
        'DS COP Dist. (cm.)','Stance COP Dist. %','SS COP Dist. %', ...
        'Stance COP Path Eff. %','SS COP Path Eff. %', ...
        'DS COP Path Eff. %'};

% Add to smpl some features that I know only as mean values and that I want
% to read from MAT_param.
% vel                   Velocity [cm/s]
% amb_T                 Ambulation time [sec]
% CAD                   Cadence [steps/min]
% FAP                   FAP
% mean_GVI              mean GVI
% L_GVI                 Left GVI
% R_GVI                 Right GVI
% WR                    Walking ratio [cm/(steps/min)]

matp = [smpl,'vel','amb_T','CAD','FAP','mean_GVI','L_GVI','R_GVI','WR'];

mapp = [smpp,'Velocity (cm./sec.)','Ambulation Time (sec.)', ...
       'Cadence (steps/min.)','FAP','Mean GVI','Left GVI','Right GVI', ...
       'Walk Ration (cm./(steps/min.))'];

% Final statistics

% mean                  mean value
% std                   standard deviation
% CV                    coefficient of variation

stat = {'mean','mean_L','mean_R','ratio_RL','ASI','std','std_L',...
        'std_R','CV','CV_L','CV_R'};

% ________________________________________________________________________
%% Right and Left Sides
idx_L = [find(strcmp(stat,'mean_L')),find(strcmp(stat,'std_L')),find(strcmp(stat,'CV_L'))];
idx_R = [find(strcmp(stat,'mean_R')),find(strcmp(stat,'std_R')),find(strcmp(stat,'CV_R'))];

stt = {'mean','std','CV'};
% Sampled features and statistics
for i = 1:length(smpl)
    % Right side
    mat_R = MAT_smpls(strcmp(RL,'Right'),:);
    OUT.R.(smpl{i}) = mat_R(:,strcmp(FEAT,(smpp{i})))';
    
    % Left side
    mat_L = MAT_smpls(strcmp(RL,'Left'),:);
    OUT.L.(smpl{i}) = mat_L(:,strcmp(FEAT,(smpp{i})))';
    
    % Statistics
    for j = 1:length(stt)
        OUT.R_stat.(stt{j}).(smpl{i}) = MAT_param(idx_R(j),strcmp(FEAT,(smpp{i})))';
        OUT.L_stat.(stt{j}).(smpl{i}) = MAT_param(idx_L(j),strcmp(FEAT,(smpp{i})))';
    end
end

% ________________________________________________________________________
%% Overall
stat_O = {'mean','ratio_RL','ASI','std','CV'};
idx_O  = [];
for i = 1:length(stat_O)
    idx_O = [idx_O,find(strcmp(stat,(stat_O{i})))];
end

for i = 1:length(matp)
    for j = 1:length(stat_O)
        OUT.TOT.(stat_O{j}).(matp{i}) = MAT_param(idx_O(j),strcmp(FEAT,(mapp{i})))';
    end
end

% ________________________________________________________________________
%% Laps, Time stamp and synchronization
N_LAPS = max(LAPS);
for i = 1:N_LAPS
    % first contact
    OUT.LAP.starts(i) = min(MAT_smpls(LAPS == i,1));
    
    % last contact
    OUT.LAP.ends(i)   = max(MAT_smpls(LAPS == i,2));
end

[MAT_time,CLL_time] = xlsread(filename,'Time');
CLL_time = CLL_time(6,1:end);

OUT.t    = round(MAT_time(:,strcmp(CLL_time,'Time (sec.)')),3)'; 
OUT.fc   = round(length(OUT.t)/max(OUT.t));
OUT.Sync = MAT_time(:,strcmp(CLL_time,'Sync. In'))'; 


end

