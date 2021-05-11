function [OUT] = importVR_data(filename)

% Import data from .csv files created during VR tests in Unity environment.
% The output will be a structure.
% 
% INPUT:
%   filename = .csv file containing VR data
%
% OUTPUT: 
%   OUT = output structure, with all VR data nested inside

% ________________________________________________________________________
%% Import data
Mat = readmatrix(filename);
Opt = detectImportOptions(filename);
Opt = Opt.VariableNames;

% ________________________________________________________________________
%% Create output structure
STR = {'CAM','TR1','TR2','TR3','TR4'};

% If csv file contains "Tracker_0X" with X going from 1 to 4, instead of
% "Controller_right" and "Controller_left" labels, we need to change that.
if ~isempty(find(strcmp(Opt,'Tracker_03_x_')))
    TRK = {'Tracker_01','Tracker_02','Tracker_03','Tracker_04'};
    CNT = {'Controller_right_','Controller_left_','Tracker_01','Tracker_02'};
    for i = 1:length(TRK)
        Opt(strcmp(Opt,[TRK{i} '_x_']))  = {[CNT{i} '_x_']};
        Opt(strcmp(Opt,[TRK{i} '_y_']))  = {[CNT{i} '_y_']};
        Opt(strcmp(Opt,[TRK{i} '_z_']))  = {[CNT{i} '_z_']};
        Opt(strcmp(Opt,[TRK{i} '_ow_'])) = {[CNT{i} '_ow_']};
        Opt(strcmp(Opt,[TRK{i} '_ox_'])) = {[CNT{i} '_ox_']};
        Opt(strcmp(Opt,[TRK{i} '_oy_'])) = {[CNT{i} '_oy_']};
        Opt(strcmp(Opt,[TRK{i} '_oz_'])) = {[CNT{i} '_oz_']};
    end
end

CVS = {'Camera_eye_','Controller_right_','Controller_left_', ...
       'Tracker_01','Tracker_02'};

for i = 1:length(STR)
    % Position
    OUT.(STR{i}).p(:,1) = Mat(:,strcmp(Opt,[CVS{i} '_x_']));
    OUT.(STR{i}).p(:,2) = Mat(:,strcmp(Opt,[CVS{i} '_y_']));
    OUT.(STR{i}).p(:,3) = Mat(:,strcmp(Opt,[CVS{i} '_z_']));
    % Rotation (scalar as 4th term)
    OUT.(STR{i}).q(:,1) = Mat(:,strcmp(Opt,[CVS{i} '_ox_']));
    OUT.(STR{i}).q(:,2) = Mat(:,strcmp(Opt,[CVS{i} '_oy_']));
    OUT.(STR{i}).q(:,3) = Mat(:,strcmp(Opt,[CVS{i} '_oz_']));
    OUT.(STR{i}).q(:,4) = Mat(:,strcmp(Opt,[CVS{i} '_ow_']));
    OUT.(STR{i}).q = OUT.(STR{i}).q';
end

% ________________________________________________________________________
%% Synchronization and timestamp
OUT.Sync = Mat(:,strcmp(Opt,'Sync'));
OUT.t    = Mat(:,strcmp(Opt,'T'));

% ________________________________________________________________________
%%  Constant sampling frequency
% Sampling frequency is not constant at now: having it constant would be
% better for further analysis. 
% Sampling frequency will be imposed at 90Hz and data will be resampled. 
% Moreover, data are cut when Sync is on.
s = {'CAM','TR1','TR2','TR3','TR4'};

[OUT.t, index] = unique(OUT.t);
t_exp = (min(OUT.t):(1/90):max(OUT.t))';

OUT.Sync = interp1(OUT.t,OUT.Sync(index),t_exp)';
for i = 1:length(s)
    % p 
    OUT.(s{i}).p = interp1(OUT.t,OUT.(s{i}).p(index,:),t_exp)';
    OUT.(s{i}).p = OUT.(s{i}).p(:,OUT.Sync == 1);
    
    % q
    OUT.(s{i}).q = interp1(OUT.t,OUT.(s{i}).q(:,index)',t_exp)';
    OUT.(s{i}).q = OUT.(s{i}).q(:,OUT.Sync == 1);
end

OUT.t = t_exp';
OUT.t = OUT.t(:,OUT.Sync == 1);

end
