function [EVT_OUT] = gait_evt_VR(IN,side)

% Function that computes gait heel strikes (HS) and toe offs (TO) using the
% vertical trajectory
%
% INPUT:
%   IN = input VR data structure
%
% OUTPUT: 
%   side = "R" for right side and "L" for left side.
%
% ________________________________________________________________________
%% Import data
% Data are selected depending on the side, then filtered and detrended
% before further analysis. Only the vertical trajectory is taken into
% account.
if (strcmp(side,'R'))
    pos = IN.TR1.p(3,:);
else
    pos = IN.TR2.p(3,:);
end

fc  = 90;
[f1,f2] = butter(3,12/(fc/2),'low');
pos = filtfilt(f1,f2,pos')';

pos = detrend(pos,'linear');

% ________________________________________________________________________
%% Normalize signal
pos = (pos - min(pos))./(max(pos) - min(pos));

% Samples that are lower then a certain threshold are levered
thL = nanmean(pos(1:100)) - 0.75*nanstd(pos);
pos(pos <= thL) = thL;

% ________________________________________________________________________
%% Event detection
% Output vectors
EVT.HS = [];   
EVT.TO = [];

% TO : TOs are the first peaks for each double-peak wave. TOs that are
% closer than 45 samples are discarded. First element is always taken, so 
% its dTO value is set higher than 45.
% First peak 
th1 = nanmean(pos) + 0.2*nanstd(pos);
[~,TO] = findpeaks(pos,'MinPeakHeight',th1);
dTO = [100 diff(TO)];

L = length(dTO);
i = 1;
while (i <= L)
    if (dTO(i) <= 45)
        TO(i) = [];
        dTO = [100 diff(TO)];
        L = length(dTO);
        i = 1;
    else
        i = i + 1;
    end
end

% HS : HSs are the second peaks for each double-peak wave, so they are the
% even detected peaks. Look for a HS event between two TOs.
for i = 2:length(TO)
    M = TO(i) - TO(i-1);
    tmp = pos(:,TO(i-1) + 10:TO(i-1) + round(M/2) - 1);
    [p_tmp,HS_tmp] = findpeaks(tmp,'MinPeakProminence',0.02);
    HS_tmp = HS_tmp(p_tmp == max(p_tmp));
    % In case no peak is found between two TOs, then a peak is assigned
    % with another method: baseline is found and the HS is halfway between
    % TO(i-1) and the next baseline.
    if (isempty(HS_tmp))
        HS_tmp = round(M/4);
    end
    %
    HS(i-1) = HS_tmp + TO(i-1) + 10 - 1;
end

% Events should always start with TO events and end with HS events. So TO
% events after last HS are discarded.
TO(TO >= HS(end)) = [];

% Output
EVT_OUT.HS = TO;
EVT_OUT.TO = HS;


end

