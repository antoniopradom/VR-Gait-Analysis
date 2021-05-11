function [FEAT] = gait_feat_VR(data_VR)

% Function to compute gait spatiotemporal parameters using VR data from
% feet trakers. Works on LAPS
%
% INPUT: 
%   data_VR = VR data structure
%
% OUTPUT:
%   FEAT: structure containing all the evaluated spatiotemporal features 
%
% ________________________________________________________________________
%% Import data
PEL   = data_VR.TR3.p;
RFT   = data_VR.TR1.p;
LFT   = data_VR.TR2.p;
HEAD  = data_VR.CAM;

PEL_q = data_VR.TR3.q;
TRK_q = data_VR.TR4.q;

EVT_R = data_VR.EVT_R;
EVT_L = data_VR.EVT_L;

% Sampling frequency
fc = 90;

% Initialize output:
feat = {'ST','SL','SW','GS','TH','MED','PWR','STC','SWG','SI', ...
        'TRKx','TRKy','TRKz','tSP','tS','tRS'};
for i = 1:length(feat)
    FEAT.(feat{i}) = [];
end

% _________________________________________________________________________
%% Stride Time/Length/Width/
% Stride Time
% Time in seconds elapsed between consecutive heel strikes.
RST = [];  LST = [];
RST = [RST diff(EVT_R.HS)/fc];
LST = [LST diff(EVT_L.HS)/fc];
FEAT.ST  = [RST LST];

% Timestamps
FEAT.tS  = data_VR.t([EVT_R.HS(2:end) EVT_L.HS(2:end)]);
FEAT.tRS = data_VR.t(EVT_R.HS(2:end));


% Stride Length
% distance between corresponding successive heel points of opposite feet, 
% measured parallel to the direction of walking.
% Table is specified on the paper for sensor to heel transformation 
RSL = [];
for i = 2:length(EVT_R.HS)
    tmp = [RFT(1:2,EVT_R.HS(i))' ; RFT(1:2,EVT_R.HS(i-1))'];
    RSL = [RSL pdist(tmp)];
end
LSL = [];
for i = 2:length(EVT_L.HS)
    tmp = [LFT(1:2,EVT_L.HS(i))' ; LFT(1:2,EVT_L.HS(i-1))'];
    LSL = [LSL pdist(tmp)];
end
FEAT.SL = [RSL LSL];

% Stride Width
% Median distance between left and right foot during each stride. Distance
% is computed as euclidean distance between feet trajectories.
% Table is specified on the paper for sensor to heel transformation 
for i = 1:length(RFT)
    d_tmp(i) = pdist([RFT(1:2,i)'; LFT(1:2,i)']);
end
for i = 2:length(EVT_R.HS)
    % For each interval between two consecutive HS take d_tmp and compute
    % the median value
    SW_tmp = d_tmp(:,EVT_R.HS(i-1):EVT_R.HS(i));
    FEAT.SW  = [FEAT.SW nanmedian(SW_tmp)];
end

% ________________________________________________________________________
%% Gait Velocity
FEAT.GS = [RSL LSL]./[RST LST];

% ________________________________________________________________________
%% FREQUENCY FEATURES
% Median frequency and power for each stride window, both for right and
% left sides
for i = 2:length(EVT_R.HS)
    XR = RFT(3,EVT_R.HS(i-1):EVT_R.HS(i));
    % 
    L = length(XR);
    F_XR2 = abs(fft(XR)/L);
    F_XR1 = F_XR2(1:floor(L/2)+1);
    F_XR1(2:end-1) = 2*F_XR1(2:end-1);
    f = fc*(0:(L/2))/L;
    %
    [MED_R(:,i-1),PWR_R(:,i-1)] = medfreq(F_XR1,f);
end
for i = 2:length(EVT_L.HS)
    XL = LFT(3,EVT_L.HS(i-1):EVT_L.HS(i));
    % 
    L = length(XL);
    F_XL = fft(XL);
    F_XL2 = abs(F_XL/L);
    F_XL1 = F_XL2(1:floor(L/2)+1);
    F_XL1(2:end-1) = 2*F_XL1(2:end-1);
    f = fc*(0:(L/2))/L;
    %
    [MED_L(:,i-1),PWR_L(:,i-1)] = medfreq(F_XL1,f);
end
FEAT.MED = [MED_R MED_L];
FEAT.PWR = [PWR_R PWR_L];

% ________________________________________________________________________
%% Gait Symmetry
% Computed on SL and ST for now. Start from the first event, right or left.
% Then, search for the nearest nex element on the other side.
R_mean = mean(RSL);
L_mean = mean(LSL);
FEAT.SI = 100*abs(R_mean - L_mean)/(0.5*(R_mean + L_mean));
if (isempty(FEAT.SI))
    FEAT.SI = NaN;
end

% ________________________________________________________________________
%% Heading
for i = 2:length(EVT_R.HS)
    tmp = PEL(:,EVT_R.HS(i-1):EVT_R.HS(i));
    A = pdist([tmp(1,1) tmp(2,1); tmp(1,end) tmp(2,end)]);
    B = pdist([tmp(1,1) tmp(2,1); tmp(1,end) tmp(2,1)]);
    %
    TH_R(i-1) = acosd(B/A);
end
for i = 2:length(EVT_L.HS)
    tmp = PEL(:,EVT_L.HS(i-1):EVT_L.HS(i));
    A = pdist([tmp(1,1) tmp(2,1); tmp(1,end) tmp(2,end)]);
    B = pdist([tmp(1,1) tmp(2,1); tmp(1,end) tmp(2,1)]);
    %
    TH_L(i-1) = acosd(B/A);
end
FEAT.TH = [TH_R TH_L];

% ________________________________________________________________________
%% Stance and Swing Percentages
% Right Cycle
RSTC = 100*(EVT_R.TO(2:end) - EVT_R.HS(1:end-1))./diff(EVT_R.TO);
RSWG = 100*(EVT_R.HS(1:end-1) - EVT_R.TO(1:end-1))./diff(EVT_R.TO);

% Left Cycle
LSTC = 100*(EVT_L.TO(2:end) - EVT_L.HS(1:end-1))./diff(EVT_L.TO);
LSWG = 100*(EVT_L.HS(1:end-1) - EVT_L.TO(1:end-1))./diff(EVT_L.TO);

FEAT.STC = [RSTC LSTC];   
FEAT.SWG = [RSWG LSWG];

% ________________________________________________________________________
%% Trunk Sway
% Difference, within each stride event between TRK inclination and PEL
% inclination. Delta between minimum and maximum sway is taken as feature.
R_TRKx = [];    R_TRKy = [];    R_TRKz = [];
for i = 2:length(EVT_R.HS)
    tmp = PEL_q.select(EVT_R.HS(i-1):EVT_R.HS(i)) - ...
          TRK_q.select(EVT_R.HS(i-1):EVT_R.HS(i));
    tmp = 180/pi*EA(quaternion(tmp),'YXZ').getValue;
    R_TRKx = [R_TRKx abs(max(tmp(:,2)))];
    R_TRKy = [R_TRKy abs(max(tmp(:,1)))];
    R_TRKz = [R_TRKz abs(max(tmp(:,3)))];
end
L_TRKx = [];    L_TRKy = [];    L_TRKz = [];
for i = 2:length(EVT_L.HS)
    TRK_tmp = TRK_q.select(EVT_L.HS(i-1):EVT_L.HS(i));
    PEL_tmp = PEL_q.select(EVT_L.HS(i-1):EVT_L.HS(i));
    SWAY = PEL_tmp*TRK_tmp.inv;
    
    tmp = 180/pi*EA(quaternion(SWAY),'YXZ').getValue;
    L_TRKx = [L_TRKx abs(max(tmp(:,2)))];
    L_TRKy = [L_TRKy abs(max(tmp(:,1)))];
    L_TRKz = [L_TRKz abs(max(tmp(:,3)))];
end
FEAT.TRKx = [R_TRKx L_TRKx];
FEAT.TRKy = [R_TRKy L_TRKy];
FEAT.TRKz = [R_TRKz L_TRKz];

end
