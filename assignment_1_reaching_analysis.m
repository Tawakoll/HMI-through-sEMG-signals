%% Assignment 1 - MCI reaching task: analysis
%
% The pipeline, step by step:
%   Part 1  raw EMG -> your own envelopes (conventional filtering),
%           compared with the envelopes the MCI computed live
%   Part 2  the MCI map: muscles -> cursor velocity (given in the .mat)
%   Part 3  rebuild the cursor from the envelopes and compare it with the
%           logged one
%   Part 4  performance per reach
%   BONUS   drive a simulated 2-link robot arm with your cursor
%
% Input: your group's .mat file. It contains:
%   centerout          table, one row per sample: t, raw_<muscle>,
%                      env_<muscle>, cursor_x_px, cursor_y_px, cursor_x_cm, cursor_y_cm
%   centerout_targets  table, one row per target: target_x_px, target_y_px,
%                      spawn_time, hit_time, reach_time, overshoots
%   random, random_targets   the same for the random task
%   info               sampling frequency, muscles, screen, calibration, MCI map
% (double-click a table in the Workspace to see it like a spreadsheet)
%
% Toolbox needed: Signal Processing Toolbox (butter, filtfilt, xcorr).
%
% Wherever you see "% YOUR CODE HERE", complete the code. Keep the variable
% names written in the instructions: the code given after them uses those names.

clear; close all; clc;
% To be done once, then commented
%  make_group_mat('D:\university\AY-2 2026-2027\Semester 1 fall\Biomedical Robotics\week 5\assignment1\ID_02_Session_4', 'Cometa', 'Group_1.mat');


%% 0. Settings
% -------------------------------------------------------------------------
matFile = 'Group_1.mat';   % your group's data file
task    = 'centerout';     % which task to analyze: 'centerout' or 'random'
fsDown  = 100;             % rate we use after the envelope is computed [Hz]

animateCursor     = true;  % Part 3: animate the logged and the rebuilt cursor on the screen
resetAtEachTarget = true;  % Part 3: restart the rebuilt cursor from the logged one at each new target


%% 1. Load the data
% -------------------------------------------------------------------------
% The .mat already contains only the game part of the recording, with a
% regular time vector starting at 0. (In the original txt files the time
% column is uneven, because samples reach the PC in bursts, but there are
% exactly 2000 samples per second, so the time was rebuilt from the
% sample index.)
load(matFile);                        % info, centerout, centerout_targets, random, random_targets

% Pick the task we analyze
if strcmp(task, 'centerout')
    D = centerout;          targets = centerout_targets;
elseif strcmp(task,'random')
    D = random;             targets = random_targets;

else

end

fs       = info.fs;                   % EMG sampling frequency [Hz]
muscles  = info.muscles;              % 'R biceps', 'L biceps', ...
screenPx = info.screen_px;            % screen resolution [px]
cmPerPx  = info.cm_per_px;            % [cm per px]
centerPx = info.center_px;            % center target [px]
MAX      = info.calib_max;            % value of each muscle that the MCI calls 100%

% The table has one column per muscle (e.g. raw_R_biceps). For filtering
% it is handier to have the 4 muscles side by side in one matrix, so take
% those columns out of the table by name (column order = info.muscles).
muscleCols = strrep(muscles, ' ', '_');               % 'R biceps' -> 'R_biceps'

% Extract these variables from the table D and from targets:
%   t         time [s]                                     (N x 1)
%   raw       raw EMG, one column per muscle               (N x 4)
%   envMCI   MCI envelopes, one column per muscle (0-1)  (N x 4)
%   cursorPx  cursor position [px], y grows DOWNWARD       (N x 2)
%   cursorCm  cursor position [cm]                         (N x 2)
%   targetPx  target positions [px], one row per target    (nTargets x 2)
% Hint: D{:, strcat('raw_', muscleCols)} takes the 4 raw columns at once.

% YOUR CODE HERE
%amb = strcat('raw_', muscleCols);
raw = D{:, strcat('raw_', muscleCols)};

t= D{:,'t'};
envMCI =  D{:, strcat('env_', muscleCols)};
cursorPx = D{:, {'cursor_x_px','cursor_y_px'}};
cursorCm = D{:, {'cursor_x_cm','cursor_y_cm'}}; 

% target_D =strcat(task,'_targets'); 
% not important :D
% to generalize the task once, for example if we enter task above as
% centerout or random to have this done once above :D , to improve later
% because not working
% 
targetPx = centerout_targets{:,{'target_x_px','target_y_px'}};

% Units: the raw EMG is stored in ADC counts. The sensors use a signed
% 16-bit converter (-32768 ... +32767) over an input range of +-2.5 mV, so
% one count is 2500 uV / 32768 = 0.0763 uV. We convert to microvolts so the
% plots have physical units. MAX comes from the same raw signal, so it is
% converted too (the normalized envelopes, 0-1, don't change).
uVperCount = 2500 / 32768;            % [uV per ADC count]
raw = raw * uVperCount;               % raw EMG [uV]
MAX = MAX * uVperCount;               % calibration max [uV]
N = numel(t);
fprintf('Loaded %s, %s task: %d samples = %.1f s, %d targets\n', ...
        info.system, task, N, t(end), height(targets));
fprintf('Muscles: %s\n', strjoin(muscles, ', '));


%% PART 1 - Preprocessing (conventional pipeline, as in Assignment 1)
% -------------------------------------------------------------------------
% 1) Band-pass 30-450 Hz: keeps the EMG band, removes movement artifacts
%    and offset (<30 Hz) and high-frequency noise (>450 Hz).
% 2) Rectify: EMG goes positive and negative around zero; its size (how
%    hard the muscle works) is in the absolute value.
% 3) Low-pass 5 Hz: smooths the rectified signal into an envelope that
%    follows the muscle activation, not the single spikes.
% 4) Normalize by MAX: puts every muscle on the same 0-1 scale, like the MCI.
%
% We use filtfilt (forward + backward): it has NO delay. The MCI works
% live, so it can only use past samples and its envelope has a delay.
% This is the main difference you will see in the comparison below.
%
% Steps (all 4 muscles at once: filtfilt works column by column):
%   a) design the band-pass and the low-pass filters
%   b) remove the DC offset from raw (subtract the mean of each column)
%   c) band-pass             -> emgFilt  (N x 4) [uV]
%   d) rectify               -> emgRect  (N x 4) [uV]
%   e) low-pass              -> envOurs  (N x 4) [uV]
%   f) normalize: divide each column by its MAX  -> envOurs (0-1)
%   g) filtering can dip slightly below 0: set negative values to 0

% YOUR CODE HERE
    % from matrix D  -> 1 R_bicep ,  2 L_bicep , 3 L_trap , 4 r_trap
    raw_R_bicep = raw (:,1);
    raw_L_bicep = raw (:,2);
    raw_L_trap = raw (:,3);
    raw_R_trap = raw (:,4);

%   a) design the band-pass and the low-pass filters
    [b, a] = butter(4, [30 450]/(fs/2), 'bandpass');
    [d,c] = butter(4,5/(fs/2),'low');


%   b) remove the DC offset from raw (subtract the mean of each column)
    % remove the mean for each raw muscle signal 

    dc_bias= mean(raw,1); % returns a row vector each entry is the mean of the column in the input raw
    raw_remove_DC = raw - dc_bias; 
    

%   c) band-pass             -> emgFilt  (N x 4) [uV]
    emgFilt = filtfilt(b,a,raw_remove_DC);
    filt_R_bicep = emgFilt (:,1);
    filt_L_bicep = emgFilt (:,2);
    filt_L_trap = emgFilt (:,3);
    filt_R_trap = emgFilt (:,4);

%   d) rectify               -> emgRect  (N x 4) [uV]
    emgRect = abs(emgFilt);

%   e) low-pass              -> envOurs  (N x 4) [uV]
    envOurs = filtfilt(d,c,emgRect);

%   g) filtering can dip slightly below 0: set negative values to 0
    % for any matrix A , doing A(A = 0) for example indexes all elements of matrix A and checks if = 0 
    % then you can choose to replace this value as done in the line bellow. For negative values.
    envOurs(envOurs<0)= 0;
    
%   f) normalize: divide each column by its MAX  -> envOurs (0-1)
    envOurs_norm= envOurs./max(envOurs); % takes column1/its max and so on..
    


% --- Figure 1: the preprocessing steps for each muscle, on a 3 s window
% with some activity (easier to see). For each muscle plot:
%   - raw and filtered EMG [uV]
%   - rectified EMG with the envelope on top [uV]
% Label the axes with units.

% YOUR CODE HERE
figure('Name','Pre-processing steps: raw , filtered, rectified, Envelop')  

subplot(4,4,1)
plot(t,raw_R_bicep,'y');
ylabel('right bicep')
title('Raw EMG')
xlim([0 3]);

subplot(4,4,5)
plot(t,raw_L_bicep,'g');
ylabel('left bicep');
xlim([0 3]);

subplot(4,4,9)
plot(t,raw_L_trap,'c');
ylabel('left trapezius');
xlim([0 3]);

subplot(4,4,13)
plot(t,raw_R_trap,'color', [1 0.7 0]);
ylabel('right trapezius');
xlabel('Time');
xlim([0 3]);

%%%%%%%%%%%% after band pass filtered
subplot(4,4,2)
plot(t,filt_R_bicep,'y');
ylabel('right bicep')
title('Filtered EMG 30 to 450 Hz')
xlim([0 3]);

subplot(4,4,6)
plot(t,filt_L_bicep,'g');
ylabel('left bicep');
xlim([0 3]);

subplot(4,4,10)
plot(t,filt_L_trap,'c');
ylabel('left trapezius');
xlim([0 3]);

subplot(4,4,14)
plot(t,filt_R_trap,'color', [1 0.7 0]);
ylabel('right trapezius');
xlabel('Time');
xlim([0 3]);

%%%%%%%%%%%% after Full wave rectifier, aka absolute value :D
subplot(4,4,3)
plot(t,emgRect(:,1),'y');
ylabel('right bicep')
title('Rectified')
xlim([0 3]);

subplot(4,4,7)
plot(t,emgRect(:,2),'g');
ylabel('left bicep');
xlim([0 3]);

subplot(4,4,11)
plot(t,emgRect(:,3),'c');
ylabel('left trapezius');
xlim([0 3]);

subplot(4,4,15)
plot(t,emgRect(:,4),'color', [1 0.7 0]);
ylabel('right trapezius');
xlabel('Time');
xlim([0 3]);

%%%%%%%%%%%% after linear envelop , by 5 Hz low pass filter
subplot(4,4,4)
plot(t,envOurs(:,1),'y');
ylabel('right bicep')
title('Linear Envelop')
xlim([0 3]);

subplot(4,4,8)
plot(t,envOurs(:,2),'g');
ylabel('left bicep');
xlim([0 3]);

subplot(4,4,12)
plot(t,envOurs(:,3),'c');
ylabel('left trapezius');
xlim([0 3]);

subplot(4,4,16)
plot(t,envOurs(:,4),'color', [1 0.7 0]);
ylabel('right trapezius');
xlabel('Time');
xlim([0 3]);

% --- Compare your envelopes with the MCI ones (envMCI)
% For each muscle:
%   - correlation between the two envelopes (corrcoef): same shape?
%   - delay between them: the lag of the peak of xcorr (convert it to ms)
% Then plot the two envelopes on top of each other (a 30 s window is enough).

% YOUR CODE HERE
a1=corrcoef(envMCI(:,1),envOurs_norm(:,1));
alpha1=corrcoef(envMCI(:,1),envOurs_norm(:,1));


[amb1,alags1] =xcorr(envMCI(:,1),envOurs_norm(:,1));
[maxValue1, maxIndex1] = max(amb1);

delay1 = (alags1(maxIndex1)/fs)*1000;
% to get seconds divide by fs, delay1 is in samples so samples/ (samples/ second) we get seconds. 
% Then multiply by 1000 to get milliseconds


[amb2,alags2] =xcorr(envMCI(:,2),envOurs_norm(:,2));
[maxValue2, maxIndex2] = max(amb2);

delay2 = (alags2(maxIndex2)/fs)*1000; 

[amb3,alags3] =xcorr(envMCI(:,3),envOurs_norm(:,3));
[maxValue3, maxIndex3] = max(amb3);

delay3 = (alags3(maxIndex3)/fs)*1000; 

[amb4,alags4] =xcorr(envMCI(:,4),envOurs_norm(:,4));
[maxValue4, maxIndex4] = max(amb4);

delay4 = (alags4(maxIndex4)/fs)*1000; 


figure('Name','Recorded MCI Envelopes VS our Envelopes (Normalized)'); 
subplot(4,1,1);
plot(t,envMCI(:,1),t,envOurs_norm(:,1));
legend('MCI Env', 'Our Env Norm');
xlabel('Time');
ylabel('Right bicep');
xlim([0 30]);

subplot(4,1,2);
plot(t,envMCI(:,2),t,envOurs_norm(:,2));
legend('MCI Env', 'Our Env Norm');
xlabel('Time');
ylabel('left bicep');
xlim([0 30]);

subplot(4,1,3);
plot(t,envMCI(:,3),t,envOurs_norm(:,3));
legend('MCI Env', 'Our Env Norm');
xlabel('Time');
ylabel('left trap');
xlim([0 30]);

subplot(4,1,4);
plot(t,envMCI(:,4),t,envOurs_norm(:,4));
legend('MCI Env', 'Our Env Norm');
xlabel('Time');
ylabel('right trap');
xlim([0 30]);

%%% comments after comparing both envelopes. 
% First It can be seen clearly from the figures
% that our signal where we reached the envelope by a lowpass filter of
% 5Hz was noticeably smoother than the signal recorded ( envelope at 9 Hz) which makes
% sense. As higher frequency envelopes tends to be more jittery "Jumpy"
% while lower frequency LP filter to get the envelope is more smooth. 
% About the lag / leading of one signal w.r.to the other , there is a lead
% between our envelope and the recorded MCI as expected and to be confirmed
% by exact number from the correlation analysis.
% From the xcorr value of delay1. Right bicep delay was 36.5 ms. Similarly
% the other muscle signals follow.

% --- Downsample to 100 Hz
% WHY only now: the envelope is a slow signal (below 5 Hz), so 100  Hz is
% plenty. Doing it BEFORE filtering would alias the EMG (which goes up to
% 450 Hz) and destroy it. The cursor also only updates ~60 times per second.
step = fs / fsDown;                        % keep 1 sample every 20 , 2000/100
dt   = 1 / fsDown;                         % time between downsampled samples [s]
% Keep one sample every 'step' and create:
%   tD, envOursD, envMCID, cursorCmD
% YOUR CODE HERE
tD = 0 : dt:259;



% %% PART 2 - The MCI map (given)
% % -------------------------------------------------------------------------
% % The MCI is in velocity control: the muscles set the cursor SPEED.
% % Each axis is driven by a pair of muscles, and the speed is proportional
% % to their DIFFERENCE, so contracting both together cancels out:
% %
% %     x speed = gain * (muscle A - muscle B)      (biceps pair)
% %     y speed = gain * (muscle C - muscle D)      (trapezius pair)
% %
% % gain = Max Cursor Speed (50 cm/s): a muscle at 100% with the other one
% % relaxed moves the cursor at 50 cm/s. The .mat stores the map (info.W) as a
% % matrix W, one row per axis and one column per muscle, so that for all
% % samples at once:
% %
% %     velocity (N x 2) = gain * envelopes (N x 4) * W'
% gain = info.gain_cm_s;                     % [cm/s]
% W    = info.W;                             % 2 x 4
% fprintf('\nPart 2 - MCI map\n');
% fprintf('  x speed        = %g * (%s)\n', gain, info.map_x);
% fprintf('  y speed (down) = %g * (%s)\n', gain, info.map_y);
% 
% % Dead zone: the MCI ignores any muscle below 5% of its max, so that
% % resting noise does not move the cursor. Without this step the rebuilt
% % cursor creeps away from the real one even when the subject is relaxed.
% deadZone = 0.05;
% 
% % 1) Apply the dead zone to envOursD and to envMCID (values below
% %    deadZone become 0)
% % 2) Compute the cursor velocity [cm/s] with the formula above:
% %      velOurs  from your envelopes      (N x 2)
% %      velMCI  from the MCI envelopes  (N x 2)
% 
% % YOUR CODE HERE
% 
% 
% 
% %% PART 3 - Rebuild the cursor from the envelopes
% % -------------------------------------------------------------------------
% % Velocity from the map (Part 2), then integrate (position = previous position +
% % velocity * dt). Like the real cursor, ours can't leave the screen, so we
% % clip it at the edges at every step.
% %
% % Integrating adds up every small error, so the rebuilt cursor slowly
% % drifts away from the real one (try resetAtEachTarget = false to see it).
% % With the reset, we restart from the logged position at each new target,
% % so the error only builds up within one reach.
% screenCm = info.screen_cm;                       % screen size [cm]
% spawnT   = targets.spawn_time;                   % when each target appears [s]
% 
% % Rebuild two cursors, posOurs and posMCI [cm] (nD x 2, nD = numel(tD)),
% % from velOurs and velMCI:
% %   1) start both from the logged position: cursorCmD(1,:)
% %   2) for each sample k = 2..nD: previous position + velocity * dt
% %   3) clip to the screen: between [0 0] and screenCm
% %   4) if resetAtEachTarget is true: at the first sample after each
% %      spawnT, set both positions back to the logged one
% % Then compute, for both, the distance from the logged cursor at each
% % sample (errOurs, errMCI [cm]) and print its mean.
% 
% % YOUR CODE HERE
% 
% 
% % Targets, converted to cm -> targCm (nTargets x 2)
% 
% % YOUR CODE HERE
% 
% 
% % --- Animation: the cursor moving on the screen, target by target (given)
% % Black = the real (logged) cursor, red = the one we rebuilt from the EMG.
% % The current target is the green circle. It plays in real time.
% if animateCursor
%     nD = numel(tD);
%     % which target is on screen at each sample = the last one spawned
%     curTarget = zeros(nD, 1);
%     for s = 1:numel(spawnT)
%         curTarget(tD >= spawnT(s)) = s;
%     end
% 
%     figure('Name', 'Part 3 - cursor animation');
%     rectangle('Position', [0 0 screenCm]); hold on;            % the screen
%     hTarg = plot(nan, nan, 'go', 'MarkerSize', 25, 'LineWidth', 2);
%     hLogTrail = plot(nan, nan, 'Color', [0.6 0.6 0.6]);
%     hLog  = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
%     hOurs = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%     set(gca, 'YDir', 'reverse');                               % screen y grows downward
%     axis equal; xlim([0 screenCm(1)]); ylim([0 screenCm(2)]);
%     xlabel('x [cm]'); ylabel('y [cm]');
%     legend([hLog hOurs hTarg], {'logged cursor', 'rebuilt cursor', 'target'}, 'Location', 'eastoutside');
%     tic; k = 1;
%     while k <= nD
%         if curTarget(k) > 0
%             set(hTarg, 'XData', targCm(curTarget(k),1), 'YData', targCm(curTarget(k),2));
%         end
%         k0 = max(1, k - 200);              % trail = last 2 s of the logged cursor
%         set(hLogTrail, 'XData', cursorCmD(k0:k,1), 'YData', cursorCmD(k0:k,2));
%         set(hLog,  'XData', cursorCmD(k,1), 'YData', cursorCmD(k,2));
%         set(hOurs, 'XData', posOurs(k,1),   'YData', posOurs(k,2));
%         title(sprintf('t = %.1f s   target %d', tD(k), curTarget(k)));
%         drawnow;
%         k = max(k + 1, round(toc * fsDown) + 1);   % jump to the sample of "now"
%         pause(max(0, (k - 1) / fsDown - toc));     % ahead of the clock? wait
%     end
% end
% 
% 
% 
% %% PART 4 - Performance per reach
% % -------------------------------------------------------------------------
% % Each trial goes from the moment the target appears (spawn_time) to the
% % moment it is hit (hit_time). In center-out, every other target is the
% % center (the way back): we look at the outward reaches only. In the
% % random task all targets count.
% %   1) select the reaches to use: target hit (hit_time is not NaN) and
% %      target not at the center (targetPx differs from centerPx)
% %   2) for each of them, plot in order (a learning curve):
% %      - the reach time [s]            (targets.reach_time)
% %      - the overshoots [count]        (targets.overshoots)
% 
% % YOUR CODE HERE
