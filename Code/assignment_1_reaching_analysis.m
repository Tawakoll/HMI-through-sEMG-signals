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
matFile = 'sEMG live recording data/Group_1.mat';   % your group's data file
task    = 'centerout';     % which task to analyze: 'centerout' or 'random'
fsDown  = 100;             % rate we use after the envelope is computed [Hz]

animateCursor     = true;  % Part 3: animate the logged and the rebuilt cursor on the screen
resetAtEachTarget = true;  % Part 3: restart the rebuilt cursor from the logged one at each new target
animateArm        = false;  % Part 5: animate the robotic arm on the screen

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
targetPx = targets{:,{'target_x_px','target_y_px'}};

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
    [b, a] = butter(4, [10 450]/(fs/2), 'bandpass');
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
    % for any matrix A , doing A(A = 0)  indexes all elements of matrix A and checks if aij = 0 
    % then you can choose to replace this value as done in the line bellow. For negative values.
    envOurs(envOurs<0)= 0;
%   f) normalize: divide each column by its MAX  -> envOurs (0-1)
    envOurs_norm= (envOurs./MAX)*1.1 ; % takes column1/its max and so on..
     
% --- Figure 1: the preprocessing steps for each muscle, on a 3 s window
% with some activity (easier to see). For each muscle plot:
%   - raw and filtered EMG [uV]
%   - rectified EMG with the envelope on top [uV]
% Label the axes with units.

% YOUR CODE HERE
figure('Name','Pre-processing steps: raw , filtered, rectified, Envelop')  

subplot(4,3,1)
plot(t,raw_R_bicep,'y');
ylabel('right bicep [uV]')
title('Raw EMG')
xlim([1.5 4.5]);

subplot(4,3,4)
plot(t,raw_L_bicep,'g');
ylabel('left bicep [uV]');
xlim([1.5 4.5]);

subplot(4,3,7)
plot(t,raw_L_trap,'c');
ylabel('left trapezius [uV]');
xlim([1.5 4.5]);

subplot(4,3,10)
plot(t,raw_R_trap,'color', [1 0.7 0]);
ylabel('right trapezius [uV]');
xlabel('Time [s]');
xlim([1.5 4.5]);

%%%%%%%%%%%% after band pass filtered
subplot(4,3,2)
plot(t,filt_R_bicep,'y');
title('Filtered EMG 10 to 450 Hz')
xlim([1.5 4.5]);

subplot(4,3,5)
plot(t,filt_L_bicep,'g');
xlim([1.5 4.5]);

subplot(4,3,8)
plot(t,filt_L_trap,'c');
ylabel('left trapezius');
xlim([1.5 4.5]);

subplot(4,3,11)
plot(t,filt_R_trap,'color', [1 0.7 0]);
xlabel('Time [s]');
xlim([1.5 4.5]);

%%%%%%%%%%%% after Full wave rectifier, aka absolute value :D
subplot(4,3,3)
plot(t, emgRect(:,1), 'Color', [0.8 0.8 0.8]); 
hold on;
plot(t, envOurs(:,1), 'Color', 'r');
hold off;
legend({'Rectified EMG', 'Envelope'});

title('Rectified EMG , 5Hz Envelope')
xlim([1.5 4.5]);

subplot(4,3,6)
plot(t,emgRect(:,2),'Color',[0.8 0.8 0.8]);
hold on;
plot(t, envOurs(:,2), 'Color', 'r');
hold off;
xlim([1.5 4.5]);

subplot(4,3,9)
plot(t,emgRect(:,3),'Color',[0.8 0.8 0.8]);
hold on;
plot(t, envOurs(:,3), 'Color', 'r');
hold off;
xlim([1.5 4.5]);

subplot(4,3,12)
plot(t,emgRect(:,4),'Color',[0.8 0.8 0.8]);
hold on;
plot(t, envOurs(:,4), 'Color', 'r');
hold off;
xlabel('Time [s]');
xlim([1.5 4.5]);

%%%%%%%%%%%% after linear envelop , by 5 Hz low pass filter
% subplot(4,4,4)
% plot(t,envOurs(:,1),'r');
% ylabel('right bicep [uV]')
% title('Linear Envelop')
% xlim([0 3]);
% 
% subplot(4,4,8)
% plot(t,envOurs(:,2),'r');
% ylabel('left bicep [uV]');
% xlim([0 3]);
% 
% subplot(4,4,12)
% plot(t,envOurs(:,3),'r');
% ylabel('left trapezius [uV]');
% xlim([0 3]);
% 
% subplot(4,4,16)
% plot(t,envOurs(:,4),'r');
% ylabel('right trapezius [uV]');
% xlabel('Time [s]');
% xlim([0 3]);

% --- Compare your envelopes with the MCI ones (envMCI)
% For each muscle:
%   - correlation between the two envelopes (corrcoef): same shape?
%   - delay between them: the lag of the peak of xcorr (convert it to ms)
% Then plot the two envelopes on top of each other (a 30 s window is enough).

% YOUR CODE HERE
% alpha1=corrcoef(envMCI(:,1),envOurs_norm(:,1));
% 
% 


% to get seconds divide by fs, delay1 is in samples so samples/ (samples/ second) we get seconds. 
% Then multiply by 1000 to get milliseconds
fprintf('\nPart 1 -  Correlation and lag values\n');
 delay = zeros(size(muscles,2),1);
for m = 1 : size(muscles,2) % 4 = columns for each muscle
    alpha=corrcoef(envMCI(:,m),envOurs_norm(:,m));


    [R,alags] =xcorr(envMCI(:,m),envOurs_norm(:,m));

    %R is basically calculated by multiplying both signals, when they are
    %most similar the value is larger and vise versa.
    [maxValue1, maxIndex1] = max(R); % find the maximum value , correlation value between the two muscles, return this value with its index. We then use the max correlation value where the signals are most similar to calculate the lag.
  
    delay(m) = (alags(maxIndex1)/fs)*1000;
    fprintf('Correlation coefficient for %s [0 - 1] = %.2f \n',muscles{m}, alpha(2,1));
    
    
    

    if delay(m) < 0 
        fprintf('MCI leads our %s signal by  %.2f milliseconds\n',muscles{m}, abs(delay(m)));
    else
        fprintf('MCI lags our %s signal by  %.2f milliseconds\n',muscles{m},delay(m));
    
    end
 fprintf('\n');
%fprintf('  MCI error mean = %.2f cm\n', mean_errMCI);
end



% [amb2,alags2] =xcorr(envMCI(:,2),envOurs_norm(:,2));
% [maxValue2, maxIndex2] = max(amb2);
% 
% delay2 = (alags2(maxIndex2)/fs)*1000; 
% 
% [amb3,alags3] =xcorr(envMCI(:,3),envOurs_norm(:,3));
% [maxValue3, maxIndex3] = max(amb3);
% 
% delay3 = (alags3(maxIndex3)/fs)*1000; 
% 
% [amb4,alags4] =xcorr(envMCI(:,4),envOurs_norm(:,4));
% [maxValue4, maxIndex4] = max(amb4);
% 
% delay4 = (alags4(maxIndex4)/fs)*1000; 


figure('Name','Live BoMI Envelopes VS our Envelopes (Normalized)'); 

subplot(4,1,1);
plot(t,envMCI(:,1),t,envOurs_norm(:,1));
subtitle('Right biceps');
legend('MCI Env', 'Our Env Norm');
ylabel('Activation[0-1]');
xlim([15 45]);

subplot(4,1,2);
plot(t,envMCI(:,2),t,envOurs_norm(:,2));
subtitle('Left biceps');
ylabel('Activation[0-1]');
xlim([15 45]);

subplot(4,1,3);
plot(t,envMCI(:,3),t,envOurs_norm(:,3));
subtitle('Left trapezius');
ylabel('Activation[0-1]');
xlim([15 45]);

subplot(4,1,4);

plot(t,envMCI(:,4),t,envOurs_norm(:,4));
subtitle('Right trapezius');

ylabel('Activation[0-1]');
xlim([15 45]);

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



%% vectorization implementation, % Keep one sample every 'step' and create:  tD, envOursD, envMCID, cursorCmD
%  basically the idea is to create vector index, from 1 to
% reach all number of samples ND. And  for all columns ':' index all the
% rows from 1 to N take these values into our matrices.
index = 1 : step :N;
envMCID   = envMCI(index, :) ;
envOursD  = envOurs_norm (index, :) ;
cursorCmD = cursorCm (index, :) ;

nD = size(envOursD,1); % new number of samples
tD = (0 : dt:(nD-1)*dt)';

%% Archived the C/C++ :D 'implementation' to Keep one sample every 'step' and create: tD, envOursD, envMCID, cursorCmD
% envMCID(ND,4)= 0 ;
% envOursD(ND,4)= 0 ;
% cursorCmD(ND,4)= 0 ;
%
% for j = 1 : 1: 4 % to loop the column
%     for i = 1 : step : ND  % to loop the rows
% 
%         envOursD(i,j) = envOurs_norm(i,j);
%         envMCID(i,j)  = envMCI(i,j);
%         if (i <= 2) ; cursorCmD(i,j) = cursorCm(i,j); end
% 
%     end
% end

%% PART 2 - The MCI map (given)
% -------------------------------------------------------------------------
% The MCI is in velocity control: the muscles set the cursor SPEED.
% Each axis is driven by a pair of muscles, and the speed is proportional
% to their DIFFERENCE, so contracting both together cancels out:
%
%     x speed = gain * (muscle A - muscle B)      (biceps pair)
%     y speed = gain * (muscle C - muscle D)      (trapezius pair)
%
% gain = Max Cursor Speed (50 cm/s): a muscle at 100% with the other one
% relaxed moves the cursor at 50 cm/s. The .mat stores the map (info.W) as a
% matrix W, one row per axis and one column per muscle, so that for all
% samples at once:
%
%     velocity (N x 2) = gain * envelopes (N x 4) * W'
gain = info.gain_cm_s;                     % [cm/s]
W    = info.W;                             % 2 x 4
fprintf('\nPart 2 - MCI map\n');
fprintf('  x speed        = %g * (%s)\n', gain, info.map_x);
fprintf('  y speed (down) = %g * (%s)\n', gain, info.map_y);

% Dead zone: the MCI ignores any muscle below 5% of its max, so that
% resting noise does not move the cursor. Without this step the rebuilt
% cursor creeps away from the real one even when the subject is relaxed.
deadZone = 0.05;

% 1) Apply the dead zone to envOursD and to envMCID (values below deadZone become 0)
% 2) Compute the cursor velocity [cm/s] with the formula above:
%      velOurs  from your envelopes      (N x 2)
%      velMCI  from the MCI envelopes  (N x 2)

% YOUR CODE HERE
% 1) Apply the dead zone to envOursD and to envMCID (values below deadZone become 0)
% for any matrix A , doing A(A = 0)  indexes all elements of matrix A and checks if aij = 0 
 envOursD(envOursD < deadZone) = 0;
 envMCID(envMCID < deadZone) = 0;
% figure
% plot(t,envMCI(:,1),t,envOurs_norm(:,1));

% 2) Compute the cursor velocity [cm/s] with the formula : velocity (N x 2) = gain * envelopes (N x 4) * W':
%      velOurs  from your envelopes      (N x 2)
%      velMCI  from the MCI envelopes  (N x 2)
    velOurs = gain * envOursD * W';
    velMCI  = gain * envMCID *  W';

%% PART 3 - Rebuild the cursor from the envelopes
% -------------------------------------------------------------------------
% Velocity from the map (Part 2), then integrate (position = previous position +
% velocity * dt). Like the real cursor, ours can't leave the screen, so we
% clip it at the edges at every step.
%
% Integrating adds up every small error, so the rebuilt cursor slowly
% drifts away from the real one (try resetAtEachTarget = false to see it).
% With the reset, we restart from the logged position at each new target,
% so the error only builds up within one reach.
screenCm = info.screen_cm;                       % screen size [cm]
spawnT   = targets.spawn_time;                   % when each target appears [s]

% Rebuild two cursors, posOurs and posMCI [cm] (nD x 2, nD = numel(tD)),
% from velOurs and velMCI:
%   1) start both from the logged position: cursorCmD(1,:)
%   2) for each sample k = 2..nD: previous position + velocity * dt
%   3) clip to the screen: between [0 0] and screenCm
%   4) if resetAtEachTarget is true: at the first sample after each
%      spawnT, set both positions back to the logged one
%      Then compute, for both, the distance from the logged cursor at each
%      sample (errOurs, errMCI [cm]) and print its mean.

% YOUR CODE HERE

 posOurs(1,:) = cursorCmD(1,:);
 posMCI(1,:)  = cursorCmD(1,:);
 errOurs = zeros(nD,2);
 errMCI = zeros(nD,2);

 current_time= 0; % to calculate time 
 prev_time = 0;
 % loop from k= 2 : step = 1 : end = nD (number of samples) 
 for idx = 2 : 1 : nD 

    % 1) calculated positions from velocities success :D
    posOurs(idx,:) = dt*velOurs(idx,:) + posOurs(idx-1,:);
    posMCI (idx,:) = dt*velMCI (idx,:) + posMCI (idx-1,:);
    
    % 2) clip values to min [0 0] and max screenCm
    % implementation comments : screenCm is a row vector 1x2 , and at each
    % iteration of idx, our positions posOurs and posMCI are also row vectors
    % of 1x2 size. We want to bound the values between 0 and max screen
    % pixles. 
    % assuming our position is pos = [ 37 -1] and screenCm [30 20]. So we
    % want to stay bounded from [0 0] to screenCm
     % doing pos = minimum ( pos , screenCm) , returns [30 -1]
     % then post maximum ( pos,[0 0]) , returns [30 0] and we achieved our
     % goal!

    posOurs(idx,:) = min(posOurs(idx,:), screenCm);
    posOurs(idx,:) = max(posOurs(idx,:), [0 0]);
    
    posMCI(idx,:) = min(posMCI(idx,:), screenCm);
    posMCI(idx,:) = max(posMCI(idx,:), [0 0]);

    % % some old implementation ideas.for 3)
    % temp_ours =posOurs(idx,:); temp_MCI =posMCI(idx,:);
    % 
    % xposOurs= temp_ours(1,1); xposMCI = temp_MCI (1,1); yposOurs=
    % temp_ours(1,2); yposMCI = temp_MCI (1,2);
    % if (posOurs(idx,:) > screenCm ) ; posOurs(idx,:) = screenCm; end if
    % (posMCI(idx,:) >screenCm);posMCI(idx,:)= screenCm ; end
    
%   4) if resetAtEachTarget is true: at the first sample after each
%      spawnT, set both positions back to the logged one

        if resetAtEachTarget
           current_time = tD(idx);
           prev_time= tD(idx-1);

           % comments on implementation 
           % My goal is to reset the cursor positions to their recorded positions after each spawn
           % time. So we want to compare our elapsed time and the spawn
           % times. If the spawn time is more or equal to my elapsed time.
           % then we know a spawn time has passed.
           % the implementation in details 
           % the statement : spawnT(:,1)> prev_time , returns an array of size 48, since spawnT is 48x1. 
           % This array holds logical 1 or 0 values for each element of spawnT in the rows
           % compared with prev_time. We use statement
           % spawnT(:,1)<current_time. to see if any of the spawn time
           % stamps fall within our current and previous times. Finally we
           % use the & between statements to enforce our logic and then we
           % use any (). Because they are arrays of logical values if we
           % don't use any we get a true FOR THE WHOLE if statement if any
           % of the element wise & is true. But we want to analyse for each
           % element from 1 to 48 if our spawnT is in between the two
           % values previous and current time.
          
           if any(spawnT(:,1)>= prev_time & spawnT(:,1)<current_time) 
                 posOurs(idx,:)  = cursorCmD(idx,:);
                 posMCI (idx,:)  = cursorCmD(idx,:);
            end
        end

%      Then compute, for both, the distance from the logged cursor at each
%      sample (errOurs, errMCI [cm]) and print its mean.
            errOurs(idx,:) =  posOurs(idx,:)  - cursorCmD(idx,:);
            errMCI(idx,:) =  posMCI(idx,:)  - cursorCmD(idx,:);

           



 end
 errOurs_norm = vecnorm(errOurs, 2, 2); 
 errMCI_norm  = vecnorm(errMCI, 2, 2);
 mean_errOurs =  mean(errOurs_norm);
 mean_errMCI = mean(errMCI_norm);
fprintf('\nPart 3 - Cursor mean Errors\n');
fprintf('  Our error mean = %.2f cm\n', mean_errOurs);
fprintf('  MCI error mean = %.2f cm\n', mean_errMCI);

% Targets, converted to cm -> targCm (nTargets x 2)

% YOUR CODE HERE
%   cmPerPx  = info.cm_per_px;            % [cm per px]
%   targetPx  target positions [px], one row per target    (nTargets x 2)
    targCm = cmPerPx * targetPx ;

% --- Animation: the cursor moving on the screen, target by target (given)
% Black = the real (logged) cursor, red = the one we rebuilt from the EMG.
% The current target is the green circle. It plays in real time.
if animateCursor
    nD = numel(tD);
    % which target is on screen at each sample = the last one spawned
    curTarget = zeros(nD, 1);
    for s = 1:numel(spawnT)
        curTarget(tD >= spawnT(s)) = s;
    end

    figure('Name', 'Part 3 - cursor animation');
    rectangle('Position', [0 0 screenCm]); hold on;            % the screen
    hTarg = plot(nan, nan, 'go', 'MarkerSize', 25, 'LineWidth', 2);
    hLogTrail = plot(nan, nan, 'Color', [0.6 0.6 0.6]);
    hLog  = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    hOurs = plot(nan, nan, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    set(gca, 'YDir', 'reverse');                               % screen y grows downward
    axis equal; xlim([0 screenCm(1)]); ylim([0 screenCm(2)]);
    xlabel('x [cm]'); ylabel('y [cm]');
    legend([hLog hOurs hTarg], {'logged cursor', 'rebuilt cursor', 'target'}, 'Location', 'eastoutside');
    tic; k = 1;
    while k <= nD
        if curTarget(k) > 0
            set(hTarg, 'XData', targCm(curTarget(k),1), 'YData', targCm(curTarget(k),2));
        end
        k0 = max(1, k - 200);              % trail = last 2 s of the logged cursor
        set(hLogTrail, 'XData', cursorCmD(k0:k,1), 'YData', cursorCmD(k0:k,2));
        set(hLog,  'XData', cursorCmD(k,1), 'YData', cursorCmD(k,2));
        set(hOurs, 'XData', posOurs(k,1),   'YData', posOurs(k,2));
        title(sprintf('t = %.1f s   target %d', tD(k), curTarget(k)));
        drawnow;
        k = max(k + 1, round(toc * fsDown) + 1);   % jump to the sample of "now"
        pause(max(0, (k - 1) / fsDown - toc));     % ahead of the clock? wait
    end
end



%% PART 4 - Performance per reach
% -------------------------------------------------------------------------
% Each trial goes from the moment the target appears (spawn_time) to the
% moment it is hit (hit_time). In center-out, every other target is the
% center (the way back): we look at the outward reaches only. In the
% random task all targets count.
%   1) select the reaches to use: target hit (hit_time is not NaN) and
%      target not at the center (targetPx differs from centerPx)
%   2) for each of them, plot in order (a learning curve):
%      - the reach time [s]            (targets.reach_time)
%      - the overshoots [count]        (targets.overshoots)

% YOUR CODE HERE
% is the task centerout or random?
plot_x_axis_targets =  1: 1: size(targets,1);

reach_time =  targets.reach_time;
overshoots = targets.overshoots;

count=1;
valid_reach = []';
valid_overshoot = []';
validIndex = []';
x_px =centerPx(1);
y_px = centerPx(2);

% to identify valid targets for center out we want to exclude NaN values
% for reach time and to exclude targets at the center.
for n1=1: 1 : size(targets,1)
    if ( or(isnan(targets.hit_time(n1)), targets.target_x_px(n1) == x_px && targets.target_y_px(n1) == y_px))

    else
        % target valid
        validIndex(count) = n1;
        valid_reach(count)=reach_time(n1);
        valid_overshoot(count) = overshoots(n1);
        count = count+1;

    end
end


                               
figure('Name',"Reach Time and Overshoot");
subplot(2,1,1);
plot(validIndex, valid_reach, ...
    'Color', 'b', ...           % Line color
    'Marker', 'o', ...          % Circle marker
    'MarkerFaceColor', 'r', ... % Inside of the circle is red
    'MarkerEdgeColor', 'r');    % Border of the circle is redylabel('Reach Time [s]');
legend({'Reach Time'});
subtitle('Reach Time plotted with Target number');
ylabel('Time [s] ')
xlabel('Target')

subplot(2,1,2);
plot(validIndex, valid_overshoot, ...
    'Color', 'b', ...           % Line color
    'Marker', 'o', ...          % Circle marker
    'MarkerFaceColor', 'r', ... % Inside of the circle is red
    'MarkerEdgeColor', 'r');    % Border of the circle is red
subtitle('Overshoots plotted with Target number');

xlabel('Target')
legend({'Overshoots'});


    %% PART 5 - Robot reaching: drive a 2-link robot arm with your cursor
% -------------------------------------------------------------------------
% Picture a planar arm with two links in front of the screen: the shoulder
% is at the origin, the screen center is 0.40 m in front of it. The screen
% is shrunk to 0.40 m wide, so that every point of it is within reach of
% the arm, whatever the size of the monitor used in the recording. The
% rebuilt cursor tells the HAND where to go; the robot needs JOINT angles
% q = [q1 q2] (shoulder, elbow). The Jacobian J links the two speeds:
%
%       hand velocity = J(q) * joint velocity
%
% so at every step we invert it to get the joint velocity, and integrate it
% to get the angles (resolved-rate control):
%
%       joint velocity = J(q) \ hand velocity
%
% Formulas you need (s1 = sin(q1), c1 = cos(q1), s12 = sin(q1+q2), c12 = cos(q1+q2)):
%   forward kinematics (hand position):
%       x = L1*c1 + L2*c12,   y = L1*s1 + L2*s12
%   Jacobian:
%       J = [ -L1*s1 - L2*s12 , -L2*s12 ;
%              L1*c1 + L2*c12 ,  L2*c12 ]
%   inverse kinematics (only for the first point, "elbow down"):
%       q2 = acos((x^2 + y^2 - L1^2 - L2^2) / (2*L1*L2))
%       q1 = atan2(y, x) - atan2(L2*sin(q2), L1 + L2*cos(q2))

L1 = 0.30;  L2 = 0.30;                     % link lengths [m]
centerRobot = [0 0.40];                    % screen center in the robot frame [m]
scale = 0.40 / (screenCm(1) / 100);        % screen -> robot: the screen becomes 0.40 m wide
Kcorr = 10;                                % correction gain [1/s], see step 4
nD = numel(tD);

% 1) Screen frame -> robot frame. From posOurs [cm, y down] compute the
%    hand reference handRef [m, y up] (nD x 2):
%      - subtract the screen center in cm (centerPx * cmPerPx)
%      - convert cm -> m and multiply by scale
%      - flip y (screen y grows downward, robot y grows upward)
%      - add centerRobot
%    Do the same for the targets: targCm -> targRobot (nTargets x 2).

% YOUR CODE HERE
handRef = posOurs - (centerPx * cmPerPx);
handRef = (handRef ./100); % posOurs has a larger scale since we want to scale down the screen width from X approx = 59 cm to 40 cm we multiply by scale
handRef = handRef *scale; 
handRef(:,2) = -handRef(:,2);
handRef = handRef + centerRobot;

targRobot = targCm - (centerPx * cmPerPx);
targRobot = (targRobot ./100); 
targRobot = targRobot *scale; 
targRobot(:,2) = -targRobot(:,2);
targRobot = targRobot + centerRobot;

% 2) Hand velocity handVel [m/s] (nD x 2): derivative of handRef
%    (gradient of each column, divided by dt).

% YOUR CODE HERE
handVel = zeros(nD, 2);
handVel(:, 1) = gradient(handRef(:, 1), dt);
handVel(:, 2) = gradient(handRef(:, 2), dt);
% 3) Starting angles: inverse kinematics of the first point handRef(1,:)
%    (formulas above). Store them in q(1,:), with q = zeros(nD, 2).

% YOUR CODE HERE
%   inverse kinematics (only for the first point, "elbow down"):
%       q2 = acos((x^2 + y^2 - L1^2 - L2^2) / (2*L1*L2))
%       q1 = atan2(y, x) - atan2(L2*sin(q2), L1 + L2*cos(q2))
q = zeros(nD, 2);
x1 = handRef(1,1);
y1 = handRef(1,2);
q2 = acos((x1^2 + y1^2 - L1^2 - L2^2) / (2*L1*L2));
q1 = atan2(y1, x1) - atan2(L2*sin(q2), L1 + L2*cos(q2));
q(1,:) = [q1 q2]; 
% 4) For each sample k = 1..nD:
%      - compute the hand position with forward kinematics -> handArm(k,:)
%      - compute J at q(k,:) and its determinant -> detJ(k)
%      - if k < nD: joint velocity with the Jacobian, then integrate:
%            qdot = J \ (handVel(k,:) + Kcorr*(handRef(k,:) - handArm(k,:)))'
%            q(k+1,:) = q(k,:) + qdot' * dt
%    The Kcorr term pulls the hand back toward the reference: without it,
%    small integration errors (and the jumps when the cursor is reset at a
%    new target) would never be undone.


% YOUR CODE HERE
% Formulas you need (s1 = sin(q1), c1 = cos(q1), s12 = sin(q1+q2), c12 = cos(q1+q2)):
%   forward kinematics (hand position):
%       x = L1*c1 + L2*c12,   y = L1*s1 + L2*s12
%   Jacobian:
%       J = [ -L1*s1 - L2*s12 , -L2*s12 ;
%              L1*c1 + L2*c12 ,  L2*c12 ]

handArm = zeros(nD,2);

%initializing positions of the robot xr 'x robot'.
xr = zeros(nD, 1);
yr = zeros(nD, 1); 
detJ = zeros(nD, 1); 
for kr = 1:nD
    q1 = q(kr, 1);
    q2 = q(kr, 2);

    s1  = sin(q1);
    c1  = cos(q1);
    s12 = sin(q1 + q2);
    c12 = cos(q1 + q2);

    % Forward kinematics (hand position at step k)
    xr(kr) = L1 * c1 + L2 * c12;
    yr(kr) = L1 * s1 + L2 * s12;

    % Jacobian matrix at step k
    J = [ -L1 * s1 - L2 * s12, -L2 * s12 ;
        L1 * c1 + L2 * c12,  L2 * c12 ];
    handArm(kr,:) = [xr(kr) yr(kr)];
  

    detJ(kr) = det(J);
    if kr < nD
        qdot = J \ (handVel(kr,:).' + Kcorr * (handRef(kr,:) - handArm(kr,:)).');
        q(kr+1,:) = q(kr,:) + qdot.' * dt;
    end
end
% 5) Check and plot:
%      - distance between handArm and handRef [mm]: print its mean and max
%      - the smallest |detJ|: near 0 = singularity (arm fully stretched)
%      - figure: hand path with the targets, and the joint angles [deg] in time

trackErr = vecnorm(handArm - handRef, 2, 2) * 1000;
fprintf('\nPart 5 - Robot arm tracking\n');
fprintf('  Mean tracking error = %.2f mm\n', mean(trackErr));
fprintf('  Max  tracking error = %.2f mm\n', max(trackErr));
fprintf('  Min |det(J)| = %.4f\n', min(abs(detJ)));

figure('Name','Part 5 - hand path and joint angles');
subplot(1,2,1);
plot(handRef(:,1), handRef(:,2), 'r'); hold on;
plot(handArm(:,1), handArm(:,2), 'b');
plot(targRobot(:,1), targRobot(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
axis equal; xlabel('x [m]'); ylabel('y [m]');
legend('Reference','Arm hand','Targets');

subplot(1,2,2);
plot(tD, rad2deg(q(:,1)), tD, rad2deg(q(:,2)));
xlabel('Time [s]'); ylabel('Angle [deg]');
legend('q1 (shoulder)','q2 (elbow)');


% --- Animation of the arm (given): plays in real time
if animateArm
    figure('Name', 'Part 5 - robot arm animation');
    plot(targRobot(:,1), targRobot(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 2); hold on;
    hArm  = plot([0 0 0], [0 0 0], 'b-o', 'LineWidth', 3);
    hPath = plot(nan, nan, 'r');
    axis equal; axis([-0.35 0.35 -0.05 0.65]); xlabel('x [m]'); ylabel('y [m]');
    tic; k = 1;
    while k <= nD
        elbow = [L1*cos(q(k,1)), L1*sin(q(k,1))];
        set(hArm,  'XData', [0 elbow(1) handArm(k,1)], 'YData', [0 elbow(2) handArm(k,2)]);
        set(hPath, 'XData', handArm(1:k,1), 'YData', handArm(1:k,2));
        title(sprintf('t = %.1f s', tD(k)));
        drawnow;
        k = max(k + 1, round(toc * fsDown) + 1);   % jump to the sample of "now"
        pause(max(0, (k - 1) / fsDown - toc));     % ahead of the clock? wait
    end
end
