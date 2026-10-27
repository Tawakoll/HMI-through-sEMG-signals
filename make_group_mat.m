function make_group_mat(sessionFolder, systemName, outFile)
% MAKE_GROUP_MAT  Turn one BoMI session folder into a clean .mat for students.
%
%   make_group_mat(sessionFolder, systemName, outFile)
%
%   sessionFolder  the session folder created by the BoMI software, e.g.
%                  'C:\...\ID_01_Session_1' (it contains 'Acquisizioni'
%                  and 'Calibrazione')
%   systemName     'Cometa' or 'ZeroWire' (only written into info)
%   outFile        name of the .mat to create, e.g. 'Group_A.mat'
%
% The .mat contains these variables:
%   centerout          table, one row per sample (2000 per second), columns:
%                      t | raw_<muscle> x4 | env_<muscle> x4 | cursor_x_px,
%                      cursor_y_px | cursor_x_cm, cursor_y_cm
%   centerout_targets  table, one row per target: target_x_px, target_y_px,
%                      spawn_time, hit_time, reach_time, overshoots
%   random, random_targets   the same for the random task (if recorded)
%   info               struct: sampling frequency, muscles, screen size,
%                      calibration values and the BoMI map
%
% What it does:
%   - finds the center-out and random reaching files, their GameEvents
%     files and the calibration file inside the session folder
%   - keeps only the game part (state G) and rebuilds a regular time
%     vector starting at 0 (the Time column in the txt is uneven)
%   - puts event times on the same clock
%   - takes MAX and MEANRAW from the calibration file
%   - stores the BoMI map (gain and which muscle drives which direction)
%   - translates the muscle names to English
%   - leaves out all personal data (date of birth, etc.)
%
% Example:
%   make_group_mat('C:\Users\fedes\Downloads\ID_01_Cometa\ID_01\ID_01_Session_1', ...
%                  'Cometa', 'Group_A.mat');

fs = 2000;                                 % sampling frequency of the BoMI files [Hz]

%% Find the files
% The reaching files are in 'Acquisizioni' (or its subfolders) when you
% give the whole session folder; a folder with just the files works too.
% You can also give the center-out file itself: then only that task is read.
% (File names vary: 'Reaching_Internal_CenterOut' or 'Reaching_Internal
% CenterOut', so the search only looks for 'Reaching' ... 'CenterOut'.)
if exist(sessionFolder, 'file') == 2               % a file, not a folder
    coFile  = sessionFolder;
    rndFile = '';
    sessionFolder = fileparts(coFile);
else
    acq = fullfile(sessionFolder, 'Acquisizioni');
    if ~exist(acq, 'dir'), acq = sessionFolder; end
    coFile  = findOne(acq, '*Reaching*CenterOut*.txt', true);
    rndFile = findOne(acq, '*Reaching*Random*.txt', false);       % optional
end
calFile = '';                                                      % optional (older recordings need it)
if exist(fullfile(sessionFolder, 'Calibrazione'), 'dir')
    calFile = findOne(fullfile(sessionFolder, 'Calibrazione'), '*.txt', false);
end

%% Reaching files
co  = readTask(coFile, fs);
if isempty(rndFile)
    warning('No random reaching file found: the .mat will only have the center-out task.');
    rnd = [];
else
    rnd = readTask(rndFile, fs);
end

%% Muscles and calibration values
% Default muscle order of the BoMI. The calibration file, when there is
% one, names the muscles explicitly.
muscles      = {'R biceps', 'L biceps', 'L trapezius', 'R trapezius'};
calibMax     = [];
calibMeanRaw = nan(1, 4);
if ~isempty(calFile)
    txt = fileread(calFile);
    tok = regexp(txt, 'Muscles-->\s*([^\r\n]*)', 'tokens', 'once');
    if ~isempty(tok)
        musclesIta = strtrim(strsplit(tok{1}, ';'));
        musclesIta = musclesIta(~cellfun(@isempty, musclesIta));
        muscles = cellfun(@toEnglish, musclesIta, 'UniformOutput', false);
    end
    iA = strfind(txt, 'A seguire');
    if ~isempty(iA)
        lines = strsplit(txt(iA:end), {'\n', '\r'});
        lines = lines(~cellfun(@isempty, strtrim(lines)));
        calibMax     = sscanf(lines{2}, '%f')';    % row 1 after the text: MAX per muscle
        calibMeanRaw = sscanf(lines{4}, '%f')';    % row 3: MEANRAW per muscle (row 2 is all zeros)
    end
end
% Newer recordings write the MVC the BoMI used in the reaching file itself:
% that is exactly the value the envelopes were divided by, so it wins.
if ~isempty(co.mvc)
    calibMax = co.mvc;
end
if isempty(calibMax)
    error('No MVC found: give a session folder with a Calibrazione file, or a newer recording.');
end
co  = rmfield(co, 'mvc');
if ~isempty(rnd), rnd = rmfield(rnd, 'mvc'); end

%% BoMI map: muscles -> cursor velocity
% The BoMI is in velocity control:
%     cursor velocity [cm/s] = gain * envelopes * W'
% gain = Max Cursor Speed from the file header. W (2x4) says which muscle
% pushes which way: x = biceps pair (cols 1-2), y = trapezius pair (cols 3-4).
% Which side goes + and which goes - is a BoMI setting that is not written
% in the files (and the shoulders were swapped in the test sessions), so we
% read the direction from the recording: compare each muscle difference
% with how the cursor actually moved.
W = [1 -1 0 0;                             % x: muscle 1 minus muscle 2
     0  0 1 -1];                           % y: muscle 3 minus muscle 4
dt  = 1 / 100;                             % check at 100 Hz, plenty for this
k   = 1 : fs*dt : numel(co.t);
vel = [gradient(co.cursor_cm(k,1)) gradient(co.cursor_cm(k,2))] / dt;
for ax = 1:2
    drive = co.env_bomi(k,:) * W(ax,:)';   % muscle difference for this axis
    if sum(drive .* vel(:,ax)) < 0         % cursor moved the opposite way
        W(ax,:) = -W(ax,:);
    end
end

%% Info
info.system         = systemName;
info.fs             = fs;                  % [Hz]
info.muscles        = muscles;             % muscle order used in the column names
info.hardware_filter_hz = [10 500];        % band-pass applied by the EMG hardware
info.screen_px      = [2560 1600];
info.cm_per_px      = co.cm_per_px;
info.screen_cm      = (info.screen_px - 1) * info.cm_per_px;
info.center_px      = info.screen_px / 2;
info.calib_max      = calibMax;            % MVC: value of each muscle the BoMI calls 100% (1.0) [ADC counts]
info.calib_mean_raw = calibMeanRaw;        % mean raw value at rest (NaN if not recorded)
info.gain_cm_s      = co.max_speed;        % BoMI gain = max cursor speed [cm/s]
info.W              = W;                   % velocity = gain * env * W'
info.map_x          = mapText(W(1,:), muscles);
info.map_y          = mapText(W(2,:), muscles);
fprintf('BoMI map: x speed = %g * (%s), y speed (down) = %g * (%s)\n', ...
        info.gain_cm_s, info.map_x, info.gain_cm_s, info.map_y);

%% Tables and save
centerout         = samplesTable(co, muscles);
centerout_targets = targetsTable(co.events);
vars = {'info', 'centerout', 'centerout_targets'};
if ~isempty(rnd)
    random         = samplesTable(rnd, muscles); %#ok<NASGU>
    random_targets = targetsTable(rnd.events);   %#ok<NASGU>
    vars = [vars, {'random', 'random_targets'}];
end
save(outFile, vars{:}, '-v7');

fprintf('Saved %s: center-out %.0f s (%d targets)', outFile, co.t(end), numel(co.events.spawn_time));
if ~isempty(rnd)
    fprintf(', random %.0f s (%d targets)', rnd.t(end), numel(rnd.events.spawn_time));
end
fprintf('\n');
end


%% ------------------------------------------------------------------------
function T = samplesTable(task, muscles)
% One row per sample, one named column per signal.
names = strrep(muscles, ' ', '_');         % 'R biceps' -> 'R_biceps'
T = table(task.t, 'VariableNames', {'t'});
for m = 1:numel(muscles)
    T.(['raw_' names{m}]) = task.raw(:, m);
end
for m = 1:numel(muscles)
    T.(['env_' names{m}]) = task.env_bomi(:, m);
end
T.cursor_x_px = task.cursor_px(:, 1);
T.cursor_y_px = task.cursor_px(:, 2);
T.cursor_x_cm = task.cursor_cm(:, 1);
T.cursor_y_cm = task.cursor_cm(:, 2);
end


function T = targetsTable(ev)
% One row per target.
T = table(ev.target_px(:,1), ev.target_px(:,2), ev.spawn_time, ev.hit_time, ...
          ev.reach_time, ev.overshoots, 'VariableNames', ...
          {'target_x_px', 'target_y_px', 'spawn_time', 'hit_time', 'reach_time', 'overshoots'});
end


%% ------------------------------------------------------------------------
function task = readTask(file, fs)
% Read one reaching txt + its GameEvents file, keep only the game part.

% Header: read the lines up to the column names (the line starting with
% 'Time'). We keep only the max cursor speed and the MVC: the rest holds
% personal data. Older and newer recordings have different headers and
% columns, so columns are found by name.
fid = fopen(file, 'r');
header = '';
line = fgetl(fid);
while ischar(line) && isempty(regexp(line, '^Time\s', 'once'))
    header = [header line newline]; %#ok<AGROW>
    line = fgetl(fid);
end
names = strsplit(strtrim(line));                   % column names
iState = find(strcmp(names, 'State'));
fmt = repmat({'%f'}, 1, numel(names));  fmt{iState} = '%s';
C = textscan(fid, strjoin(fmt, ' '), 'MultipleDelimsAsOne', true);
fclose(fid);
col = @(name) C{strcmp(names, name)};

tok = regexp(header, 'Max Cursor Speed:\s*([\d.]+)', 'tokens', 'once');
task.max_speed = str2double(tok{1});

timeFile = col('Time');
isGame = strcmp(col('State'), 'G');
iGame = find(isGame);

n = numel(iGame);
task.t         = (0:n-1)' / fs;            % regular time: there are exactly fs samples per second
task.raw       = [col('RawM1') col('RawM2') col('RawM3') col('RawM4')];
task.raw       = task.raw(iGame, :);
task.env_bomi  = [col('%M1') col('%M2') col('%M3') col('%M4')];
task.env_bomi  = task.env_bomi(iGame, :);
task.cursor_px = [col('Mouse.X') col('Mouse.Y')];
task.cursor_px = task.cursor_px(iGame, :);
task.cursor_cm = [col('Mouse.X[cm]') col('Mouse.Y[cm]')];
task.cursor_cm = task.cursor_cm(iGame, :);

% MVC used by the BoMI (newer recordings only): one column per muscle.
% It should stay fixed; if the BoMI raised it during the run, warn.
task.mvc = [];
if any(strcmp(names, 'MvcM1'))
    mvc = [col('MvcM1') col('MvcM2') col('MvcM3') col('MvcM4')];
    mvc = mvc(iGame, :);
    if any(max(mvc) ~= min(mvc))
        warning('The MVC changed during %s: using the first value.', file);
    end
    task.mvc = mvc(1, :);
end

moving = task.cursor_px(:,1) > 0;
task.cm_per_px = median(task.cursor_cm(moving,1) ./ task.cursor_px(moving,1));

% Events: put them on the new clock by looking up the sample they fall on
evFile = strrep(file, '.txt', ' GameEvents.txt');
if ~exist(evFile, 'file')                  % same name with '_' instead of ' '
    evFile = strrep(file, '.txt', '_GameEvents.txt');
end
if ~exist(evFile, 'file')
    error('GameEvents file not found. Put it next to the reaching file, named:\n  %s', ...
          strrep(file, '.txt', ' GameEvents.txt'));
end
fid = fopen(evFile, 'r');
E = textscan(fid, '%f %s %f %s %s %s %s', 'HeaderLines', 1);
fclose(fid);
evTime = E{1}; evName = E{2}; evTrial = E{3};
tGameFile = timeFile(iGame);

isSpawn = strcmp(evName, 'TargetSpawned');
trials  = evTrial(isSpawn);
ev.target_px  = [str2double(E{4}(isSpawn)) str2double(E{5}(isSpawn))];
ev.spawn_time = nan(numel(trials), 1);
ev.hit_time   = nan(numel(trials), 1);
ev.overshoots = nan(numel(trials), 1);
for k = 1:numel(trials)
    ts = evTime(isSpawn & evTrial == trials(k));
    ev.spawn_time(k) = task.t(toSample(tGameFile, ts(1)));   % spawned before the game starts -> 0
    iHit = find(strcmp(evName, 'TargetHit') & evTrial == trials(k), 1);
    if ~isempty(iHit)                              % the last random target may never be hit
        ev.hit_time(k)   = task.t(toSample(tGameFile, evTime(iHit)));
        ev.overshoots(k) = str2double(E{7}{iHit});
    end
end
ev.reach_time = ev.hit_time - ev.spawn_time;
task.events = ev;
end


%% ------------------------------------------------------------------------
function f = findOne(folder, pattern, required)
% Find exactly one file matching pattern in folder or its subfolders
% (GameEvents and Commands files are skipped). If it isn't required and
% there is none, return ''.
found = {};
folders = {folder};
sub = dir(folder);
for s = 1:numel(sub)
    if sub(s).isdir && ~any(strcmp(sub(s).name, {'.', '..'}))
        folders{end+1} = fullfile(folder, sub(s).name); %#ok<AGROW>
    end
end
for k = 1:numel(folders)
    d = dir(fullfile(folders{k}, pattern));
    for j = 1:numel(d)
        if isempty(strfind(d(j).name, 'GameEvents')) && isempty(strfind(d(j).name, 'Commands'))
            found{end+1} = fullfile(folders{k}, d(j).name); %#ok<AGROW>
        end
    end
end
if isempty(found) && ~required
    f = '';
    return
end
if numel(found) ~= 1
    error('Expected one file matching %s in %s, found %d', pattern, folder, numel(found));
end
f = found{1};
end


function i = toSample(tFile, tEvent)
% First game sample at or after an event (file clock). Events before the
% game map to the first sample, events after it to the last one.
i = find(tFile >= tEvent, 1, 'first');
if isempty(i), i = numel(tFile); end
end


function s = mapText(w, muscles)
% e.g. w = [1 -1 0 0] -> 'R biceps - L biceps'
plus  = muscles{w == 1};
minus = muscles{w == -1};
s = [plus ' - ' minus];
end


function en = toEnglish(ita)
switch lower(strtrim(ita))
    case 'bicipite dx',  en = 'R biceps';
    case 'bicipite sx',  en = 'L biceps';
    case 'trapezio dx',  en = 'R trapezius';
    case 'trapezio sx',  en = 'L trapezius';
    otherwise,           en = ita;         % unknown name: keep it as it is
end
end
