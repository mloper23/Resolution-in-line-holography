%% resolution_target_demo_visible_gui.m
% Demo-style MATLAB code for resolution estimation with a polynomial-gap target.
% @ 06th July, 2026
% @ mloperaa, mloper23@eafit.edu.co; maria.josef.lopera.acosta@vub.be

% This version is intentionally GUI-simple and repository-friendly:
%   1) Background/DC removal with tunable Gaussian sigma.
%   2) Contrast enhancement with tunable imadjust limits.
%   3) Canny edge detection and Hough line detection.
%   4) Automatic center from ONLY TWO selected Hough lines, or manual center click.
%   5) User-defined maximum circular cross-section radius.
%   6) Circular profile visualization for each radius.
%   7) Adjacent-peak-pair detection ONLY: no non-adjacent peaks are compared.
%   8) Resolution from s(x)=A*x^3+B*x^2+C*x+1e-6.
%
% IMPORTANT GEOMETRY
%   Mag = L/z
%   p_x = sensor_pixel / Mag
%   x = radius_px * p_x
%


close all; clc;

%% ---------------- USER INPUTS ----------------
cfg = struct();

% Data location
cfg.name = '20251002_095542';
cfg.imagePath = fullfile('./data/100225/IC_L10/Phase', [cfg.name '_2.png']);
cfg.metadataPath = './data/100225/IC_L10/Resolution_recording.csv';

% Sensor pixel pitch [m]
cfg.sensor_pixel_m = 1.85e-6;

% Polynomial-gap coefficients, SI units
% s(x) = A*x^3 + B*x^2 + C*x + 1e-6, x and s in meters.
cfg.poly.A = 290682.4;    % [1/m^3]
cfg.poly.B = -100.5;      % [1/m^2]
cfg.poly.C = 0.0305;      % [1/m]
cfg.poly.offset = 1e-6;   % [m]

% Tunable preprocessing defaults
cfg.backgroundSigmaPx = 70;
cfg.contrastLow = 0.2;
cfg.contrastHigh = 1;
cfg.cannyLow = 0.04;
cfg.cannyHigh = 0.18;

% Hough defaults. The code can detect many lines, but it only SELECTS TWO
% lines for the center estimate and overlay.
cfg.houghThetaStepDeg = 0.5;
cfg.houghPeakNumber = 25;
cfg.houghThresholdFraction = 0.22;
cfg.houghNeighborhood = [31 31];
cfg.houghFillGapPx = 120;
cfg.houghMinLengthFraction = 0.12;
cfg.minAngleBetweenCenterLinesDeg = 8;

% Circular profiles
cfg.nAngles = 1440;
cfg.minRadiusPx = 20;
cfg.radiusStepPx = 2;
cfg.defaultInspectRadiusPx = 120;
cfg.defaultMaxRadiusPx = 350;
cfg.profileSmoothingSigma = 2;

% Adjacent peak-pair criteria
cfg.peakMinProminence = 0.2;
cfg.peakMinSeparationDeg = 0.1;
cfg.minPairSeparationUm = 0.25;
cfg.maxPairSeparationUm = 25;
cfg.useExpectedGapGate = true;
cfg.expectedGapLowFactor = 0.25;
cfg.expectedGapHighFactor = 4.0;
cfg.minValleyDepthFraction = 0.10;
cfg.requireConsecutiveResolved = 3;

% Demo behavior
cfg.livePlotDuringScan = true;
cfg.saveEachProfileDuringScan = false;  % set true if you want one PNG per radius
cfg.outputFolder = fullfile('./results_resolution_demo', cfg.name);

%% ---------------- LOAD DATA ----------------
imgRaw = im2double(imread(cfg.imagePath));
if ndims(imgRaw) == 3
    imgRaw = rgb2gray(imgRaw);
end
imgRaw = squeeze(imgRaw);

info = readtable(cfg.metadataPath, 'VariableNamingRule','preserve');
nameVar = findVar(info, ["name","Name","filename","Filename"]);
row = info(strcmp(string(info.(nameVar)), cfg.name), :);
assert(~isempty(row), 'Name %s not found in metadata table.', cfg.name);

L = getNum(row, "L (mm)");      % helper converts mm -> m
z = getNum(row, "z (mm)");      % helper converts mm -> m
Mag = L/z;
px = cfg.sensor_pixel_m / Mag;

fprintf('\nImage: %s\n', cfg.name);
fprintf('L = %.6g m, z = %.6g m, Mag = %.4f\n', L, z, Mag);
fprintf('Object-plane pixel size p_x = %.4f um\n\n', px*1e6);

if ~exist(cfg.outputFolder, 'dir'), mkdir(cfg.outputFolder); end
if ~exist(fullfile(cfg.outputFolder,'circular_profiles'), 'dir')
    mkdir(fullfile(cfg.outputFolder,'circular_profiles'));
end

%% ---------------- INITIAL STATE ----------------
state = struct();
state.imgRaw = imgRaw;
state.cleanPhase = [];
state.response = [];
state.enhanced = [];
state.edges = [];
state.allLines = struct([]);
state.centerLines = struct([]);
state.center = [size(imgRaw,2)/2, size(imgRaw,1)/2]; % [x,y] = [col,row]
state.centerMode = "initial image center";
state.maxRadiusPx = min(cfg.defaultMaxRadiusPx, floor(min(size(imgRaw))/2)-5);
state.inspectRadiusPx = min(cfg.defaultInspectRadiusPx, state.maxRadiusPx);
state.lastThetaDeg = [];
state.lastProfile = [];
state.lastResult = [];
state.scanTable = table();
state.allProfiles = [];
state.scanRadii = [];
state.scanThetaDeg = [];
state.resolvedRadiusPx = NaN;
state.resolutionUm = NaN;

state = runPreprocessing(state, cfg);
state.allLines = detectHoughLines(state.edges, cfg);
state.centerLines = selectTwoCenterLines(state.allLines, cfg);

%% ---------------- GUI ----------------
ui = createVisibleGUI(cfg, state, px);
setappdata(ui.fig, 'cfg', cfg);
setappdata(ui.fig, 'state', state);
setappdata(ui.fig, 'ui', ui);
setappdata(ui.fig, 'px', px);
redrawEverything(ui.fig);
setStatus(ui.fig, 'Ready. Tune preprocessing, then use Auto center or Manual center.');

%% ========================= CALLBACKS =========================
function onUpdateProcessing(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);

cfg = readControlsIntoCfg(cfg, ui);
state = runPreprocessing(state, cfg);
state.allLines = detectHoughLines(state.edges, cfg);
state.centerLines = selectTwoCenterLines(state.allLines, cfg);
state.lastProfile = [];
state.lastResult = [];

setAll(fig, cfg, state, ui);
redrawEverything(fig);
setStatus(fig, sprintf('Processing updated. Edge pixels = %d. Hough lines found = %d. Selected center lines = %d.', ...
    nnz(state.edges), numel(state.allLines), numel(state.centerLines)));
end

function onAutoCenter(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);

state.allLines = detectHoughLines(state.edges, cfg);
state.centerLines = selectTwoCenterLines(state.allLines, cfg);

if numel(state.centerLines) < 2
    setAll(fig, cfg, state, ui);
    redrawEverything(fig);
    setStatus(fig, 'Auto center failed: fewer than two usable Hough lines. Tune preprocessing or click the center manually.');
    return;
end

[center, ok] = intersectionOfTwoHoughLines(state.centerLines(1), state.centerLines(2));
if ~ok || any(~isfinite(center))
    setStatus(fig, 'Auto center failed: the two selected lines are nearly parallel. Click center manually.');
else
    state.center = center;
    state.centerMode = "automatic, two Hough lines";
    setStatus(fig, sprintf('Auto center set to x = %.2f px, y = %.2f px.', center(1), center(2)));
end

setAll(fig, cfg, state, ui);
redrawEverything(fig);
end

function onManualCenter(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);

axes(ui.axImage);
title(ui.axImage, 'Click the target center');
[xc, yc] = ginput(1);
state.center = [xc, yc];
state.centerMode = "manual click";

setAll(fig, cfg, state, ui);
redrawEverything(fig);
setStatus(fig, sprintf('Manual center set to x = %.2f px, y = %.2f px.', xc, yc));
end

function onClickMaxRadius(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);
px = getappdata(fig, 'px');

axes(ui.axImage);
title(ui.axImage, 'Click the maximum circular-section radius');
[xc, yc] = ginput(1);
state.maxRadiusPx = hypot(xc-state.center(1), yc-state.center(2));
state.inspectRadiusPx = min(state.inspectRadiusPx, state.maxRadiusPx);

set(ui.edMaxR, 'String', sprintf('%.1f', state.maxRadiusPx));
set(ui.edInspectR, 'String', sprintf('%.1f', state.inspectRadiusPx));
updateRadiusSlider(ui, state);

setAll(fig, cfg, state, ui);
redrawEverything(fig);
setStatus(fig, sprintf('Max radius = %.1f px = %.2f um in the object plane.', state.maxRadiusPx, state.maxRadiusPx*px*1e6));
end

function onPlotRadius(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);

state.maxRadiusPx = max(1, readNumber(ui.edMaxR, state.maxRadiusPx));
state.inspectRadiusPx = max(1, readNumber(ui.edInspectR, state.inspectRadiusPx));
state.inspectRadiusPx = min(state.inspectRadiusPx, state.maxRadiusPx);
set(ui.edInspectR, 'String', sprintf('%.1f', state.inspectRadiusPx));
updateRadiusSlider(ui, state);

setAll(fig, cfg, state, ui);
plotOneRadius(fig, state.inspectRadiusPx, false);
end

function onRadiusSlider(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig); %#ok<ASGLU>
state.inspectRadiusPx = get(src, 'Value');
set(ui.edInspectR, 'String', sprintf('%.1f', state.inspectRadiusPx));
setAll(fig, cfg, state, ui);
plotOneRadius(fig, state.inspectRadiusPx, false);
end

function onScanRadii(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);
px = getappdata(fig, 'px');

cfg = readControlsIntoCfg(cfg, ui);
state.maxRadiusPx = max(1, readNumber(ui.edMaxR, state.maxRadiusPx));
set(ui.edMaxR, 'String', sprintf('%.1f', state.maxRadiusPx));
updateRadiusSlider(ui, state);

radii = cfg.minRadiusPx:cfg.radiusStepPx:state.maxRadiusPx;
if isempty(radii)
    setStatus(fig, 'No radii to scan. Increase Max radius px.');
    return;
end

profileImage = state.enhanced;  % use the same processed image used for Canny/Hough
[theta, thetaDeg] = makeTheta(cfg.nAngles);
allProfiles = NaN(numel(radii), cfg.nAngles);

T = table('Size',[numel(radii), 14], ...
    'VariableTypes', {'double','double','double','logical','double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'radius_px','x_um','resolution_um','resolved','n_peaks','selected_peak_1_deg','selected_peak_2_deg', ...
                      'selected_pair_sep_um','expected_gap_um','valley_depth_fraction','valley_curvature','selection_score','n_adjacent_candidates','n_all_adjacent_pairs'});

for ii = 1:numel(radii)
    r = radii(ii);
    vals = circularProfile(profileImage, state.center, r, theta);
    vals = normalize01(vals);
    vals = imgaussfilt(vals, cfg.profileSmoothingSigma);
    allProfiles(ii,:) = vals;

    result = evaluateAdjacentPeakPair(thetaDeg, vals, r, px, cfg);
    x_m = radiusToTargetX(r, px, cfg);
    res_um = polynomialGap(x_m, cfg.poly)*1e6;

    T.radius_px(ii) = r;
    T.x_um(ii) = x_m*1e6;
    T.resolution_um(ii) = res_um;
    T.resolved(ii) = result.resolved;
    T.n_peaks(ii) = result.nPeaks;
    T.selected_peak_1_deg(ii) = result.selectedPeakAnglesDeg(1);
    T.selected_peak_2_deg(ii) = result.selectedPeakAnglesDeg(2);
    T.selected_pair_sep_um(ii) = result.selectedPairSepUm;
    T.expected_gap_um(ii) = result.expectedGapUm;
    T.valley_depth_fraction(ii) = result.valleyDepthFraction;
    T.valley_curvature(ii) = result.valleyCurvature;
    T.selection_score(ii) = result.selectionScore;
    T.n_adjacent_candidates(ii) = result.nAdjacentCandidates;
    T.n_all_adjacent_pairs(ii) = result.nAllAdjacentPairs;

    if cfg.livePlotDuringScan
        state.inspectRadiusPx = r;
        state.lastThetaDeg = thetaDeg;
        state.lastProfile = vals;
        state.lastResult = result;
        setAll(fig, cfg, state, ui);
        drawImagePanel(fig);
        drawProfilePanel(fig);
        setStatus(fig, sprintf('Scanning radius %.1f / %.1f px. Adjacent candidates: %d. Resolved: %d.', ...
            r, state.maxRadiusPx, result.nAdjacentCandidates, result.resolved));
        drawnow;
    end

    if cfg.saveEachProfileDuringScan
        saveSingleProfileFigure(thetaDeg, vals, result, r, px, fullfile(cfg.outputFolder,'circular_profiles'));
    end
end

idx = firstStableResolvedIndex(T.resolved, cfg.requireConsecutiveResolved);
if isnan(idx)
    state.resolvedRadiusPx = NaN;
    state.resolutionUm = NaN;
    msg = sprintf('Scan complete. No stable resolved radius found up to %.1f px.', state.maxRadiusPx);
else
    state.resolvedRadiusPx = T.radius_px(idx);
    state.resolutionUm = T.resolution_um(idx);
    msg = sprintf('Scan complete. First stable resolved radius = %.1f px. Resolution = %.3f um.', ...
        state.resolvedRadiusPx, state.resolutionUm);
end

state.scanTable = T;
state.allProfiles = allProfiles;
state.scanRadii = radii;
state.scanThetaDeg = thetaDeg;

setAll(fig, cfg, state, ui);
redrawEverything(fig);
setStatus(fig, msg);

writetable(T, fullfile(cfg.outputFolder, 'radius_scan_results.csv'));
writetable(makeSummaryTable(cfg, state, px), fullfile(cfg.outputFolder, 'resolution_summary.csv'));
end

function onSaveOutputs(src, ~)
fig = ancestor(src, 'figure');
[cfg, state, ui] = getAll(fig);
px = getappdata(fig, 'px');

if ~exist(cfg.outputFolder, 'dir'), mkdir(cfg.outputFolder); end
saveFigureCompat(ui.fig, fullfile(cfg.outputFolder, '00_gui_snapshot.png'));

if ~isempty(state.scanTable)
    writetable(state.scanTable, fullfile(cfg.outputFolder, 'radius_scan_results.csv'));
    writetable(makeSummaryTable(cfg, state, px), fullfile(cfg.outputFolder, 'resolution_summary.csv'));
    saveStackAndResolutionFigures(cfg, state, px);
end

setStatus(fig, sprintf('Saved outputs to %s', cfg.outputFolder));
end

%% ========================= GUI CREATION =========================
function ui = createVisibleGUI(cfg, state, px)
% Compact GUI: all important buttons stay visible on a laptop screen.
% Controls are split into two columns on the right. No scrolling is needed.
ui = struct();
ui.fig = figure('Name','Resolution target demo - compact controls', ...
    'Color','w', 'Units','pixels', 'Position',[20 40 1520 760], ...
    'MenuBar','none', 'Toolbar','figure', 'NumberTitle','off', ...
    'KeyPressFcn', @onKeyPress);

% Four diagnostic axes. Keep these in fixed positions and leave the right
% side for controls.
ui.axImage   = axes('Parent',ui.fig, 'Units','pixels', 'Position',[40 420 500 300]);
ui.axEdges   = axes('Parent',ui.fig, 'Units','pixels', 'Position',[570 420 500 300]);
ui.axProfile = axes('Parent',ui.fig, 'Units','pixels', 'Position',[40 65 500 285]);
ui.axStack   = axes('Parent',ui.fig, 'Units','pixels', 'Position',[570 65 500 285]);

% Right-side compact control panel.
xPanel = 1100;
wPanel = 390;
yTop = 720;
addTitle(ui.fig, 'Controls', xPanel, yTop, wPanel, 22);

% --- Big action buttons: always visible ---
y = yTop - 38;
ui.btnUpdate = addButton(ui.fig, 'Update processing', xPanel, y, 185, 30, @onUpdateProcessing);
ui.btnAuto   = addButton(ui.fig, 'Auto center',        xPanel+200, y, 185, 30, @onAutoCenter); y = y-40;
ui.btnManual = addButton(ui.fig, 'Manual center click', xPanel, y, 185, 30, @onManualCenter);
ui.btnPlot   = addButton(ui.fig, 'Plot selected radius',xPanel+200, y, 185, 30, @onPlotRadius); y = y-40;
ui.btnScan   = addButton(ui.fig, 'SCAN RADII',          xPanel, y, 185, 34, @onScanRadii);
ui.btnSave   = addButton(ui.fig, 'Save outputs',        xPanel+200, y, 185, 34, @onSaveOutputs); y = y-50;

% --- Two-column parameters ---
x1 = xPanel;
x2 = xPanel + 200;
w  = 185;
y1 = y;
y2 = y;
rowGap = 43;

addTitle(ui.fig, 'Preprocessing', x1, y1, w, 20); y1 = y1-28;
ui.edSigma  = addCompactEdit(ui.fig, 'Background sigma', cfg.backgroundSigmaPx, x1, y1, w); y1 = y1-rowGap;
ui.edCLow   = addCompactEdit(ui.fig, 'Contrast low', cfg.contrastLow, x1, y1, w); y1 = y1-rowGap;
ui.edCHigh  = addCompactEdit(ui.fig, 'Contrast high', cfg.contrastHigh, x1, y1, w); y1 = y1-rowGap;
ui.edCanLow = addCompactEdit(ui.fig, 'Canny low', cfg.cannyLow, x1, y1, w); y1 = y1-rowGap;
ui.edCanHigh= addCompactEdit(ui.fig, 'Canny high', cfg.cannyHigh, x1, y1, w); y1 = y1-rowGap;
ui.edFillGap= addCompactEdit(ui.fig, 'Hough fill gap', cfg.houghFillGapPx, x1, y1, w); y1 = y1-rowGap;
ui.edMinLen = addCompactEdit(ui.fig, 'Hough min length', cfg.houghMinLengthFraction, x1, y1, w); y1 = y1-rowGap;

addTitle(ui.fig, 'Circular scan', x2, y2, w, 20); y2 = y2-28;
ui.edMaxR = addCompactEdit(ui.fig, 'Max radius px', state.maxRadiusPx, x2, y2, w); y2 = y2-rowGap;
ui.btnMaxR = addButton(ui.fig, 'Click max radius', x2, y2+8, w, 28, @onClickMaxRadius); y2 = y2-rowGap;
ui.edInspectR = addCompactEdit(ui.fig, 'Inspect radius px', state.inspectRadiusPx, x2, y2, w); y2 = y2-rowGap;
ui.sliderR = uicontrol(ui.fig, 'Style','slider', 'Units','pixels', ...
    'Position',[x2 y2+12 w 20], 'Min',1, 'Max',state.maxRadiusPx, 'Value',state.inspectRadiusPx, ...
    'Callback',@onRadiusSlider); y2 = y2-rowGap;
ui.edPeakProm = addCompactEdit(ui.fig, 'Peak prominence', cfg.peakMinProminence, x2, y2, w); y2 = y2-rowGap;
ui.edMaxPair = addCompactEdit(ui.fig, 'Max pair sep um', cfg.maxPairSeparationUm, x2, y2, w); y2 = y2-rowGap;
ui.edValley = addCompactEdit(ui.fig, 'Min valley depth', cfg.minValleyDepthFraction, x2, y2, w); y2 = y2-rowGap;

% Put checkboxes below both columns, still above status.
yCheck = min(y1, y2) + 10;
ui.chkExpected = uicontrol(ui.fig, 'Style','checkbox', 'Units','pixels', ...
    'Position',[xPanel yCheck+30 185 22], 'String','Expected-gap gate', ...
    'Value',cfg.useExpectedGapGate, 'BackgroundColor','w');
ui.chkLive = uicontrol(ui.fig, 'Style','checkbox', 'Units','pixels', ...
    'Position',[xPanel+200 yCheck+30 185 22], 'String','Live scan plot', ...
    'Value',cfg.livePlotDuringScan, 'BackgroundColor','w');

uicontrol(ui.fig, 'Style','text', 'Units','pixels', 'Position',[xPanel yCheck-28 wPanel 38], ...
    'String',sprintf('p_x = %.4f um     x = radius_px p_x\nKeys: A auto | M manual | P plot | S scan | W save', px*1e6), ...
    'HorizontalAlignment','left', 'BackgroundColor','w', 'FontWeight','bold');

ui.txtStatus = uicontrol(ui.fig, 'Style','text', 'Units','pixels', ...
    'Position',[xPanel 25 wPanel 115], 'String','Status:', 'HorizontalAlignment','left', ...
    'BackgroundColor',[0.95 0.95 0.95]);
end

function hEdit = addCompactEdit(parent, label, value, x, y, w)
% Compact label/edit pair for the two-column control panel.
uicontrol(parent, 'Style','text', 'Units','pixels', 'Position',[x y+24 w 16], ...
    'String',label, 'BackgroundColor','w', 'HorizontalAlignment','left', 'FontSize',9);
hEdit = uicontrol(parent, 'Style','edit', 'Units','pixels', 'Position',[x y w 24], ...
    'String',num2str(value), 'BackgroundColor','w');
end

function onKeyPress(src, event)
% Keyboard shortcuts in case MATLAB clips the right-side panel on small screens.
ui = getappdata(src, 'ui');
if isempty(ui), return; end
switch lower(event.Key)
    case 'a'
        onAutoCenter(ui.btnAuto, []);
    case 'm'
        onManualCenter(ui.btnManual, []);
    case 'p'
        onPlotRadius(ui.btnPlot, []);
    case 's'
        onScanRadii(ui.btnScan, []);
    case 'w'
        onSaveOutputs(ui.btnSave, []);
end
end

function addTitle(parent, str, x, y, w, h)
uicontrol(parent, 'Style','text', 'Units','pixels', 'Position',[x y w h], ...
    'String',str, 'BackgroundColor','w', 'HorizontalAlignment','left', 'FontWeight','bold', 'FontSize',11);
end

function hEdit = addLabeledEdit(parent, label, value, x, y, w)
uicontrol(parent, 'Style','text', 'Units','pixels', 'Position',[x y+1 145 20], ...
    'String',label, 'BackgroundColor','w', 'HorizontalAlignment','left');
hEdit = uicontrol(parent, 'Style','edit', 'Units','pixels', 'Position',[x+150 y 80 23], ...
    'String',num2str(value), 'BackgroundColor','w');
end

function hBtn = addButton(parent, label, x, y, w, h, cb)
hBtn = uicontrol(parent, 'Style','pushbutton', 'Units','pixels', 'Position',[x y w h], ...
    'String',label, 'Callback',cb, 'FontWeight','bold');
end

function updateRadiusSlider(ui, state)
mx = max(1, state.maxRadiusPx);
val = min(max(1, state.inspectRadiusPx), mx);
set(ui.sliderR, 'Min',1, 'Max',mx, 'Value',val);
if mx > 1
    set(ui.sliderR, 'SliderStep', [min(1,1/(mx-1)), min(1,20/(mx-1))]);
end
end

%% ========================= DRAWING =========================
function redrawEverything(fig)
drawImagePanel(fig);
drawEdgesPanel(fig);
drawProfilePanel(fig);
drawStackPanel(fig);
end

function drawImagePanel(fig)
[cfg, state, ui] = getAll(fig); %#ok<ASGLU>
px = getappdata(fig, 'px');
axes(ui.axImage); cla(ui.axImage);
imshow(state.enhanced, [], 'Parent', ui.axImage); hold(ui.axImage, 'on');
plotSelectedLines(ui.axImage, state.centerLines, 2.0);
plot(ui.axImage, state.center(1), state.center(2), '+', 'MarkerSize', 13, 'LineWidth', 2);
drawCircle(ui.axImage, state.center, state.maxRadiusPx, '--', 1.2);
drawCircle(ui.axImage, state.center, state.inspectRadiusPx, '-', 1.3);
title(ui.axImage, sprintf('Processed image | center: %s', state.centerMode), 'Interpreter','none');
text(ui.axImage, 8, 18, sprintf('x = r p_x, p_x = %.4f um', px*1e6), ...
    'Color','w','BackgroundColor','k','Margin',2, 'FontWeight','bold');
hold(ui.axImage, 'off'); axis(ui.axImage, 'image');
end

function drawEdgesPanel(fig)
[cfg, state, ui] = getAll(fig); %#ok<ASGLU>
axes(ui.axEdges); cla(ui.axEdges);
imshow(state.edges, [], 'Parent', ui.axEdges); hold(ui.axEdges, 'on');
plotSelectedLines(ui.axEdges, state.centerLines, 2.0);
plot(ui.axEdges, state.center(1), state.center(2), '+', 'MarkerSize', 13, 'LineWidth', 2);
title(ui.axEdges, sprintf('Canny edges + ONLY TWO selected Hough lines (%d found)', numel(state.allLines)));
hold(ui.axEdges, 'off'); axis(ui.axEdges, 'image');
end

function drawProfilePanel(fig)
[cfg, state, ui] = getAll(fig); %#ok<ASGLU>
axes(ui.axProfile); cla(ui.axProfile);
if isempty(state.lastProfile)
    text(ui.axProfile, 0.08, 0.5, 'Press Plot selected radius or Scan radii.', 'Units','normalized');
    axis(ui.axProfile, 'off');
    return;
end

plot(ui.axProfile, state.lastThetaDeg, state.lastProfile, 'LineWidth',1.2); hold(ui.axProfile, 'on');
res = state.lastResult;
if ~isempty(res) && res.nPeaks > 0
    plot(ui.axProfile, res.peakAnglesDeg, res.peakValues, 'o', 'MarkerSize', 4);
end
if ~isempty(res) && all(isfinite(res.selectedPeakAnglesDeg))
    xline(ui.axProfile, res.selectedPeakAnglesDeg(1), '--', 'LineWidth', 1.2);
    xline(ui.axProfile, res.selectedPeakAnglesDeg(2), '--', 'LineWidth', 1.2);
    ttl = sprintf('r = %.1f px | selected ADJACENT pair sep = %.2f um | expected = %.2f um | resolved = %d', ...
        state.inspectRadiusPx, res.selectedPairSepUm, res.expectedGapUm, res.resolved);
else
    ttl = sprintf('r = %.1f px | no adjacent pair selected', state.inspectRadiusPx);
end
xlabel(ui.axProfile, 'Angle [deg]'); ylabel(ui.axProfile, 'Normalized profile');
title(ui.axProfile, ttl); grid(ui.axProfile, 'on'); xlim(ui.axProfile, [0 360]);
hold(ui.axProfile, 'off');
end

function drawStackPanel(fig)
[cfg, state, ui] = getAll(fig); %#ok<ASGLU>
px = getappdata(fig, 'px');
axes(ui.axStack); cla(ui.axStack);
if isempty(state.allProfiles)
    text(ui.axStack, 0.06, 0.55, 'Profile stack and resolution curve appear after Scan radii.', 'Units','normalized');
    axis(ui.axStack, 'off');
    return;
end
imagesc(ui.axStack, state.scanThetaDeg, state.scanRadii*px*1e6, state.allProfiles);
axis(ui.axStack, 'xy'); colormap(ui.axStack, gray); hold(ui.axStack, 'on');
if ~isnan(state.resolvedRadiusPx)
    yline(ui.axStack, state.resolvedRadiusPx*px*1e6, '--', sprintf('resolution %.2f um', state.resolutionUm));
end
xlabel(ui.axStack, 'Angle [deg]'); ylabel(ui.axStack, 'x = r p_x [um]');
title(ui.axStack, 'Circular profile stack');
hold(ui.axStack, 'off');
end

function plotSelectedLines(ax, lines, lw)
for k = 1:min(2,numel(lines))
    xy = [lines(k).point1; lines(k).point2];
    plot(ax, xy(:,1), xy(:,2), '-', 'LineWidth', lw);
end
end

function drawCircle(ax, center, radius, style, lw)
if radius <= 0 || any(~isfinite(center)), return; end
th = linspace(0,2*pi,600);
plot(ax, center(1)+radius*cos(th), center(2)+radius*sin(th), style, 'LineWidth', lw);
end

function plotOneRadius(fig, radiusPx, saveFlag)
[cfg, state, ui] = getAll(fig);
px = getappdata(fig, 'px');
[theta, thetaDeg] = makeTheta(cfg.nAngles);
vals = circularProfile(state.enhanced, state.center, radiusPx, theta);
vals = normalize01(vals);
vals = imgaussfilt(vals, cfg.profileSmoothingSigma);
res = evaluateAdjacentPeakPair(thetaDeg, vals, radiusPx, px, cfg);

state.inspectRadiusPx = radiusPx;
state.lastThetaDeg = thetaDeg;
state.lastProfile = vals;
state.lastResult = res;
setAll(fig, cfg, state, ui);
redrawEverything(fig);

if saveFlag
    saveSingleProfileFigure(thetaDeg, vals, res, radiusPx, px, fullfile(cfg.outputFolder,'circular_profiles'));
end
setStatus(fig, sprintf('Profile plotted at r = %.1f px = %.2f um. Peaks = %d. Adjacent candidates = %d. Resolved = %d.', ...
    radiusPx, radiusPx*px*1e6, res.nPeaks, res.nAdjacentCandidates, res.resolved));
end

%% ========================= PROCESSING =========================
function state = runPreprocessing(state, cfg)
img = state.imgRaw;
background = imgaussfilt(img, cfg.backgroundSigmaPx);
cleanPhase = img - background;
response = abs(cleanPhase);
response = normalize01(response);

lo = max(0, min(0.99, cfg.contrastLow));
hi = max(lo + eps, min(1, cfg.contrastHigh));
enhanced = imadjust(response, [lo hi], [0 1]);

cl = max(0, min(0.99, cfg.cannyLow));
ch = max(cl + eps, min(1, cfg.cannyHigh));
edges = edge(enhanced, 'canny', [cl ch]);

state.cleanPhase = cleanPhase;
state.response = response;
state.enhanced = enhanced;
state.edges = edges;
end

function lines = detectHoughLines(BW, cfg)
if nnz(BW) < 10
    lines = struct([]); return;
end
[H,T,R] = hough(BW, 'Theta', -90:cfg.houghThetaStepDeg:(90-cfg.houghThetaStepDeg));
if isempty(H) || max(H(:)) == 0
    lines = struct([]); return;
end
P = houghpeaks(H, cfg.houghPeakNumber, ...
    'Threshold', cfg.houghThresholdFraction*max(H(:)), ...
    'NHoodSize', cfg.houghNeighborhood);
minLength = max(5, round(cfg.houghMinLengthFraction * min(size(BW))));
lines = houghlines(BW, T, R, P, 'FillGap', cfg.houghFillGapPx, 'MinLength', minLength);
if isempty(lines), return; end
len = arrayfun(@(s) norm(s.point1 - s.point2), lines);
[~, idx] = sort(len, 'descend');
lines = lines(idx);
end

function selected = selectTwoCenterLines(lines, cfg)
selected = struct([]);
if isempty(lines), return; end
selected = lines(1);
if numel(lines) == 1, return; end
for k = 2:numel(lines)
    dth = abs(wrapTo180Local(lines(k).theta - selected(1).theta));
    dth = min(dth, 180-dth);
    if dth >= cfg.minAngleBetweenCenterLinesDeg
        selected(2) = lines(k); %#ok<AGROW>
        return;
    end
end
% Fallback: use the second longest line even if angle separation is poor.
selected(2) = lines(2);
end

function [center, ok] = intersectionOfTwoHoughLines(line1, line2)
% Hough convention: x*cos(theta) + y*sin(theta) = rho
A = [cosd(line1.theta), sind(line1.theta); cosd(line2.theta), sind(line2.theta)];
b = [line1.rho; line2.rho];
if rcond(A) < 1e-8
    center = [NaN NaN]; ok = false; return;
end
xy = A\b;
center = xy(:).';
ok = all(isfinite(center));
end

function vals = circularProfile(img, center, radiusPx, theta)
x = center(1) + radiusPx*cos(theta);
y = center(2) + radiusPx*sin(theta);
vals = interp2(img, x, y, 'linear', NaN);
if any(isnan(vals))
    medv = median(vals(~isnan(vals)));
    if isempty(medv) || isnan(medv), medv = 0; end
    vals(isnan(vals)) = medv;
end
end

function result = evaluateAdjacentPeakPair(thetaDeg, vals, radiusPx, px, cfg)
% This function compares ONLY adjacent peaks in angular order.
% Non-adjacent peaks are never paired.

vals = vals(:).';
thetaDeg = thetaDeg(:).';
N = numel(vals);

result = struct();
result.resolved = false;
result.nPeaks = 0;
result.peakAnglesDeg = [];
result.peakValues = [];
result.selectedPeakAnglesDeg = [NaN NaN];
result.selectedPeakValues = [NaN NaN];
result.selectedPairSepUm = NaN;
result.expectedGapUm = polynomialGap(radiusToTargetX(radiusPx, px, cfg), cfg.poly)*1e6;
result.valleyDepthFraction = NaN;
result.valleyCurvature = NaN;
result.selectionScore = NaN;
result.nAdjacentCandidates = 0;
result.nAllAdjacentPairs = 0;

minSepSamples = max(1, round(cfg.peakMinSeparationDeg / mean(diff(thetaDeg))));
[TF, promAll] = islocalmax(vals, 'MinProminence', cfg.peakMinProminence, 'MinSeparation', minSepSamples);
locs = find(TF);
pks = vals(locs);
prom = promAll(locs);

% Guarantee angular sorting.
[locs, order] = sort(locs);
pks = pks(order);
prom = prom(order);

result.nPeaks = numel(locs);
result.peakAnglesDeg = thetaDeg(locs);
result.peakValues = pks;

if numel(locs) < 2
    return;
end

locsExt = [locs, locs(1)+N];
pksExt = [pks, pks(1)];
promExt = [prom, prom(1)];

bestScore = -Inf;
expected = max(result.expectedGapUm, eps);

for ii = 1:numel(locs)
    i1 = locsExt(ii);
    i2 = locsExt(ii+1);
    dSamples = i2 - i1;
    result.nAllAdjacentPairs = result.nAllAdjacentPairs + 1;

    % Physical separation along the circular cross-section.
    dTheta = dSamples * 2*pi / N;
    sepUm = radiusPx * px * dTheta * 1e6;

    if sepUm < cfg.minPairSeparationUm || sepUm > cfg.maxPairSeparationUm
        continue;
    end

    if cfg.useExpectedGapGate
        if sepUm < cfg.expectedGapLowFactor*expected || sepUm > cfg.expectedGapHighFactor*expected
            continue;
        end
    end

    idxSeg = mod((i1:i2)-1, N) + 1;
    seg = vals(idxSeg);
    [valleyValue, valleyLocalIdx] = min(seg);
    valleyIdx = idxSeg(valleyLocalIdx);

    lowerPeak = min(pksExt(ii), pksExt(ii+1));
    baseline = median(vals);
    valleyDepth = (lowerPeak - valleyValue) / max(eps, lowerPeak - baseline);

    idxLeft = mod(valleyIdx-2, N) + 1;
    idxRight = mod(valleyIdx, N) + 1;
    curvature = vals(idxLeft) - 2*vals(valleyIdx) + vals(idxRight);

    % Selection score: adjacent pair whose physical separation is closest to the
    % expected polynomial gap, while still favoring strong/prominent peaks and a valley.
    closenessPenalty = abs(log(sepUm / expected));
    score = -closenessPenalty + 0.25*min(promExt(ii), promExt(ii+1)) + 0.25*valleyDepth;

    result.nAdjacentCandidates = result.nAdjacentCandidates + 1;
    if score > bestScore
        bestScore = score;
        result.selectedPairSepUm = sepUm;
        result.selectedPeakAnglesDeg = [thetaDeg(mod(locsExt(ii)-1,N)+1), thetaDeg(mod(locsExt(ii+1)-1,N)+1)];
        result.selectedPeakValues = [pksExt(ii), pksExt(ii+1)];
        result.valleyDepthFraction = valleyDepth;
        result.valleyCurvature = curvature;
        result.selectionScore = score;
        result.resolved = valleyDepth >= cfg.minValleyDepthFraction && curvature > 0;
    end
end
end

function idx = firstStableResolvedIndex(resolved, nConsecutive)
idx = NaN;
resolved = logical(resolved(:));
if numel(resolved) < nConsecutive, return; end
for k = 1:(numel(resolved)-nConsecutive+1)
    if all(resolved(k:k+nConsecutive-1))
        idx = k;
        return;
    end
end
end

%% ========================= GEOMETRY AND EQUATION =========================
function x_m = radiusToTargetX(radiusPx, px, cfg) %#ok<INUSD>
% Current geometry assumption:
% x is the object-plane radial coordinate of the circular cross-section.
% p_x = sensor_pixel / Mag, therefore x = radius_px * p_x.
x_m = radiusPx .* px;
end

function s_m = polynomialGap(x_m, poly)
s_m = poly.A.*x_m.^3 + poly.B.*x_m.^2 + poly.C.*x_m + poly.offset;
end

function [theta, thetaDeg] = makeTheta(nAngles)
theta = linspace(0, 2*pi, nAngles+1);
theta(end) = [];
thetaDeg = rad2deg(theta);
end

%% ========================= SAVE FUNCTIONS =========================
function saveSingleProfileFigure(thetaDeg, vals, res, radiusPx, px, outFolder)
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
f = figure('Visible','off','Color','w','Position',[100 100 900 320]);
plot(thetaDeg, vals, 'LineWidth',1.2); hold on;
if res.nPeaks > 0
    plot(res.peakAnglesDeg, res.peakValues, 'o', 'MarkerSize',4);
end
if all(isfinite(res.selectedPeakAnglesDeg))
    xline(res.selectedPeakAnglesDeg(1), '--');
    xline(res.selectedPeakAnglesDeg(2), '--');
end
xlabel('Angle [deg]'); ylabel('Normalized profile'); grid on; xlim([0 360]);
title(sprintf('r = %.1f px, x = %.2f um, selected adjacent sep = %.2f um, resolved = %d', ...
    radiusPx, radiusPx*px*1e6, res.selectedPairSepUm, res.resolved));
saveFigureCompat(f, fullfile(outFolder, sprintf('profile_r_%06.1f_px.png', radiusPx)));
close(f);
end

function saveStackAndResolutionFigures(cfg, state, px)
if isempty(state.allProfiles), return; end

f1 = figure('Visible','off','Color','w','Position',[100 100 900 500]);
imagesc(state.scanThetaDeg, state.scanRadii*px*1e6, state.allProfiles); axis xy; colormap gray;
xlabel('Angle [deg]'); ylabel('x = r p_x [um]'); title('Circular profile stack'); hold on;
if ~isnan(state.resolvedRadiusPx)
    yline(state.resolvedRadiusPx*px*1e6, '--', sprintf('resolution %.2f um', state.resolutionUm));
end
saveFigureCompat(f1, fullfile(cfg.outputFolder, 'profile_stack.png'));
close(f1);

if isempty(state.scanTable), return; end
T = state.scanTable;
f2 = figure('Visible','off','Color','w','Position',[100 100 750 450]);
plot(T.x_um, T.resolution_um, 'LineWidth',1.3); hold on; grid on;
scatter(T.x_um(T.resolved), T.resolution_um(T.resolved), 25, 'filled');
if ~isnan(state.resolutionUm)
    xr = state.resolvedRadiusPx*px*1e6;
    plot(xr, state.resolutionUm, 'p', 'MarkerSize',12, 'LineWidth',1.5);
    text(xr, state.resolutionUm, sprintf('  %.2f um', state.resolutionUm));
end
xlabel('x = r p_x [um]'); ylabel('s(x) [um]'); title('Resolution from calibrated polynomial gap');
saveFigureCompat(f2, fullfile(cfg.outputFolder, 'resolution_curve.png'));
close(f2);
end

function summary = makeSummaryTable(cfg, state, px)
summary = table();
summary.name = string(cfg.name);
summary.center_x_px = state.center(1);
summary.center_y_px = state.center(2);
summary.center_mode = string(state.centerMode);
summary.object_pixel_size_um = px*1e6;
summary.max_radius_px = state.maxRadiusPx;
summary.max_x_um = state.maxRadiusPx*px*1e6;
summary.resolved_radius_px = state.resolvedRadiusPx;
summary.resolved_x_um = state.resolvedRadiusPx*px*1e6;
summary.resolution_um = state.resolutionUm;
summary.background_sigma_px = cfg.backgroundSigmaPx;
summary.contrast_low = cfg.contrastLow;
summary.contrast_high = cfg.contrastHigh;
summary.canny_low = cfg.cannyLow;
summary.canny_high = cfg.cannyHigh;
summary.peak_min_prominence = cfg.peakMinProminence;
summary.max_pair_separation_um = cfg.maxPairSeparationUm;
summary.min_valley_depth_fraction = cfg.minValleyDepthFraction;
summary.require_consecutive_resolved = cfg.requireConsecutiveResolved;
end

function saveFigureCompat(fig, pathOut)
try
    exportgraphics(fig, pathOut, 'Resolution', 200);
catch
    saveas(fig, pathOut);
end
end

%% ========================= UTILITIES =========================
function cfg = readControlsIntoCfg(cfg, ui)
cfg.backgroundSigmaPx = max(1, readNumber(ui.edSigma, cfg.backgroundSigmaPx));
cfg.contrastLow = readNumber(ui.edCLow, cfg.contrastLow);
cfg.contrastHigh = readNumber(ui.edCHigh, cfg.contrastHigh);
cfg.cannyLow = readNumber(ui.edCanLow, cfg.cannyLow);
cfg.cannyHigh = readNumber(ui.edCanHigh, cfg.cannyHigh);
cfg.houghFillGapPx = max(1, round(readNumber(ui.edFillGap, cfg.houghFillGapPx)));
cfg.houghMinLengthFraction = max(0.01, readNumber(ui.edMinLen, cfg.houghMinLengthFraction));
cfg.peakMinProminence = max(0, readNumber(ui.edPeakProm, cfg.peakMinProminence));
cfg.maxPairSeparationUm = max(cfg.minPairSeparationUm, readNumber(ui.edMaxPair, cfg.maxPairSeparationUm));
cfg.minValleyDepthFraction = max(0, readNumber(ui.edValley, cfg.minValleyDepthFraction));
cfg.useExpectedGapGate = logical(get(ui.chkExpected, 'Value'));
cfg.livePlotDuringScan = logical(get(ui.chkLive, 'Value'));

% Keep limits valid.
cfg.contrastLow = max(0, min(0.99, cfg.contrastLow));
cfg.contrastHigh = max(cfg.contrastLow + eps, min(1, cfg.contrastHigh));
cfg.cannyLow = max(0, min(0.99, cfg.cannyLow));
cfg.cannyHigh = max(cfg.cannyLow + eps, min(1, cfg.cannyHigh));

set(ui.edCLow, 'String', num2str(cfg.contrastLow));
set(ui.edCHigh, 'String', num2str(cfg.contrastHigh));
set(ui.edCanLow, 'String', num2str(cfg.cannyLow));
set(ui.edCanHigh, 'String', num2str(cfg.cannyHigh));
end

function val = readNumber(h, fallback)
val = str2double(get(h, 'String'));
if isnan(val) || ~isfinite(val)
    val = fallback;
    set(h, 'String', num2str(fallback));
end
end

function out = normalize01(in)
in = double(in);
mn = min(in(:)); mx = max(in(:));
if mx > mn
    out = (in - mn) ./ (mx - mn);
else
    out = zeros(size(in));
end
end

function [cfg, state, ui] = getAll(fig)
cfg = getappdata(fig, 'cfg');
state = getappdata(fig, 'state');
ui = getappdata(fig, 'ui');
end

function setAll(fig, cfg, state, ui)
setappdata(fig, 'cfg', cfg);
setappdata(fig, 'state', state);
setappdata(fig, 'ui', ui);
end

function setStatus(fig, msg)
ui = getappdata(fig, 'ui');
if isfield(ui, 'txtStatus') && ishandle(ui.txtStatus)
    set(ui.txtStatus, 'String', sprintf('Status:\n%s', msg));
end
drawnow;
end

function a = wrapTo180Local(a)
a = mod(a + 180, 360) - 180;
end

function vname = findVar(T, candidates)
vname = "";
for c = string(candidates)
    idx = strcmpi(T.Properties.VariableNames, c);
    if any(idx)
        vname = T.Properties.VariableNames{find(idx,1)};
        return;
    end
end
error('None of the variables %s found in table.', strjoin(string(candidates), ', '));
end

function x = getNum(row, candidates)
for c = string(candidates)
    idx = strcmpi(row.Properties.VariableNames, c);
    if any(idx)
        varName = row.Properties.VariableNames{find(idx,1)};
        x = row.(varName)(1);
        if iscell(x), x = x{1}; end
        if ischar(x) || isstring(x), x = str2double(x); end
        if contains(lower(varName), 'mm')
            x = x * 1e-3;
        end
        return;
    end
end
error('Could not find variables: %s', strjoin(string(candidates), ', '));
end
