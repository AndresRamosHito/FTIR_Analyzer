function RawFtirAnalyzer()
% RAWFTIRANALYZER  FTIR batch analyzer with readable labels & color control (R2025b+)
%
% - Smooths (Savitzky–Golay), finds peaks, labels using getMoleculeLabel.
% - Plot styling: light theme, staggered labels, legend toggle.
% - Color control: name→RGB map with fuzzy matching.
% - Parallelization: computation only (no graphics in parfor).
%
% Author: (you)

% ---------- Tunables ----------
windowSize          = 11;     % SG frame (odd, > polyOrder)
polyOrder           = 2;      % SG poly (< windowSize)
minPeakProminence   = 0.01;   % relative to y-range
minAbsorbance       = 0.05;   % absolute absorbance floor
xlim_cm1            = [600 3500];
ylim_abs            = [];     % [] => auto
useParallel         = false;  % computation only
outDir              = fullfile(pwd,'FTIR_Figures');
saveSVG             = true;
dpi                 = 300;

% %% << you can edit >> Label text style: 'name' | 'wavenumber' | 'both'
labelStyle          = 'both';

% %% << you can edit >> Force colors per sample (keys are substrings, case-insensitive)
colorMap = containers.Map( ...
    lower([ ...
        "cattleya"; ...
        "coelogyne stigma (cristata)"; ...
        "coelogyne stigma (flaccida)"; ...
        "tormentosa" ...
    ]), ...
    { ...
        [0.00 0.00 0.00];      % black
        [0.12 0.47 0.71];      % blue
        [0.84 0.15 0.16];      % red
        [0.17 0.63 0.17]       % green
    } ...
);

% ---------- File pick ----------
[fileNames, pathName] = uigetfile({'*.csv;*.CSV','CSV files (*.csv)'}, ...
    'Select FTIR Data Files','MultiSelect','on');
if isequal(fileNames,0), disp('Canceled.'); return; end
if ~iscell(fileNames), fileNames = {fileNames}; end
filePaths = fullfile(pathName, fileNames);

% ---------- Output dir ----------
if ~exist(outDir,'dir'), mkdir(outDir); end

% ---------- Validate SG params ----------
validateattributes(windowSize,{'numeric'},{'scalar','odd','>=',5});
validateattributes(polyOrder, {'numeric'},{'scalar','>=',0,'<',windowSize});

% ---------- Combined figure (light theme) ----------
combinedFig = figure('Name','Combined FTIR Spectra','NumberTitle','off','Color','w');
tl = tiledlayout(combinedFig, 1, 1, 'TileSpacing','compact','Padding','compact');
axC = nexttile(tl); hold(axC,'on');

axC.Color      = 'w';
axC.GridColor  = [0.85 0.85 0.85];
axC.GridAlpha  = 1.0;
axC.MinorGridColor = [0.92 0.92 0.92];
axC.MinorGridAlpha = 1.0;
axC.XMinorGrid = 'on'; axC.YMinorGrid = 'on';
axC.Box        = 'on'; axC.LineWidth  = 1.0;

set(axC,'XDir','reverse');
xlabel(axC,'Wavenumber (cm^{-1})');
ylabel(axC,'Absorbance (a.u.)');
title(axC,'Combined FTIR Spectra');

set(groot,'DefaultAxesFontName','Arial','DefaultAxesFontSize',11,...
    'DefaultLineLineWidth',1.4);
colororder(axC, [...
    0.121 0.466 0.705;  1.000 0.498 0.054;  0.172 0.627 0.172; ...
    0.839 0.152 0.156;  0.580 0.404 0.741;  0.549 0.337 0.294; ...
    0.890 0.467 0.761;  0.498 0.498 0.498;  0.737 0.741 0.133; ...
    0.090 0.745 0.811]);

% ---------- Storage for interactivity ----------
plotData = struct('line',{},'text',{},'peaks',{});

% ---------- Optional parallel (no graphics in parfor) ----------
if useParallel && license('test','Distrib_Computing_Toolbox')
    p = gcp('nocreate'); if isempty(p), parpool; end
else
    useParallel = false;
end

% ---------- Compute first (no graphics) ----------
n = numel(filePaths);
names       = cell(n,1);
xCell       = cell(n,1);   % wavenumber (desc)
ySCell      = cell(n,1);   % smoothed absorbance
peakTables  = cell(n,1);   % table of peaks/labels

if useParallel
    parfor i = 1:n
        [names{i}, xCell{i}, ySCell{i}, peakTables{i}] = processOneFile( ...
            filePaths{i}, windowSize, polyOrder, minPeakProminence, ...
            minAbsorbance, xlim_cm1, outDir);
    end
else
    for i = 1:n
        [names{i}, xCell{i}, ySCell{i}, peakTables{i}] = processOneFile( ...
            filePaths{i}, windowSize, polyOrder, minPeakProminence, ...
            minAbsorbance, xlim_cm1, outDir);
    end
end

% ---------- Plot & save (serial for stability) ----------
for i = 1:n
    if isempty(xCell{i}), continue; end

    % color control
    thisColor = lookupColor(names{i}, colorMap);  % [] -> use colororder
    if isempty(thisColor)
        hLine = plot(axC, xCell{i}, ySCell{i}, 'DisplayName', names{i});
    else
        hLine = plot(axC, xCell{i}, ySCell{i}, 'Color', thisColor, 'DisplayName', names{i});
    end

    % build label text according to preference
    T = peakTables{i};
    if isempty(T)
        labText = strings(0,1);
    else
        switch lower(labelStyle)
            case 'name'
                labText = T.Label;
            case 'wavenumber'
                labText = string(compose('%.1f cm^{-1}', T.Wavenumber_cm1));
            otherwise
                % both
                labText = string(compose('%s (%.1f)', T.Label, T.Wavenumber_cm1));
        end
    end

    % Labels/stems on combined
    [thC, lhC] = drawPeakLabels(axC, T.Wavenumber_cm1, T.Absorbance, ...
                    'FontSize', 9, 'StemStyle', ':', 'LabelText', labText);

    plotData(end+1).line  = hLine; %#ok<SAGROW>
    plotData(end).text    = thC;
    plotData(end).peaks   = lhC;

    % Individual figure
    f = figure('Visible','off','Color','w','Name',['FTIR - ' names{i}]);
    axI = axes('Parent',f); hold(axI,'on');
    set(axI,'XDir','reverse'); box(axI,'on');
    grid(axI,'on'); axI.GridColor = [0.85 0.85 0.85]; axI.MinorGridColor = [0.92 0.92 0.92];

    if isempty(thisColor)
        plot(axI, xCell{i}, ySCell{i}, 'DisplayName', names{i});
    else
        plot(axI, xCell{i}, ySCell{i}, 'Color', thisColor, 'DisplayName', names{i});
    end
    xlabel(axI,'Wavenumber (cm^{-1})'); ylabel(axI,'Absorbance (a.u.)');
    title(axI, ['FTIR Spectrum - ' strrep(names{i},'_','\_')]);
    if ~isempty(xlim_cm1), xlim(axI, xlim_cm1); end
    if ~isempty(ylim_abs), ylim(axI, ylim_abs); end

    drawPeakLabels(axI, T.Wavenumber_cm1, T.Absorbance, ...
                   'FontSize',10, 'StemStyle','--', 'LabelText', labText);

    % Save individual
    indPNG = fullfile(outDir, sprintf('%s.png', names{i}));
    exportgraphics(f, indPNG, 'Resolution', dpi);
    if saveSVG
        indSVG = fullfile(outDir, sprintf('%s.svg', names{i}));
        exportgraphics(f, indSVG);
    end
    close(f);
end

% ---------- Finish combined ----------
if ~isempty(xlim_cm1), xlim(axC,xlim_cm1); end
if ~isempty(ylim_abs), ylim(axC,ylim_abs); end
grid(axC,'on');
lgd = legend(axC,'show','Location','northeastoutside'); title(lgd,'Samples');
lgd.TextColor = [0.1 0.1 0.1]; lgd.Color = [1 1 1];
lgd.Box = 'on'; lgd.EdgeColor = [0.85 0.85 0.85];
lgd.ItemHitFcn = @(~,evt) onLegendClick(evt, plotData);

% ---------- Save combined ----------
combPNG = fullfile(outDir,'Combined_FTIR_Spectra.png');
exportgraphics(combinedFig, combPNG, 'Resolution', dpi);
if saveSVG
    combSVG = fullfile(outDir,'Combined_FTIR_Spectra.svg');
    exportgraphics(combinedFig, combSVG);
end

% ---------- Optional HTML report ----------
try
    reportPath = fullfile(outDir,'FTIR_Analysis_Report.html');
    writeSimpleHTMLReport(reportPath, outDir, fileNames, peakTables, names);
    fprintf('Report: %s\n', reportPath);
catch ME
    warning('Report generation skipped: %s', ME.message);
end

disp('Done.');
end % ===== main =====


% ================== Helper: fuzzy color lookup ==================
function c = lookupColor(sampleName, colorMap)
% Finds the first key (substring) present in sampleName, case-insensitive.
c = [];
if isempty(colorMap), return; end
nm = lower(string(sampleName));
keys = string(colorMap.keys);
for k = 1:numel(keys)
    if contains(nm, keys(k))
        c = colorMap(lower(keys(k)));
        return;
    end
end
end


% ================== Helper: legend toggle ==================
function onLegendClick(evt, plotData)
clicked = evt.Peer;
isOn = strcmp(clicked.Visible,'on');
newVis = tern(isOn,'off','on');
clicked.Visible = newVis;
for k = 1:numel(plotData)
    if isequal(plotData(k).line, clicked)
        if ~isempty(plotData(k).text), set(plotData(k).text, 'Visible', newVis); end
        if ~isempty(plotData(k).peaks), set(plotData(k).peaks,'Visible', newVis); end
        break;
    end
end
end

function out = tern(cond,a,b)
if cond, out=a; else, out=b; end
end


% ================== Helper: simple HTML report ==================
function writeSimpleHTMLReport(htmlPath, outDir, fileNames, peakTables, names)
fid = fopen(htmlPath,'w'); assert(fid>0,'Cannot write report.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'<html><head><meta charset="utf-8"><title>FTIR Analysis Report</title>');
fprintf(fid,'<style>body{font-family:Arial;max-width:1000px;margin:24px} img{max-width:100%%;height:auto;border:1px solid #ddd} table{border-collapse:collapse} td,th{border:1px solid #ccc;padding:4px 8px}</style>');
fprintf(fid,'</head><body><h1>FTIR Analysis Report</h1>\n');

fprintf(fid,'<h2>Combined Spectra</h2>\n');
if exist(fullfile(outDir,'Combined_FTIR_Spectra.png'),'file')
    fprintf(fid,'<img src="%s" alt="Combined">\n','Combined_FTIR_Spectra.png');
end

fprintf(fid,'<h2>Samples</h2>\n');
for i = 1:numel(fileNames)
    nm = names{i}; if isempty(nm), continue; end
    png = sprintf('%s.png', nm);
    fprintf(fid,'<h3>%s</h3>\n', nm);
    if exist(fullfile(outDir,png),'file')
        fprintf(fid,'<img src="%s" alt="%s">\n', png, nm);
    end
    T = peakTables{i};
    if ~isempty(T)
        fprintf(fid,'<details><summary>Detected Peaks</summary>\n<table><tr><th>#</th><th>Wavenumber (cm^-1)</th><th>Label</th><th>Absorbance</th></tr>\n');
        for r=1:height(T)
            fprintf(fid,'<tr><td>%d</td><td>%.2f</td><td>%s</td><td>%.4f</td></tr>\n', ...
                T.Peak(r), T.Wavenumber_cm1(r), T.Label{r}, T.Absorbance(r));
        end
        fprintf(fid,'</table></details>\n');
    end
end
fprintf(fid,'</body></html>');
end


% ================== Peak label drawing with overlap avoidance ==================
function [txtHandles, stemHandles] = drawPeakLabels(ax, px, py, varargin)
% Draw peak labels with simple overlap avoidance and custom text.
% [txt, stems] = drawPeakLabels(ax, xPeaks, yPeaks, 'LabelText', stringArray, ...)
p = inputParser;
addParameter(p,'FontSize',9,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'StemStyle',':',@(s)ischar(s)||isstring(s));
addParameter(p,'LabelText',[],@(s) isempty(s) || isstring(s) || iscellstr(s));
parse(p,varargin{:});
fs   = p.Results.FontSize;
lsty = char(p.Results.StemStyle);
labText = p.Results.LabelText;

txtHandles  = gobjects(0);
stemHandles = gobjects(0);
if isempty(px) || isempty(py), return; end

px = px(:); py = py(:);
[px, idx] = sort(px,'descend'); py = py(idx);
if ~isempty(labText)
    labText = string(labText(:));
    labText = labText(idx);
else
    labText = string(compose('P%d (%.1f)', (1:numel(px)).', px));
end

yl = get(ax,'YLim'); dy = 0.03*range(yl);  % vertical spacing
dx_min = 15;                               % cm^-1 proximity to stagger

ypos = py; flip = 1;
for k = 2:numel(px)
    if abs(px(k) - px(k-1)) < dx_min && abs(ypos(k) - ypos(k-1)) < 2*dy
        ypos(k) = ypos(k) + flip*dy; flip = -flip;
    end
end

ybase = yl(1);
for k = 1:numel(px)
    stemHandles(end+1,1) = line(ax,[px(k) px(k)],[ybase py(k)], ...
        'Color',[0.25 0.25 0.25],'LineStyle',lsty,'HandleVisibility','off'); %#ok<AGROW>
    txtHandles(end+1,1) = text(ax, px(k), ypos(k)+0.01*range(yl), ...
        labText(k), 'Rotation',90, ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'FontSize',fs, 'Color',[0.1 0.1 0.1], 'Clipping','on'); %#ok<AGROW>
end
end


% ================== Core processing (no graphics; parfor-safe) ==================
function [name, xDesc, ySmooth, T] = processOneFile(fp, windowSize, polyOrder, ...
    minPeakProminence, minAbsorbance, xlim_cm1, outDir)
[~, name] = fileparts(fp);
xDesc = []; ySmooth = []; T = table();

raw = readmatrix(fp,'OutputType','double');
if isempty(raw) || size(raw,2) < 2
    raw = readmatrix(fp,'OutputType','double','NumHeaderLines',1);
end
if isempty(raw) || size(raw,2) < 2
    warning('Skipping "%s": could not read two numeric columns.', name);
    return;
end

x = raw(:,1); y = raw(:,2);
mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);

if ~isempty(xlim_cm1)
    mask2 = x <= xlim_cm1(2);
    x = x(mask2); y = y(mask2);
end

[x, order] = sort(x,'descend'); y = y(order);

if numel(y) >= windowSize
    yS = sgolayfilt(y, polyOrder, windowSize);
else
    yS = y;
end

prom = max(minPeakProminence * range(yS), eps);
[pks, locs] = findpeaks(yS, 'MinPeakProminence', prom);
if ~isempty(pks)
    valid = pks > minAbsorbance;
    pks = pks(valid); locs = locs(valid);
end
px = x(locs);

labels = strings(numel(px),1);
for j = 1:numel(px)
    labels(j) = getMoleculeLabel(px(j)); % <-- your mapping
end

T = table((1:numel(px)).', px(:), labels, pks(:), ...
    'VariableNames', {'Peak','Wavenumber_cm1','Label','Absorbance'});

if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir, sprintf('%s_peaks.csv', name)));

xDesc   = x;
ySmooth = yS;
end

% ================== Label map (your original mapping, preserved) ==================
function label = getMoleculeLabel(wavenumber)
% Assigns labels to wavenumbers based on the updated IR Spectrum Table
tolerance = 3;

% O–H Stretching and Bending Vibrations
if wavenumber >= 3710 - tolerance && wavenumber <= 3710 + tolerance
    label = 'Water in dilute solution';
elseif wavenumber <= 3600 && wavenumber >= 3100
    label = 'Water, H20';
elseif wavenumber <= 3650 && wavenumber >= 3590
    label = 'Free O–H (very sharp)';
elseif wavenumber <= 3600 && wavenumber >= 3200
    label = 'H-bonded O–H (strong)';
% =C–H (alkenyl) just above 3000
elseif wavenumber >= 3025 - tolerance && wavenumber <= 3005 + tolerance
    label = '=C–H stretch (alkenyl, unsaturation)';

% Aliphatic C–H stretching region (3000–2850): split CH3 vs CH2 (lipid acyl chains)
elseif wavenumber >= 2962 - tolerance && wavenumber <= 2950 + tolerance
    label = 'CH3 asymmetric stretch (acyl chains, lipids)';
elseif wavenumber >= 2935 - tolerance && wavenumber <= 2915 + tolerance
    label = 'CH2 asymmetric stretch (methylene, lipids)';
elseif wavenumber >= 2886 - tolerance && wavenumber <= 2866 + tolerance
    label = 'CH3 symmetric stretch (acyl chains, lipids)';
elseif wavenumber >= 2860 - tolerance && wavenumber <= 2848 + tolerance
    label = 'CH2 symmetric stretch (methylene, lipids)';

% Deformation modes that help confirm lipid chains
elseif wavenumber >= 1470 - tolerance && wavenumber <= 1460 + tolerance
    label = 'CH2 scissoring (methylene)';
elseif wavenumber >= 1456 - tolerance && wavenumber <= 1440 + tolerance
    label = 'CH3 asymmetric deformation';
elseif wavenumber >= 1385 - tolerance && wavenumber <= 1368 + tolerance
    label = 'CH3 symmetric deformation';

% Lipid ester carbonyl (triacylglycerols / phospholipids)
elseif wavenumber >= 1750 - tolerance && wavenumber <= 1730 + tolerance
    label = 'Ester C=O (lipids: TAGs/phospholipids)';

% Phosphate (lipid headgroups, nucleic acids)
elseif wavenumber >= 1246 - tolerance && wavenumber <= 1220 + tolerance
    label = 'PO2^- asymmetric stretch (phospholipids / nucleic acids)';
elseif wavenumber >= 1090 - tolerance && wavenumber <= 1078 + tolerance
    label = 'PO2^- symmetric stretch (phospholipids / nucleic acids)';

% Keep your existing O–H bending / C–O fallback rules below
elseif wavenumber <= 1410 && wavenumber >= 1260
    label = 'O–H bending (strong)';
elseif wavenumber <= 1150 && wavenumber >= 1040
    label = 'C–O stretching (strong)'

% N–H Stretching and Bending Vibrations
elseif wavenumber <= 3500 && wavenumber >= 3300
    label = 'N–H primary amines (medium)';
elseif wavenumber <= 3460 && wavenumber >= 3400
    label = '–CONH– (medium)';
elseif wavenumber <= 3100 && wavenumber >= 3070
    label = 'Additional band in solid state N–H (weak)';
elseif wavenumber <= 3130 && wavenumber >= 3030
    label = 'Amino acids –NH3+ (medium)';

% More from new data
elseif wavenumber >= 2349 - tolerance && wavenumber <= 2349 + tolerance
    label = 'Carbon dioxide O=C=O stretching (strong)';
elseif wavenumber <= 2305 && wavenumber >= 2280
    label = 'Nitrile oxides –C≡N+–O– (medium)';
elseif wavenumber <= 2275 && wavenumber >= 2250
    label = 'Isocyanates –N=C=O (strong)';
elseif wavenumber <= 2300 && wavenumber >= 2150
    label = 'Internal acetylenes –C≡C– (variable)';
elseif wavenumber <= 2175 && wavenumber >= 2140
    label = 'Thiocyanates –S–C≡N (strong)';
elseif wavenumber <= 2160 && wavenumber >= 2120
    label = 'Azides –N=N+=N– (strong)';
elseif wavenumber <= 2155 && wavenumber >= 2130
    label = 'Carbodiimides and Ketenes (strong)';
elseif wavenumber >= 2040 - tolerance && wavenumber <= 2040 + tolerance
    label = 'R3Si –C≡CH';
elseif wavenumber <= 2260 && wavenumber >= 2200
    label = 'Nitriles –C≡N (variable)';
elseif wavenumber <= 2180 && wavenumber >= 2120
    label = 'Isonitriles –N+≡C– (strong)';
elseif wavenumber <= 1870 && wavenumber >= 1820 || (wavenumber <= 1800 && wavenumber >= 1750)
    label = 'Acid Anhydrides: Saturated five-ring';
elseif wavenumber >= 1170 - tolerance && wavenumber <= 1170 + tolerance
    label = '-O-CS-N, Hypothiocyanite';
elseif wavenumber >= 1110 - tolerance && wavenumber <= 1110 + tolerance
    label = 'Uronic acid';
elseif wavenumber >= 1018 - tolerance && wavenumber <= 1018 + tolerance
    label = 'Uronic acid';
elseif wavenumber >= 961 && wavenumber <= 980
    label = 'r(CH2), beta-sheets';

% Acid Chlorides
elseif wavenumber <= 1815 && wavenumber >= 1790
    label = 'Acid Chlorides: Saturated';
elseif wavenumber <= 1790 && wavenumber >= 1750
    label = 'Acid Chlorides: Aryl and α,β-unsaturated';

% Acid Peroxides
elseif (wavenumber <= 1820 && wavenumber >= 1810) || (wavenumber <= 1800 && wavenumber >= 1780)
    label = 'Acid Peroxides: Saturated';
elseif (wavenumber <= 1805 && wavenumber >= 1780) || (wavenumber <= 1785 && wavenumber >= 1755)
    label = 'Acid Peroxides: Aryl and α,β-unsaturated';

% Esters and Lactones
elseif wavenumber <= 1750 && wavenumber >= 1735
    label = 'Esters and Lactones: Saturated';
elseif wavenumber <= 1730 && wavenumber >= 1715
    label = 'Esters and Lactones: Aryl and α,β-unsaturated';
elseif wavenumber <= 1800 && wavenumber >= 1750
    label = 'Aryl and vinyl esters (C=C–O–CO–)';
elseif wavenumber <= 1770 && wavenumber >= 1745
    label = 'Esters with electronegative α-substituents';
elseif wavenumber <= 1755 && wavenumber >= 1740
    label = 'α-Keto esters';
elseif wavenumber <= 1780 && wavenumber >= 1760
    label = 'Five-ring lactones';
elseif wavenumber <= 1770 && wavenumber >= 1740
    label = 'α,β-Unsaturated five-ring lactones';
elseif wavenumber >= 1800 - tolerance && wavenumber <= 1800 + tolerance
    label = 'β,γ-Unsaturated five-ring lactones';
elseif wavenumber >= 1650 - tolerance && wavenumber <= 1740 + tolerance
    label = 'β-Keto ester in H-bonding enol form';

% Aldehydes
elseif wavenumber <= 1740 && wavenumber >= 1720
    label = 'Aldehydes: Saturated';
elseif wavenumber <= 1715 && wavenumber >= 1695
    label = 'Aldehydes: Aryl';
elseif wavenumber <= 1765 && wavenumber >= 1695
    label = 'Aldehydes: α-Chloro or bromo';
elseif wavenumber <= 1705 && wavenumber >= 1680
    label = 'α,β-Unsaturated aldehyde';
elseif wavenumber <= 1680 && wavenumber >= 1660
    label = 'α, β, γ δ-Unsaturated aldehyde';
elseif wavenumber <= 1670 && wavenumber >= 1645
    label = 'β-Keto aldehyde in enol form (H-bonding)';

% Ketones
elseif wavenumber <= 1725 && wavenumber >= 1705
    label = 'Saturated ketone';
elseif wavenumber <= 1685 && wavenumber >= 1665
    label = 'α,β-Unsaturated ketone (often two bands)';
elseif wavenumber <= 1670 && wavenumber >= 1660
    label = 'α,β,α,β-Unsaturated and diaryl ketone';
elseif wavenumber <= 1705 && wavenumber >= 1685
    label = 'Cyclopropyl ketone';
elseif wavenumber <= 1750 && wavenumber >= 1740
    label = 'Five-ring ketone';
elseif wavenumber <= 1745 && wavenumber >= 1725
    label = 'α-Chloro or bromo ketone';
elseif wavenumber <= 1765 && wavenumber >= 1745
    label = 'α,α-Dichloro or dibromo ketone';
elseif (wavenumber >= 1760 - tolerance && wavenumber <= 1760 + tolerance || wavenumber >= 1730 - tolerance && wavenumber <= 1730 + tolerance)
    label = '1,2-Diketones s-cis six-ring';
elseif (wavenumber >= 1775 - tolerance && wavenumber <= 1775 + tolerance || wavenumber >= 1760 - tolerance && wavenumber <= 1760 + tolerance)
    label = '1,2-Diketones s-cis five-ring';
elseif (wavenumber >= 1650 - tolerance && wavenumber <= 1650 + tolerance || wavenumber >= 1615 - tolerance && wavenumber <= 1615 + tolerance)
    label = '1,3-Diketones enol form';
elseif wavenumber >= 1630 - tolerance && wavenumber <= 1640 + tolerance
    label = 'H20, O-H secondary peak';
elseif wavenumber <= 1655 && wavenumber >= 1635
    label = 'o-Hydroxy- or o-aminoaryl ketones';
elseif wavenumber <= 1645 && wavenumber >= 1615
    label = 'Diazoketones';
elseif wavenumber <= 1690 && wavenumber >= 1660
    label = 'Quinones';

% Carboxylic Acids
elseif wavenumber <= 3000 && wavenumber >= 2500
    label = 'Carboxylic acids R–CO2H (O–H stretching, H-bonded)';
elseif wavenumber <= 1725 && wavenumber >= 1700
    label = 'Saturated carboxylic acid';
elseif wavenumber <= 1715 && wavenumber >= 1690
    label = 'α,β-Unsaturated carboxylic acid';
elseif wavenumber <= 1700 && wavenumber >= 1680
    label = 'Aryl carboxylic acid';
elseif wavenumber >= 1780 - tolerance && wavenumber <= 1780 + tolerance
    label = 'Carbonate –O–CO–Cl';
elseif wavenumber >= 1740 - tolerance && wavenumber <= 1740 + tolerance
    label = 'Carbonate –O–CO–O–';
elseif wavenumber >= 1785 - tolerance && wavenumber <= 1785 + tolerance
    label = 'Aromatic Carbonate Ar–O–CO–O–Ar';
elseif wavenumber >= 1820 - tolerance && wavenumber <= 1820 + tolerance
    label = 'Five-ring Carbonate';
elseif wavenumber >= 1645 - tolerance && wavenumber <= 1645 + tolerance
    label = 'Thiocarbonate –S–CO–S–';
elseif wavenumber >= 1715 - tolerance && wavenumber <= 1715 + tolerance
    label = 'Aromatic Thiocarbonate Ar–S–CO–S–Ar';
elseif wavenumber >= 1640 - tolerance && wavenumber <= 1640 + tolerance
    label = 'Acylsilanes –CO–SiR3 (saturated)';
elseif wavenumber >= 1590 - tolerance && wavenumber <= 1590 + tolerance
    label = 'Acylsilanes –CO–SiR3 (α,β-unsaturated)';

% C=N; Imines, oximes, etc.
elseif wavenumber <= 3400 && wavenumber >= 3300
    label = 'N–H stretching (imines)';
elseif wavenumber <= 1690 && wavenumber >= 1640
    label = 'C=N stretching (imines, variable)';
elseif wavenumber <= 1600 && wavenumber >= 1630
    label = 'α,β-Unsaturated imines';

% Nitro compounds and related groups
elseif wavenumber <= 1570 && wavenumber >= 1540 || wavenumber <= 1390 && wavenumber >= 1340
    label = 'Nitro compounds C–NO2';
elseif wavenumber <= 1650 && wavenumber >= 1600 || wavenumber <= 1270 && wavenumber >= 1250
    label = 'Nitrates O–NO2';
elseif wavenumber <= 1630 && wavenumber >= 1550 || wavenumber <= 1300 && wavenumber >= 1250
    label = 'Nitramines N–NO2';
elseif wavenumber <= 1585 && wavenumber >= 1540
    label = 'Nitroso compounds (saturated)';
elseif wavenumber <= 1510 && wavenumber >= 1490
    label = 'Nitroso compounds (aryl)';
elseif wavenumber <= 1680 && wavenumber >= 1650
    label = 'Nitrites O–N=O (s-trans)';
elseif wavenumber <= 1625 && wavenumber >= 1610
    label = 'Nitrites O–N=O (s-cis)';
elseif wavenumber <= 1500 && wavenumber >= 1430
    label = 'N-Nitroso compounds N–N=O';
elseif wavenumber <= 970 && wavenumber >= 950
    label = 'N-Oxides (aliphatic)';
elseif wavenumber <= 1410 && wavenumber >= 1340 || wavenumber <= 860 && wavenumber >= 800
    label = 'Nitrate ions NO3–';

% Carbonates and Epoxides
elseif wavenumber >= 1041 && wavenumber <= 1060
    label = 'Cb-O C-OH Serine';
elseif wavenumber >= 1070 && wavenumber <= 1090
    label = 'C-O N-C Serine';
elseif wavenumber <= 1150 && wavenumber >= 1070
    label = 'C–O stretching';
elseif (wavenumber <= 1275 && wavenumber >= 1200) || (wavenumber <= 1075 && wavenumber >= 1020)
    label = 'C–O stretching in different environments';
elseif wavenumber <= 2850 && wavenumber >= 2810
    label = 'C–O–CH3 C–H stretching';
elseif wavenumber >= 1250 - tolerance && wavenumber <= 1250 + tolerance || wavenumber >= 900 - tolerance && wavenumber <= 900 + tolerance || wavenumber >= 800 - tolerance && wavenumber <= 800 + tolerance
    label = 'Epoxide bands';
elseif wavenumber <= 1550 && wavenumber >= 1480
    label = 'Aryl or B-N sretching';

% Boron Compounds
elseif wavenumber <= 2640 && wavenumber >= 2200
    label = 'CO2 or B–H stretching';
elseif wavenumber <= 1380 && wavenumber >= 1310
    label = 'B–O stretching (very strong)';
elseif wavenumber <= 1550 && wavenumber >= 1330
    label = 'B–N stretching (very strong)';

% Silicon Compounds
elseif wavenumber <= 2360 && wavenumber >= 2150 || (wavenumber <= 890 && wavenumber >= 860)
    label = 'Si–H stretching sensitive to electronegativity';
elseif wavenumber >= 2135 - tolerance && wavenumber <= 3690 + tolerance || (wavenumber <= 890 && wavenumber >= 860)
    label = 'Silicon hydride stretching';
elseif wavenumber <= 1275 && wavenumber >= 1245
    label = '–SiMen stretching';
elseif wavenumber >= 3690 - tolerance && wavenumber <= 3690 + tolerance
    label = 'Si–OH free O–H';
elseif wavenumber <= 3400 && wavenumber >= 3200
    label = 'Si–OH H-bonded O–H';
elseif wavenumber <= 1110 && wavenumber >= 1000
    label = 'Si–OR stretching';
elseif wavenumber <= 1080 && wavenumber >= 1040
    label = 'R3Si–O–SiR3 stretching';

% Phosphorus Compounds
elseif wavenumber <= 2440 && wavenumber >= 2350
    label = 'P–H stretching (sharp)';
elseif wavenumber >= 1440 - tolerance && wavenumber <= 1440 + tolerance
    label = 'P–Ph stretching (sharp)';
elseif wavenumber <= 1050 && wavenumber >= 1030
    label = 'P–O–alkyl stretching';
elseif wavenumber <= 1240 && wavenumber >= 1190
    label = 'P–O–aryl stretching';
elseif wavenumber <= 1300 && wavenumber >= 1250
    label = 'P=O stretching';
elseif wavenumber <= 750 && wavenumber >= 580
    label = 'P=S stretching';
elseif wavenumber <= 970 && wavenumber >= 910
    label = 'P–O–P stretching';
elseif wavenumber <= 2700 && wavenumber >= 2560
    label = 'P=O H-bonded O–H';

% Halogen and Inorganic Ions
elseif wavenumber <= 800 && wavenumber >= 600
    label = 'C–Cl stretching';
elseif wavenumber <= 750 && wavenumber >= 500
    label = 'C–Br stretching';
elseif wavenumber == 500
    label = 'C–I stretching';
elseif wavenumber <= 3300 && wavenumber >= 3030
    label = 'NH4+ bands';
elseif wavenumber <= 2200 && wavenumber >= 2000
    label = 'CN–, –SCN, –OCN bands';
elseif wavenumber <= 1450 && wavenumber >= 1410
    label = 'CO32– bands';
elseif wavenumber <= 1130 && wavenumber >= 1080
    label = 'SO42– bands';
elseif wavenumber <= 1380 && wavenumber >= 1350
    label = 'NO3– bands';
elseif wavenumber <= 1250 && wavenumber >= 1230
    label = 'NO2– bands';
elseif wavenumber <= 1100 && wavenumber >= 1000
    label = 'PO42– bands';
else
    label = 'to be added';
end
end

