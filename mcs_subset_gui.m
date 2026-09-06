function app = mcs_subset_gui(ncPath)
%MCS_SUBSET_GUI Viewer for the Mores Creek Summit modeling-team subset file.
%   app = MCS_SUBSET_GUI() opens mcs_top_pairs_subset.nc (looked for in ZARR/
%   and the current directory). Standalone: needs only this file, the .nc, and
%   (optionally) mcs_annotations.json for the Highway 21 / peak overlays.
%
%   Four linked panels: DEM, local incidence angle (follows the selected pair's
%   flight line), UAVSAR products (dSWE / coherence / unwrapped phase / atm
%   delay difference), and QSI lidar snow depth. Global pair dropdown uses the
%   letters A-F from the summary report; per-panel min/max color sliders.

% ---------------- locate inputs ----------------
if nargin < 1
    cands = {fullfile('ZARR', 'mcs_top_pairs_subset.nc'), 'mcs_top_pairs_subset.nc'};
    ncPath = cands{find(cellfun(@isfile, cands), 1)};
    assert(~isempty(ncPath), ['mcs_top_pairs_subset.nc not found in ZARR/ or the current ' ...
        'directory (regenerate with export_subset_nc.py)']);
end

% ---------------- metadata ----------------
S.nc = ncPath;
S.x = double(ncread(ncPath, 'x'));                               % UTM easting, ascending
yRaw = double(ncread(ncPath, 'y'));                              % stored north -> south
S.flipY = yRaw(1) > yRaw(end);
S.y = sort(yRaw, 'ascend');
S.pol = cellstr(strtrim(string(ncread(ncPath, 'pol'))));
S.line = cellstr(strtrim(string(ncread(ncPath, 'line'))));
S.pairLine = cellstr(strtrim(string(ncread(ncPath, 'line_id'))));
S.pairD1 = cellstr(strtrim(string(ncread(ncPath, 'pair_date1'))));
S.pairD2 = cellstr(strtrim(string(ncread(ncPath, 'pair_date2'))));
S.rankNote = cellstr(strtrim(string(ncread(ncPath, 'rank_note'))));
S.sdTime = cellstr(strtrim(string(ncread(ncPath, 'sd_time'))));
S.cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
S.ann = [];                                                      % annotations are optional
if isfile('mcs_annotations.json')
    S.ann = jsondecode(fileread('mcs_annotations.json'));
end

nPair = numel(S.pairD1);
pairLabels = arrayfun(@(i) sprintf('%s — %s %s → %s', S.rankNote{i}, S.pairLine{i}, ...
    S.pairD1{i}, S.pairD2{i}), 1:nPair, 'UniformOutput', false);

% ---------------- figure and layout ----------------
S.fig = uifigure('Name', 'MCS dSWE subset — top pairs for snow-model input', ...
                 'Position', [40 40 1500 950]);
main = uigridlayout(S.fig, [2 1], 'RowHeight', {40, '1x'}, 'Padding', [6 6 6 6]);
top = uigridlayout(main, [1 4], 'ColumnWidth', {45, 420, 110, '1x'}, 'Padding', [0 0 0 0]);
uilabel(top, 'Text', 'Pair:', 'HorizontalAlignment', 'right');
S.ddPair = uidropdown(top, 'Items', pairLabels, 'ItemsData', num2cell(1:nPair), ...
    'ValueChangedFcn', @onPair);
uilabel(top, 'Text', 'Polarization:', 'HorizontalAlignment', 'right');
S.ddPol = uidropdown(top, 'Items', S.pol, 'ItemsData', num2cell(1:numel(S.pol)), ...
    'ValueChangedFcn', @(varargin) refreshSAR());

grid2 = uigridlayout(main, [2 2], 'Padding', [0 0 0 0]);
[S.pDEM, S.axDEM] = makePanel(grid2);
[S.pLIA, S.axLIA] = makePanel(grid2);
[S.pSAR, S.axSAR] = makePanel(grid2);
[S.pLid, S.axLid] = makePanel(grid2);
linkaxes([S.axDEM S.axLIA S.axSAR S.axLid]);                     % shared pan/zoom

% panel controls: UAVSAR product dropdown, lidar date dropdown (DEM/LIA need none)
S.ddSARProd = uidropdown(S.pSAR.UserData.ctrlRow, ...
    'Items', {'dSWE (m w.e.)', 'Coherence', 'Unwrapped phase', 'Atm delay difference'}, ...
    'ItemsData', {'dswe', 'coherence', 'unwrapped_phase', 'atm_delay_diff'}, ...
    'ValueChangedFcn', @(varargin) refreshSAR());
sdLabels = cellfun(@(d) sprintf('QSI snow depth %s', d), S.sdTime, 'UniformOutput', false);
S.ddLidDate = uidropdown(S.pLid.UserData.ctrlRow, 'Items', sdLabels, ...
    'ItemsData', num2cell(1:numel(S.sdTime)), 'ValueChangedFcn', @(varargin) refreshLid());

% ---------------- initial draw ----------------
show(S.axDEM, getSlice('dem', {}), 'DEM — 2023-02-09 lidar DTM (m)', cmapTerrain());
onPair(); refreshLid();

app = S;
app.cb = struct('pair', @onPair, 'sar', @(varargin) refreshSAR(), ...
                'lid', @(varargin) refreshLid());

% ================= nested callbacks =================

    function onPair(varargin)
        pi = S.ddPair.Value;
        li = find(strcmp(S.line, S.pairLine{pi}));               % LIA follows the pair's line
        show(S.axLIA, getSlice('lia', {li}), ...
            sprintf('Local incidence angle — line %s (deg)', S.pairLine{pi}), cmapSeq());
        refreshSAR();
    end

    function refreshSAR()
        pi = S.ddPair.Value; prod = S.ddSARProd.Value;
        isPairLevel = strcmp(prod, 'atm_delay_diff');            % no polarization dimension
        S.ddPol.Enable = onoff(~isPairLevel);
        if isPairLevel
            img = getSlice(prod, {pi}); polTxt = '';
        else
            img = getSlice(prod, {pi, S.ddPol.Value});
            polTxt = [' ' S.pol{S.ddPol.Value}];
        end
        names = containers.Map({'dswe', 'coherence', 'unwrapped_phase', 'atm_delay_diff'}, ...
            {'dSWE (m w.e., density 250, SNOTEL-anchored)', 'Coherence', ...
             'Unwrapped phase (rad)', 'Atm delay diff (rad)'});
        maps = containers.Map({'dswe', 'coherence', 'unwrapped_phase', 'atm_delay_diff'}, ...
            {cmapDiverging(), cmapSeq(), cmapDiverging(), cmapDiverging()});
        letter = regexp(S.rankNote{pi}, 'Pair \w', 'match', 'once');
        % UAVSAR products auto-span the 5-95% quantiles of the displayed slice
        show(S.axSAR, img, sprintf('%s: %s%s — %s → %s', letter, names(prod), polTxt, ...
            S.pairD1{pi}, S.pairD2{pi}), maps(prod), ...
            ismember(prod, {'dswe', 'unwrapped_phase', 'atm_delay_diff'}), [5 95]);
    end

    function refreshLid()
        ti = S.ddLidDate.Value;
        show(S.axLid, getSlice('sd', {ti}), ...
            sprintf('QSI snow depth (m) — %s', S.sdTime{ti}), cmapSeq());
    end

% ================= data access =================

    function img = getSlice(var, idx)
        key = [var sprintf('_%d', idx{:})];
        if S.cache.isKey(key), img = S.cache(key); return; end
        % MATLAB ncread reverses dim order vs the file: (..., y, x) -> (x, y, ...)
        nExtra = numel(idx);
        start = [1 1 fliplr(cell2mat(idx))];
        count = [Inf Inf ones(1, nExtra)];
        A = ncread(S.nc, var, start(1:2 + nExtra), count(1:2 + nExtra));
        img = A';
        if S.flipY, img = flipud(img); end
        S.cache(key) = img;
    end

% ================= drawing (same conventions as mcs_gui.m) =================

    function [panel, ax] = makePanel(parent)
        panel = uigridlayout(parent, [3 1], 'RowHeight', {30, '1x', 34}, 'Padding', [2 2 2 2]);
        ctrl = uigridlayout(panel, [1 2], 'Padding', [0 0 0 0], 'ColumnSpacing', 4);
        ctrl.ColumnWidth = repmat({'fit'}, 1, 2);
        panel.UserData = struct('ctrlRow', ctrl);
        ax = uiaxes(panel); ax.Layout.Row = 2;
        axis(ax, 'equal'); ax.Color = [0.88 0.88 0.88];
        ax.XLim = [S.x(1) S.x(end)]; ax.YLim = [S.y(1) S.y(end)];
        xlabel(ax, 'UTM easting (m)'); ylabel(ax, 'UTM northing (m)');
        colorbar(ax);
        srow = uigridlayout(panel, [1 5], 'Padding', [0 0 0 0], 'ColumnSpacing', 4, ...
            'ColumnWidth', {30, '1x', 30, '1x', 46});
        uilabel(srow, 'Text', 'min');
        sMin = uislider(srow, 'ValueChangingFcn', @(s, e) climDrag(ax, e, 1));
        uilabel(srow, 'Text', 'max');
        sMax = uislider(srow, 'ValueChangingFcn', @(s, e) climDrag(ax, e, 2));
        uibutton(srow, 'Text', 'reset', 'ButtonPushedFcn', @(varargin) resetClim(ax));
        ax.UserData = struct('im', [], 'sMin', sMin, 'sMax', sMax, 'noData', [], ...
                             'img', [], 'signed', false);
    end

    function show(ax, img, titleTxt, cmap, signed, qlims)
        if nargin < 5, signed = false; end
        if nargin < 6, qlims = [2 98]; end                       % default robust range
        u = ax.UserData; u.img = img; u.signed = signed; u.qlims = qlims;
        if isempty(u.im)
            u.im = imagesc(ax, [S.x(1) S.x(end)], [S.y(1) S.y(end)], img);
            ax.YDir = 'normal'; hold(ax, 'on');
            if ~isempty(S.ann)                                   % optional overlays
                for si = 1:numel(S.ann.highway21.segments)
                    seg = S.ann.highway21.segments{si};
                    plot(ax, seg(:, 1), seg(:, 2), '-', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.4);
                end
                for p = 1:numel(S.ann.peaks)
                    pk = S.ann.peaks(p);
                    plot(ax, pk.x, pk.y, '^k', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
                    text(ax, pk.x + 120, pk.y, pk.name, 'FontWeight', 'bold', 'FontSize', 10, ...
                        'BackgroundColor', [1 1 1 0.55], 'Margin', 1);
                end
            end
            u.noData = text(ax, mean(S.x([1 end])), mean(S.y([1 end])), ...
                'no data for this selection', 'FontSize', 16, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Visible', 'off');
        else
            u.im.CData = img;
        end
        u.im.AlphaData = isfinite(img);
        title(ax, titleTxt, 'Interpreter', 'none');
        colormap(ax, cmap);
        allNaN = ~any(isfinite(img), 'all');
        u.noData.Visible = onoff(allNaN);
        ax.UserData = u;
        if ~allNaN, resetClim(ax); end
    end

    function resetClim(ax)
        u = ax.UserData;
        v = u.img(isfinite(u.img));
        if isempty(v), return; end
        if numel(v) > 2e5, v = v(1:ceil(numel(v) / 2e5):end); end
        q = [2 98];
        if isfield(u, 'qlims') && ~isempty(u.qlims), q = u.qlims; end
        lims = double(prctile(v, q));
        full = double([min(v) max(v)]);
        if full(1) >= full(2), full = full(1) + [-0.5 0.5]; end
        if lims(1) >= lims(2), lims = full; end
        if u.signed                                              % zero-centered for signed products
            lims = max(abs(lims)) * [-1 1];
            full = max(abs(full)) * [-1 1];
        end
        u.sMin.Limits = full; u.sMax.Limits = full;
        u.sMin.Value = max(lims(1), full(1)); u.sMax.Value = min(lims(2), full(2));
        clim(ax, [u.sMin.Value u.sMax.Value]);
    end

    function climDrag(ax, e, which)
        u = ax.UserData;
        lo = u.sMin.Value; hi = u.sMax.Value;
        if which == 1, lo = e.Value; else, hi = e.Value; end
        if lo < hi, clim(ax, [lo hi]); end
    end
end

% ================= file-local helpers =================

function s = onoff(tf)
% logical -> 'on'/'off' (single-evaluation, unlike an inline ternary)
if tf, s = 'on'; else, s = 'off'; end
end

function m = cmapTerrain()
anchors = [0.20 0.45 0.20; 0.55 0.65 0.30; 0.80 0.70 0.45; 0.65 0.50 0.35; 0.95 0.95 0.95];
m = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 256));
end

function m = cmapDiverging()
anchors = [0.10 0.25 0.70; 1 1 1; 0.75 0.10 0.15];
m = interp1([0 0.5 1], anchors, linspace(0, 1, 256));
end

function m = cmapSeq()
m = parula(256);
end
