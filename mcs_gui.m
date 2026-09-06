function app = mcs_gui(ncPath)
%MCS_GUI Interactive viewer for the Mores Creek Summit UAVSAR/lidar datacube.
%   app = MCS_GUI() opens ZARR/mcs_gui.nc (made by export_mcs_gui_nc.py).
%   Four linked panels: Terrain/Atmosphere, Local incidence angle, UAVSAR
%   products, and lidar products, each with dropdown product selection and
%   min/max color-range sliders. Returns a struct of handles + callbacks so
%   the test harness (mcs_gui_test.m) can drive it headlessly.

if nargin < 1, ncPath = fullfile('ZARR', 'mcs_gui.nc'); end      % default companion file location
assert(isfile(ncPath), 'NetCDF companion not found: %s (run export_mcs_gui_nc.py)', ncPath);

% ---------------- metadata: coordinates and per-dimension label tables ----------------
S.nc = ncPath;                                                   % everything hangs off one shared state struct
S.x = double(ncread(ncPath, 'x'));                               % UTM easting (m), ascending
yRaw = double(ncread(ncPath, 'y'));                              % UTM northing (m), stored north->south
S.flipY = yRaw(1) > yRaw(end);                                   % remember storage order so slices can be flipped
S.y = sort(yRaw, 'ascend');                                      % display axis wants ascending northing
S.pol = cellstr(ncreadStr(ncPath, 'pol'));                       % {'HH','HV','VH','VV'}
S.line = cellstr(ncreadStr(ncPath, 'line'));                     % {'05208','23205'}
S.pairLine = cellstr(ncreadStr(ncPath, 'line_id'));              % per-pair line id
S.pairD1 = cellstr(ncreadStr(ncPath, 'pair_date1'));             % per-pair start date
S.pairD2 = cellstr(ncreadStr(ncPath, 'pair_date2'));             % per-pair end date
S.pairDt = double(ncread(ncPath, 'temporal_baseline_days'));     % per-pair baseline (days)
S.heading = double(ncread(ncPath, 'heading'));                   % per-pair aircraft heading (deg)
S.flightLine = cellstr(ncreadStr(ncPath, 'flight_line_id'));     % per-flight line id
S.flightDate = cellstr(ncreadStr(ncPath, 'flight_date'));        % per-flight date
S.sdTime = cellstr(ncreadStr(ncPath, 'sd_time'));                % snow-depth dates (QSI + MCS_Lidar)
S.sdSource = cellstr(ncreadStr(ncPath, 'sd_source'));            % snow-depth product source per date
S.lidarTime = cellstr(ncreadStr(ncPath, 'lidar_time'));          % DTM/DSM/CHM dates (MCS_Lidar only)
S.cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); % slice cache: revisit layers instantly
S.ann = jsondecode(fileread('mcs_annotations.json'));            % Highway 21 + peak annotations (UTM)

% per-line heading labels for the global dropdown (heading is constant per line)
lineHeading = zeros(1, numel(S.line));                           % look up one pair's heading per line
for k = 1:numel(S.line), lineHeading(k) = S.heading(find(strcmp(S.pairLine, S.line{k}), 1)); end
lineLabels = arrayfun(@(k) sprintf('%s — %03d° (%s-looking)', S.line{k}, lineHeading(k), ...
    ternary(lineHeading(k) < 180, 'NE', 'SW')), 1:numel(S.line), 'UniformOutput', false);

% ---------------- figure and top-level layout ----------------
S.fig = uifigure('Name', 'Mores Creek Summit — UAVSAR/lidar pilot datacube', ...
                 'Position', [40 40 1500 950]);
main = uigridlayout(S.fig, [2 1], 'RowHeight', {40, '1x'}, 'Padding', [6 6 6 6]);

top = uigridlayout(main, [1 3], 'ColumnWidth', {110, 260, '1x'}, 'Padding', [0 0 0 0]); % control strip
uilabel(top, 'Text', 'Flight heading:', 'HorizontalAlignment', 'right');
S.ddHeading = uidropdown(top, 'Items', lineLabels, 'ItemsData', S.line, ...
    'Value', S.line{end}, 'ValueChangedFcn', @onHeading);        % default 23205 (most pairs)
uilabel(top, 'Text', sprintf('%s  |  grid: 3 m UTM 11N  |  pan/zoom linked across panels', ...
    strtrim(char(ncreadatt(ncPath, '/', 'site')))), 'FontColor', [.35 .35 .35]);

grid2 = uigridlayout(main, [2 2], 'Padding', [0 0 0 0]);         % the 2x2 panel grid

% each panel: [controls row; axes; slider row], built by a shared helper
[S.pTerr, S.axTerr] = makePanel(grid2, 'Terrain / Atmosphere');
[S.pLIA,  S.axLIA]  = makePanel(grid2, 'Local incidence angle');
[S.pSAR,  S.axSAR]  = makePanel(grid2, 'UAVSAR');
[S.pLid,  S.axLid]  = makePanel(grid2, 'Lidar');
linkaxes([S.axTerr S.axLIA S.axSAR S.axLid]);                    % shared pan/zoom

% ---------------- panel-specific controls ----------------
% Terrain: product (DEM | atm delay) + flight date (enabled for atm delay only)
S.ddTerrProd = addDropdown(S.pTerr, {'DEM (2023-02-09 lidar DTM)', 'Atmospheric delay'}, ...
    {'dem', 'atm'}, @(varargin) refreshTerr());
S.ddTerrFlight = addDropdown(S.pTerr, {'-'}, {''}, @(varargin) refreshTerr());
S.ddTerrFlight.Enable = 'off';                                   % greyed out until atm delay chosen

% UAVSAR: product + pair + polarization
S.ddSARProd = addDropdown(S.pSAR, {'Coherence', 'Unwrapped phase', 'Wrapped phase (angle)', ...
    'Wrapped phase (amplitude)', 'Atm delay difference', 'dSWE (m w.e.)'}, ...
    {'coherence', 'unwrapped_phase', 'wrapped_angle', 'wrapped_amplitude', 'atm_delay_diff', 'dswe'}, ...
    @(varargin) refreshSAR());
S.ddSARPair = addDropdown(S.pSAR, {'-'}, {1}, @(varargin) refreshSAR());
S.ddSARPol = addDropdown(S.pSAR, S.pol, num2cell(1:numel(S.pol)), @(varargin) refreshSAR());

% Lidar: product + date (list depends on product: sd has 15 dates, terrain products 13)
S.ddLidProd = addDropdown(S.pLid, {'Snow depth', 'DTM', 'DSM', 'CHM'}, ...
    {'sd', 'dtm', 'dsm', 'chm'}, @onLidProd);
S.ddLidDate = addDropdown(S.pLid, {'-'}, {1}, @(varargin) refreshLid());

% ---------------- initial population and draw ----------------
onHeading();                                                     % fills pair/flight lists + draws Terr/LIA/SAR
onLidProd();                                                     % fills lidar date list + draws lidar panel

% handles + callbacks returned for the headless test harness
app = S;
app.cb = struct('heading', @onHeading, 'terr', @(varargin) refreshTerr(), ...
                'sar', @(varargin) refreshSAR(), 'lidProd', @onLidProd, ...
                'lid', @(varargin) refreshLid());

% ================= nested callbacks (share S via closure) =================

    function onHeading(varargin)
        ln = S.ddHeading.Value;                                  % selected line id
        pairIdx = find(strcmp(S.pairLine, ln));                  % pairs flown on this line
        labels = arrayfun(@(i) sprintf('%s → %s (%dd)', S.pairD1{i}, S.pairD2{i}, S.pairDt(i)), ...
            pairIdx, 'UniformOutput', false);
        setItems(S.ddSARPair, labels, num2cell(pairIdx(:)'));    % repopulate, keep selection if possible
        fIdx = find(strcmp(S.flightLine, ln));                   % flights on this line (for atm delay)
        setItems(S.ddTerrFlight, S.flightDate(fIdx), num2cell(fIdx(:)'));
        refreshLIA(); refreshTerr(); refreshSAR();               % lidar panel is heading-independent
    end

    function refreshLIA()
        li = find(strcmp(S.line, S.ddHeading.Value));            % index into the 'line' dimension
        img = getSlice('lia', {li});                             % (y,x) slice for this line
        show(S.axLIA, img, sprintf('LIA — line %s (deg)', S.ddHeading.Value), cmapSeq());
    end

    function refreshTerr()
        if strcmp(S.ddTerrProd.Value, 'dem')
            S.ddTerrFlight.Enable = 'off';                       % flight choice is meaningless for the DEM
            show(S.axTerr, getSlice('dem', {}), 'DEM — 2023-02-09 lidar DTM (m)', cmapTerrain());
        else
            S.ddTerrFlight.Enable = 'on';
            fi = S.ddTerrFlight.Value;                           % index into the 'flight' dimension
            show(S.axTerr, getSlice('atm_delay', {fi}), ...  % one-signed total delay -> sequential map
                sprintf('Atm delay — %s %s (rad)', S.flightLine{fi}, S.flightDate{fi}), cmapSeq());
        end
    end

    function refreshSAR()
        prod = S.ddSARProd.Value; pi = S.ddSARPair.Value;        % product name + pair index
        isPairLevel = strcmp(prod, 'atm_delay_diff');            % delay difference has no polarization
        S.ddSARPol.Enable = ternary(isPairLevel, 'off', 'on');
        if isPairLevel
            img = getSlice(prod, {pi});
            polTxt = '';
        else
            img = getSlice(prod, {pi, S.ddSARPol.Value});        % (pair, pol) indexed slice
            polTxt = [' ' S.pol{S.ddSARPol.Value}];
        end
        names = containers.Map({'coherence', 'unwrapped_phase', 'wrapped_angle', ...
            'wrapped_amplitude', 'atm_delay_diff', 'dswe'}, {'Coherence', 'Unwrapped phase (rad)', ...
            'Wrapped phase (rad)', 'Wrapped amplitude', 'Atm delay diff (rad)', ...
            'dSWE (m w.e., density 250, SNOTEL-anchored)'});
        maps = containers.Map({'coherence', 'unwrapped_phase', 'wrapped_angle', ...
            'wrapped_amplitude', 'atm_delay_diff', 'dswe'}, {cmapSeq(), cmapDiverging(), hsv(256), ...
            cmapSeq(), cmapDiverging(), cmapDiverging()});
        signed = ismember(prod, {'unwrapped_phase', 'atm_delay_diff', 'dswe'});  % zero-centered products
        show(S.axSAR, img, sprintf('%s%s — %s → %s', names(prod), polTxt, ...
            S.pairD1{pi}, S.pairD2{pi}), maps(prod), signed);
    end

    function onLidProd(varargin)
        if strcmp(S.ddLidProd.Value, 'sd')                       % snow depth: 15 dates, show source
            labels = cellfun(@(d, s) sprintf('%s (%s)', d, s), S.sdTime, S.sdSource, ...
                'UniformOutput', false);
            setItems(S.ddLidDate, labels, num2cell(1:numel(S.sdTime)));
        else                                                     % DTM/DSM/CHM: 13 MCS_Lidar dates
            setItems(S.ddLidDate, S.lidarTime, num2cell(1:numel(S.lidarTime)));
        end
        refreshLid();
    end

    function refreshLid()
        prod = S.ddLidProd.Value; ti = S.ddLidDate.Value;        % product + time index
        img = getSlice(prod, {ti});
        names = containers.Map({'sd', 'dtm', 'dsm', 'chm'}, ...
            {'Snow depth (m)', 'DTM (m)', 'DSM (m)', 'CHM (m)'});
        maps = containers.Map({'sd', 'dtm', 'dsm', 'chm'}, ...
            {cmapSeq(), cmapTerrain(), cmapTerrain(), cmapSeq()});
        % branch, not ternary(): both time axes get indexed otherwise, and sd has 15
        % dates while lidar_time has only 13
        if strcmp(prod, 'sd'), dateTxt = S.sdTime{ti}; else, dateTxt = S.lidarTime{ti}; end
        show(S.axLid, img, sprintf('%s — %s', names(prod), dateTxt), maps(prod));
    end

% ================= data access =================

    function img = getSlice(var, idx)
        % Read one (y,x) slice by trailing-dimension indices, cached by key.
        key = [var sprintf('_%d', idx{:})];
        if S.cache.isKey(key), img = S.cache(key); return; end
        % MATLAB ncread reverses dim order vs the file: file (…,y,x) -> MATLAB (x,y,…)
        nExtra = numel(idx);                                     % how many non-spatial dims this var has
        start = [1 1 fliplr(cell2mat(idx))];                     % x,y first, then reversed extra dims
        count = [Inf Inf ones(1, nExtra)];
        A = ncread(S.nc, var, start(1:2 + nExtra), count(1:2 + nExtra));
        img = A';                                                % (x,y) -> (y,x) for display
        if S.flipY, img = flipud(img); end                       % storage is north->south; axis is ascending
        S.cache(key) = img;
    end

% ================= drawing =================

    function [panel, ax] = makePanel(parent, title)
        % One display panel: dropdown row, image axes, and min/max clim sliders.
        panel = uigridlayout(parent, [3 1], 'RowHeight', {30, '1x', 34}, 'Padding', [2 2 2 2]);
        panel.UserData = struct('title', title, 'ctrlRow', [], 'sliderRow', []);
        ctrl = uigridlayout(panel, [1 4], 'Padding', [0 0 0 0], 'ColumnSpacing', 4);
        ctrl.ColumnWidth = repmat({'fit'}, 1, 4);
        panel.UserData.ctrlRow = ctrl;
        ax = uiaxes(panel);                                      % the image axes
        ax.Layout.Row = 2;
        axis(ax, 'equal'); ax.Color = [0.88 0.88 0.88];          % gray background shows through NaN
        ax.XLim = [S.x(1) S.x(end)]; ax.YLim = [S.y(1) S.y(end)];
        xlabel(ax, 'UTM easting (m)'); ylabel(ax, 'UTM northing (m)');
        colorbar(ax);
        srow = uigridlayout(panel, [1 5], 'Padding', [0 0 0 0], 'ColumnSpacing', 4, ...
            'ColumnWidth', {30, '1x', 30, '1x', 46});
        panel.UserData.sliderRow = srow;
        uilabel(srow, 'Text', 'min');
        sMin = uislider(srow, 'ValueChangingFcn', @(s, e) climFromSliders(ax, e, 1));
        uilabel(srow, 'Text', 'max');
        sMax = uislider(srow, 'ValueChangingFcn', @(s, e) climFromSliders(ax, e, 2));
        uibutton(srow, 'Text', 'reset', 'ButtonPushedFcn', @(varargin) resetClim(ax));
        ax.UserData = struct('im', [], 'sMin', sMin, 'sMax', sMax, 'noData', [], 'img', []);
    end

    function show(ax, img, titleTxt, cmap, signed)
        % Draw/update one panel: image, NaN transparency, robust clim, sliders, annotations.
        if nargin < 5, signed = false; end                       % signed -> clim symmetric about 0
        u = ax.UserData;
        u.img = img; u.signed = signed;
        if isempty(u.im)                                         % first draw: create image + overlays once
            u.im = imagesc(ax, [S.x(1) S.x(end)], [S.y(1) S.y(end)], img);
            ax.YDir = 'normal';                                  % ascending northing upward
            hold(ax, 'on');
            for s = 1:numel(S.ann.highway21.segments)            % Highway 21 polyline (may be cell or array)
                seg = S.ann.highway21.segments{s};
                plot(ax, seg(:, 1), seg(:, 2), '-', 'Color', [0.9 0.1 0.1], 'LineWidth', 1.4);
            end
            for p = 1:numel(S.ann.peaks)                         % peak markers + labels
                pk = S.ann.peaks(p);
                plot(ax, pk.x, pk.y, '^k', 'MarkerFaceColor', 'w', 'MarkerSize', 8);
                text(ax, pk.x + 120, pk.y, pk.name, 'FontWeight', 'bold', 'FontSize', 10, ...
                    'Color', 'k', 'BackgroundColor', [1 1 1 0.55], 'Margin', 1);
            end
            u.noData = text(ax, mean(S.x([1 end])), mean(S.y([1 end])), ...
                'no data for this selection', 'FontSize', 16, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'Visible', 'off');
        else
            u.im.CData = img;                                    % subsequent draws just swap the data
        end
        u.im.AlphaData = isfinite(img);                          % NaN -> transparent (gray background)
        title(ax, titleTxt, 'Interpreter', 'none');
        colormap(ax, cmap);
        allNaN = ~any(isfinite(img), 'all');
        u.noData.Visible = ternary(allNaN, 'on', 'off');         % explicit message instead of a blank axes
        ax.UserData = u;
        if ~allNaN, resetClim(ax); end                           % robust color limits + slider ranges
    end

    function resetClim(ax)
        % Robust 2/98 percentile color limits; sliders span the slice's full range.
        u = ax.UserData;
        v = u.img(isfinite(u.img));
        if isempty(v), return; end
        if numel(v) > 2e5, v = v(1:ceil(numel(v) / 2e5):end); end % subsample: prctile on ~200k values
        lims = double(prctile(v, [2 98]));
        full = double([min(v) max(v)]);
        if isfield(u, 'signed') && u.signed                      % diverging map: white pinned at zero
            lims = max(abs(lims)) * [-1 1];
            full = max(abs(full)) * [-1 1];
        end
        if full(1) >= full(2), full = full(1) + [-0.5 0.5]; end  % degenerate (constant) slice
        if lims(1) >= lims(2), lims = full; end
        u.sMin.Limits = full; u.sMax.Limits = full;              % sliders cover the full data range
        u.sMin.Value = max(lims(1), full(1)); u.sMax.Value = min(lims(2), full(2));
        clim(ax, [u.sMin.Value u.sMax.Value]);
    end

    function climFromSliders(ax, e, which)
        % Live clim update while dragging either slider; keep min strictly below max.
        u = ax.UserData;
        lo = u.sMin.Value; hi = u.sMax.Value;
        if which == 1, lo = e.Value; else, hi = e.Value; end     % use the in-drag value
        if lo >= hi, return; end                                 % ignore crossed sliders
        clim(ax, [lo hi]);
    end

% ================= small helpers =================

    function dd = addDropdown(panel, items, itemsData, cb)
        % Add one dropdown to a panel's control row.
        dd = uidropdown(panel.UserData.ctrlRow, 'Items', items, 'ItemsData', itemsData, ...
            'ValueChangedFcn', cb);
    end

    function setItems(dd, items, itemsData)
        % Repopulate a dropdown, preserving the current selection when still valid.
        old = dd.Value;
        dd.Items = items; dd.ItemsData = itemsData;
        if any(cellfun(@(v) isequal(v, old), itemsData)), dd.Value = old; end
    end
end

% ================= file-local helpers =================

function s = ncreadStr(nc, var)
% Read a string/char netCDF variable as a MATLAB string array.
raw = ncread(nc, var);
if ischar(raw), s = string(raw'); else, s = string(raw); end     % char matrix vs NC_STRING
s = strtrim(s(:));
end

function out = ternary(cond, a, b)
% Inline conditional (MATLAB has none).
if cond, out = a; else, out = b; end
end

function m = cmapTerrain()
% Green -> tan -> brown -> white elevation ramp.
anchors = [0.20 0.45 0.20; 0.55 0.65 0.30; 0.80 0.70 0.45; 0.65 0.50 0.35; 0.95 0.95 0.95];
m = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, 256));
end

function m = cmapDiverging()
% Blue -> white -> red, for signed quantities (phase, delay differences).
anchors = [0.10 0.25 0.70; 1 1 1; 0.75 0.10 0.15];
m = interp1([0 0.5 1], anchors, linspace(0, 1, 256));
end

function m = cmapSeq()
% Sequential default.
m = parula(256);
end
