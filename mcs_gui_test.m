function mcs_gui_test(outDir)
%MCS_GUI_TEST Headless harness for mcs_gui: drive every control, export PNGs.
%   Builds the app, walks both headings x every product (representative
%   pairs/pols/dates), asserts no callback errors, and writes one PNG per
%   state to outDir (default 'gui_test_png') for visual review.

if nargin < 1, outDir = 'gui_test_png'; end                      % where the screenshots land
if ~isfolder(outDir), mkdir(outDir); end

app = mcs_gui();                                                 % build the app on the default .nc
snap = @(name) snapRetry(app.fig, fullfile(outDir, [name '.png'])); % whole-figure screenshot (with retry)
n = 0;                                                           % screenshot counter

for ln = app.line'                                               % both flight headings
    app.ddHeading.Value = char(ln);
    app.cb.heading();                                            % repopulates pair/flight lists, redraws

    % ---- Terrain panel: DEM, then first + last atm-delay flight of this heading ----
    app.ddTerrProd.Value = 'dem'; app.cb.terr();
    n = n + 1; snap(sprintf('%02d_%s_terr_dem', n, char(ln)));
    app.ddTerrProd.Value = 'atm'; app.cb.terr();
    flights = cell2mat(app.ddTerrFlight.ItemsData);
    for fi = unique([flights(1) flights(end)])
        app.ddTerrFlight.Value = fi; app.cb.terr();
        n = n + 1; snap(sprintf('%02d_%s_terr_atm_f%02d', n, char(ln), fi));
    end

    % ---- UAVSAR panel: every product on first/middle/last pair, two pols ----
    pairs = cell2mat(app.ddSARPair.ItemsData);
    testPairs = unique([pairs(1) pairs(ceil(end / 2)) pairs(end)]);
    for prod = {'coherence', 'unwrapped_phase', 'wrapped_angle', 'wrapped_amplitude', 'atm_delay_diff', 'dswe'}
        app.ddSARProd.Value = prod{1};
        for pi = testPairs
            app.ddSARPair.Value = pi;
            pols = ternary(strcmp(prod{1}, 'atm_delay_diff'), 1, [1 4]); % HH + VV where pol applies (dswe has pols)
            for po = pols
                app.ddSARPol.Value = po; app.cb.sar();
                n = n + 1; snap(sprintf('%02d_%s_sar_%s_p%02d_%s', n, char(ln), ...
                    prod{1}, pi, app.pol{po}));
            end
        end
    end
end

% ---- Lidar panel (heading-independent): every product, first + last date ----
for prod = {'sd', 'dtm', 'dsm', 'chm'}
    app.ddLidProd.Value = prod{1}; app.cb.lidProd();
    dates = cell2mat(app.ddLidDate.ItemsData);
    for ti = unique([dates(1) dates(end)])
        app.ddLidDate.Value = ti; app.cb.lid();
        n = n + 1; snap(sprintf('%02d_lid_%s_t%02d', n, prod{1}, ti));
    end
end

fprintf('mcs_gui_test: %d states exercised, screenshots in %s/\n', n, outDir);
close(app.fig);                                                  % clean shutdown = no lingering figure
end

function out = ternary(cond, a, b)
% Inline conditional (MATLAB has none).
if cond, out = a; else, out = b; end
end

function snapRetry(fig, path)
% exportapp occasionally fails or captures a stale frame in headless batch runs --
% flush the render pipeline first, and retry once on failure.
drawnow;
try
    exportapp(fig, path);
catch
    pause(2); drawnow;
    exportapp(fig, path);
end
end
