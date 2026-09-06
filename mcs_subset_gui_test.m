function mcs_subset_gui_test(outDir)
%MCS_SUBSET_GUI_TEST Headless harness for mcs_subset_gui: drive every control state.
%   All pairs x UAVSAR products (HH + VV where pol applies) + both snow-depth
%   dates -> screenshots in outDir (default 'gui_test_png') for visual review.

if nargin < 1, outDir = 'gui_test_png'; end
if ~isfolder(outDir), mkdir(outDir); end

app = mcs_subset_gui();
n = 0;

for pi = cell2mat(app.ddPair.ItemsData)                          % all 6 pairs, A-F
    app.ddPair.Value = pi; app.cb.pair();
    for prod = {'dswe', 'coherence', 'unwrapped_phase', 'atm_delay_diff'}
        app.ddSARProd.Value = prod{1};
        pols = 1;                                                % HH; add VV for pol products
        if ~strcmp(prod{1}, 'atm_delay_diff'), pols = [1 4]; end
        for po = pols
            app.ddPol.Value = po; app.cb.sar();
            n = n + 1;
            snapRetry(app.fig, fullfile(outDir, sprintf('sub%02d_p%d_%s_%s.png', ...
                n, pi, prod{1}, app.pol{po})));
        end
    end
end
for ti = cell2mat(app.ddLidDate.ItemsData)                       % both QSI dates
    app.ddLidDate.Value = ti; app.cb.lid();
    n = n + 1;
    snapRetry(app.fig, fullfile(outDir, sprintf('sub%02d_sd_t%d.png', n, ti)));
end

fprintf('mcs_subset_gui_test: %d states exercised, screenshots in %s/\n', n, outDir);
close(app.fig);
end

function snapRetry(fig, path)
% exportapp occasionally fails or captures a stale frame headlessly -- flush + retry.
drawnow;
try
    exportapp(fig, path);
catch
    pause(2); drawnow;
    exportapp(fig, path);
end
end
