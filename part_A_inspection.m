%% Part A(a): raw data plot, sampling rates, async proof, GNSS outage
% Called from main_p1.m. Uses t_imu, u_imu, t_gnss, z_gnss already loaded.

%% Raw measurement plots (reproduces Fig. 1 style)
figure('Name','A - Raw measurements');

subplot(2,2,[1 3]);
plot(z_gnss(2,:), z_gnss(1,:), '.'); grid on; axis equal;
xlabel('East [m]'); ylabel('North [m]'); title('(a) GNSS position measurements');

subplot(2,2,2);
plot(t_imu, u_imu(1,:)); grid on;
ylabel('a_m [m/s^2]'); title('(b) IMU measurements');

subplot(2,2,4);
plot(t_imu, rad2deg(u_imu(2,:))); grid on;
xlabel('t [s]'); ylabel('\omega_m [deg/s]');

save_fig('fig_A1_raw_measurements.png');

%% Sampling rates
dt_imu_vec  = diff(t_imu);
dt_gnss_vec = diff(t_gnss);

fs_imu  = 1 / mean(dt_imu_vec);
fs_gnss = 1 / mean(dt_gnss_vec);   % average rate only - GNSS is NOT uniform

fprintf('IMU : fs = %.3f Hz (dt = %.4f s, std = %.3e s)\n', fs_imu, mean(dt_imu_vec), std(dt_imu_vec));
fprintf('GNSS: average fs = %.3f Hz (mean dt = %.4f s, std = %.4f s)\n', fs_gnss, mean(dt_gnss_vec), std(dt_gnss_vec));

%% Proof that GNSS timestamps are not on the IMU grid
% If t_gnss(j) fell exactly on the IMU grid, (t_gnss(j) - t_imu(1)) would
% be an integer multiple of dt_imu, i.e. remainder ~ 0 (mod dt_imu).
dt_imu_nom = dt_imu; % 0.01 s nominal, set in main_p1.m
offset_from_grid = mod(t_gnss - t_imu(1), dt_imu_nom);
% Fold into [-dt_imu_nom/2, dt_imu_nom/2] for a symmetric "distance to nearest grid point"
offset_from_grid(offset_from_grid > dt_imu_nom/2) = offset_from_grid(offset_from_grid > dt_imu_nom/2) - dt_imu_nom;

fprintf('GNSS offset from nearest IMU sample: mean=%.5f s, max|.|=%.5f s (IMU dt = %.3f s)\n', ...
    mean(offset_from_grid), max(abs(offset_from_grid)), dt_imu_nom);

figure('Name','A - GNSS/IMU timing');
stem(t_gnss, offset_from_grid*1000, 'filled'); grid on;
xlabel('t [s]'); ylabel('offset from nearest IMU sample [ms]');
title('GNSS time stamps are not aligned with the IMU sampling grid');
save_fig('fig_A2_timing_offset.png');

%% GNSS outage detection (timestamps only - no position/truth used)
gap = diff(t_gnss);
nominal_gap = median(gap);              % robust estimate of the "normal" 1 Hz spacing
outage_idx  = find(gap > 5*nominal_gap); % anything many times the normal spacing = outage

if ~isempty(outage_idx)
    t_outage_start = t_gnss(outage_idx(1));
    t_outage_end   = t_gnss(outage_idx(1)+1);
    outage_duration = t_outage_end - t_outage_start;
    fprintf('GNSS outage detected: start = %.3f s, end = %.3f s, duration = %.3f s\n', ...
        t_outage_start, t_outage_end, outage_duration);
else
    fprintf('No GNSS outage detected above threshold.\n');
end
