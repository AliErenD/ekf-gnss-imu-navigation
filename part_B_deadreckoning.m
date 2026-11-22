%% Part B: Dead reckoning with IMU only
% Uses pn0_hat, pe0_hat, ba_hat, bw_hat computed in part_A_calibration.m
% psi0 and V0 are given by the scenario, NOT estimated from data.

psi_rwy = deg2rad(30);  % given runway heading
V0  = 0;                % stationary at t=0
am  = u_imu(1,:);
wm  = u_imu(2,:);

%% (i) Raw IMU (biases assumed zero)
[pn_raw, pe_raw, psi_raw, V_raw] = dead_reckon(am, wm, dt_imu, pn0_hat, pe0_hat, psi_rwy, V0, 0, 0);

%% (ii) Bias-corrected IMU (static estimates from Part A)
[pn_cal, pe_cal, psi_cal, V_cal] = dead_reckon(am, wm, dt_imu, pn0_hat, pe0_hat, psi_rwy, V0, ba_hat, bw_hat);

%% Horizontal position error vs ground truth (truth used ONLY for plotting/evaluation)
err_raw = hypot(pn_raw - pn_true, pe_raw - pe_true);
err_cal = hypot(pn_cal - pn_true, pe_cal - pe_true);

%% Plot: trajectories
figure('Name','B - Dead reckoning trajectories');
plot(pe_true, pn_true, 'k-', 'LineWidth', 1.5); hold on;
plot(pe_raw, pn_raw, 'r--');
plot(pe_cal, pn_cal, 'b-.');
grid on; axis equal;
xlabel('East [m]'); ylabel('North [m]');
legend('Ground truth','Raw IMU DR','Bias-corrected DR','Location','best');
title('Part B: dead-reckoned trajectories vs ground truth');
save_fig('fig_B1_trajectories.png');

%% Plot: horizontal error vs time (log scale)
figure('Name','B - Dead reckoning error');
semilogy(t_imu, err_raw, 'r-'); hold on;
semilogy(t_imu, err_cal, 'b-');
grid on;
xlabel('t [s]'); ylabel('horizontal position error [m]');
legend('Raw IMU DR','Bias-corrected DR','Location','best');
title('Part B: dead-reckoning position error (log scale)');
save_fig('fig_B2_error.png');

%% Empirical growth-rate check: fit log(err) = p*log(t) + c on t in [30,280] s
% (avoid the stationary segment and the very end where error saturates)
fit_mask = t_imu > 30 & t_imu < 280;
p_raw = polyfit(log(t_imu(fit_mask)), log(err_raw(fit_mask)), 1);
p_cal = polyfit(log(t_imu(fit_mask)), log(err_cal(fit_mask)), 1);

fprintf('\n--- Part B: empirical error growth exponent (error ~ t^p) ---\n');
fprintf('Raw IMU DR        : p = %.2f  (theory for ~constant large bias: p ~ 2)\n', p_raw(1));
fprintf('Bias-corrected DR : p = %.2f  (theory for drifting residual bias: p > 2)\n', p_cal(1));
fprintf('Final error:  raw = %.1f m,  calibrated = %.1f m  (t = %.0f s)\n', err_raw(end), err_cal(end), t_imu(end));
fprintf('Final heading error: raw = %.1f deg, calibrated = %.1f deg\n', ...
    rad2deg(wrap_angle(psi_raw(end) - psi_true(end))), rad2deg(wrap_angle(psi_cal(end) - psi_true(end))));
dr_mask = t_imu >= 20;
fprintf('Horizontal RMSE (t>=20s): raw = %.1f m, calibrated = %.1f m\n', ...
    sqrt(mean(err_raw(dr_mask).^2)), sqrt(mean(err_cal(dr_mask).^2)));

%% --- local function ---
function [pn, pe, psi, V] = dead_reckon(am, wm, dt, pn0, pe0, psi0, V0, ba, bw)
    N = numel(am);
    psi = psi0 + [0, cumsum((wm(1:end-1) - bw) * dt)];
    V   = V0   + [0, cumsum((am(1:end-1) - ba) * dt)];
    pn  = pn0  + [0, cumsum(V(1:end-1) .* cos(psi(1:end-1)) * dt)];
    pe  = pe0  + [0, cumsum(V(1:end-1) .* sin(psi(1:end-1)) * dt)];
end
