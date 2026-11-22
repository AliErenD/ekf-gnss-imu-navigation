%% Part C: Multi-rate EKF with bias estimation (no outlier protection yet)

%% Q: process noise covariance
% (3,3),(4,4): IMU white noise entering psi, V through Eqs (5)-(6) - see derivation
Q33 = (sigma_w*dt_imu)^2;
Q44 = (sigma_a*dt_imu)^2;

% (5,5),(6,6): bias random-walk intensity. Engineering assumption (NOT from
% ground truth): over the full mission the bias may wander by about half
% of its own measurement-noise magnitude. See report for justification.
DRIFT_FRACTION = 0.5;
T_mission = t_imu(end) - t_imu(1);
target_std_a = DRIFT_FRACTION * sigma_a;
target_std_w = DRIFT_FRACTION * sigma_w;
Q55 = target_std_a^2 / N;
Q66 = target_std_w^2 / N;

Q = diag([0, 0, Q33, Q44, Q55, Q66]);

fprintf('\n--- Part C: Q diagonal ---\n');
fprintf('Q = diag([%.3e %.3e %.3e %.3e %.3e %.3e])\n', diag(Q));

%% x0, P0
x0 = [pn0_hat; pe0_hat; psi_rwy; 0; ba_hat; bw_hat];

sigma_psi0 = deg2rad(2);
sigma_V0   = 0.2;
sigma_pn0  = sigma_n_R / sqrt(n_gnss_static);
sigma_pe0  = sigma_e_R / sqrt(n_gnss_static);
sigma_ba0  = sigma_a / sqrt(n_imu_static);
sigma_bw0  = sigma_w / sqrt(n_imu_static);

P0 = diag([sigma_pn0^2, sigma_pe0^2, sigma_psi0^2, sigma_V0^2, sigma_ba0^2, sigma_bw0^2]);

fprintf('\n--- Part C: x0 / P0 ---\n');
fprintf('x0 = [%.3f %.3f %.4f %.3f %.5f %.6f]\n', x0);
fprintf('sqrt(diag(P0)) = [%.3f %.3f %.5f %.3f %.6f %.7f]\n', sqrt(diag(P0)));

%% Run EKF without any gating (Part C(d))
out_C = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Q, R, false, inf);

%% Plot estimated trajectory vs GNSS and ground truth
figure('Name','C - EKF trajectory (no gating)');
plot(pe_true, pn_true, 'k-', 'LineWidth', 1.5); hold on;
plot(z_gnss(2,:), z_gnss(1,:), 'g.', 'MarkerSize', 8);
plot(out_C.X(2,:), out_C.X(1,:), 'b-');
grid on; axis equal;
xlabel('East [m]'); ylabel('North [m]');
legend('Ground truth','GNSS measurements','EKF estimate (ungated)','Location','best');
title('Part C: 6-state EKF, no outlier protection');
save_fig('fig_C1_ekf_ungated.png');

err_ekf_ungated = hypot(out_C.X(1,:)-pn_true, out_C.X(2,:)-pe_true);
fprintf('\n--- Part C: ungated EKF performance ---\n');
fprintf('RMSE (t>=20s) = %.3f m,  max error = %.3f m (at t=%.2f s)\n', ...
    sqrt(mean(err_ekf_ungated(t_imu>=20).^2)), max(err_ekf_ungated), t_imu(find(err_ekf_ungated==max(err_ekf_ungated),1)));
