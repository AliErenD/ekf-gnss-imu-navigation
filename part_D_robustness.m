%% Part D: Robustness - Mahalanobis gate, GNSS outage, 4-state comparison

%% (a) Gate threshold
P_CONF = 0.99;
gamma = -2*log(1 - P_CONF);
fprintf('\n--- Part D(a): Mahalanobis gate ---\n');
fprintf('Confidence level p = %.3f  ->  gamma = %.4f (chi2, 2 dof)\n', P_CONF, gamma);

out_Dgated = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Q, R, true, gamma);

rejected_idx = find(~out_Dgated.accepted);
fprintf('Rejected epochs: %d / %d\n', numel(rejected_idx), M);
fprintf('%8s %10s\n', 't [s]', 'd^2');
for ii = rejected_idx
    fprintf('%8.3f %10.3f\n', t_gnss(ii), out_Dgated.d2(ii));
end

%% Classify rejected epochs as true outlier vs false alarm (TRUTH USED FOR EVALUATION ONLY)
k_gnss = out_Dgated.k_gnss;
true_outlier = false(1, M);
for ii = rejected_idx
    k = k_gnss(ii);
    meas_err = z_gnss(:,ii) - [pn_true(k); pe_true(k)];
    d2_meas = meas_err' * (R \ meas_err);  % is the raw measurement itself inconsistent with R?
    true_outlier(ii) = d2_meas > gamma;
end
false_alarm_idx = rejected_idx(~true_outlier(rejected_idx));
true_outlier_idx = rejected_idx(true_outlier(rejected_idx));

fprintf('\nClassification: %d true outliers, %d false alarms\n', numel(true_outlier_idx), numel(false_alarm_idx));
n_accepted = sum(out_Dgated.accepted);
fprintf('Observed false-alarm rate among all epochs = %.2f%%  (expected ~%.2f%% at p=%.3f)\n', ...
    100*numel(false_alarm_idx)/M, 100*(1-P_CONF), P_CONF);

%% (b) Error + 3-sigma bound, outage window, P/K behaviour
err_gated = hypot(out_Dgated.X(1,:) - pn_true, out_Dgated.X(2,:) - pe_true);
sigma3h = 3*sqrt(out_Dgated.Pdiag(1,:) + out_Dgated.Pdiag(2,:));

figure('Name','D - error and 3-sigma bound');
plot(t_imu, err_gated, 'b-'); hold on;
plot(t_imu, sigma3h, 'r--');
xline(t_outage_start,'k:'); xline(t_outage_end,'k:');
grid on;
xlabel('t [s]'); ylabel('horizontal error [m]');
legend('|error|','3\sigma_h bound','outage start/end','Location','best');
title('Part D: gated 6-state EKF error vs 3\sigma bound');
save_fig('fig_D1_error_3sigma.png');

outage_mask = t_imu >= t_outage_start & t_imu <= t_outage_end;
[max_err_outage, idx_max] = max(err_gated(outage_mask));
t_outage_vec = t_imu(outage_mask);
fprintf('\n--- Part D(b): outage behaviour ---\n');
fprintf('Max horizontal error during outage = %.3f m (at t=%.2f s)\n', max_err_outage, t_outage_vec(idx_max));

% Last fix before the gap and first fix after it (from Part A, timestamps only)
j_before = outage_idx(1);
j_after  = outage_idx(1) + 1;
Ppos_post_before = sum(out_Dgated.Pdiag(1:2, k_gnss(j_before)));  % P_k|k after the last update
Ppos_prior_after = out_Dgated.Ppos_prior(j_after);                % P_k|k-1 before the first new update
Ppos_post_after  = sum(out_Dgated.Pdiag(1:2, k_gnss(j_after)));   % P_k|k after the first new update
fprintf('P11+P22: after last pre-outage update = %.3f m^2 (3sigma_h = %.2f m)\n', Ppos_post_before, 3*sqrt(Ppos_post_before));
fprintf('         prior at first post-outage fix = %.3f m^2 (3sigma_h = %.2f m), growth %.1fx\n', ...
    Ppos_prior_after, 3*sqrt(Ppos_prior_after), Ppos_prior_after/Ppos_post_before);
fprintf('         after first post-outage update = %.3f m^2 (3sigma_h = %.2f m)\n', Ppos_post_after, 3*sqrt(Ppos_post_after));

K_before = out_Dgated.K_hist(:,:,j_before);
K_after  = out_Dgated.K_hist(:,:,j_after);
fprintf('Gain applied at last pre-outage fix : K11 = %.3f, K22 = %.3f, ||K||_F = %.3f\n', K_before(1,1), K_before(2,2), norm(K_before,'fro'));
fprintf('Gain applied at first post-outage fix: K11 = %.3f, K22 = %.3f, ||K||_F = %.3f\n', K_after(1,1), K_after(2,2), norm(K_after,'fro'));
fprintf('d^2 of first post-outage fix = %.3f (gamma = %.3f)\n', out_Dgated.d2(j_after), gamma);

%% (c) 4-state EKF with fixed static bias correction
Q4 = diag([0, 0, Q33, Q44]);
x0_4 = [pn0_hat; pe0_hat; psi_rwy; 0];
P0_4 = diag([sigma_pn0^2, sigma_pe0^2, sigma_psi0^2, sigma_V0^2]);

out_4_nogate = ekf4_run(t_imu, u_imu, t_gnss, z_gnss, x0_4, P0_4, Q4, R, ba_hat, bw_hat, false, inf);
err_4_nogate = hypot(out_4_nogate.X(1,:) - pn_true, out_4_nogate.X(2,:) - pe_true);

fprintf('\n--- Part D(c): 4-state vs 6-state, no gating, outage window ---\n');
fprintf('Max error during outage: 6-state = %.3f m, 4-state = %.3f m\n', ...
    max(err_gated(outage_mask)), max(err_4_nogate(outage_mask)));
% (use ungated 6-state for a fair apples-to-apples comparison)
err_6_nogate = hypot(out_C.X(1,:) - pn_true, out_C.X(2,:) - pe_true);
fprintf('(ungated 6-state max error during outage = %.3f m, for reference)\n', max(err_6_nogate(outage_mask)));

%% (c) continued: enable the gate on the 4-state filter
out_4_gated = ekf4_run(t_imu, u_imu, t_gnss, z_gnss, x0_4, P0_4, Q4, R, ba_hat, bw_hat, true, gamma);
err_4_gated = hypot(out_4_gated.X(1,:) - pn_true, out_4_gated.X(2,:) - pe_true);

fprintf('\n--- Part D(c): 4-state WITH gate ---\n');
fprintf('Total rejected epochs (4-state, gated) = %d / %d\n', sum(~out_4_gated.accepted), M);
fprintf('Max consecutive rejections = %d\n', max(out_4_gated.consec_rejected));
fprintf('Final horizontal error: 4-state gated = %.3f m, 6-state gated = %.3f m\n', err_4_gated(end), err_gated(end));

sigma3h_4 = 3*sqrt(out_4_gated.Pdiag(1,:) + out_4_gated.Pdiag(2,:));
fprintf('4-state 3sigma_h: max = %.2f m, at t=%.0f s = %.2f m (actual error there %.1f m)\n', ...
    max(sigma3h_4(t_imu >= 20)), t_imu(end), sigma3h_4(end), err_4_gated(end));

figure('Name','D - 4-state vs 6-state comparison');
semilogy(t_imu, err_gated, 'b-'); hold on;
semilogy(t_imu, err_4_gated, 'm-');
semilogy(t_imu, sigma3h, 'b--');
semilogy(t_imu, sigma3h_4, 'm--');
xline(t_outage_start,'k:'); xline(t_outage_end,'k:');
grid on;
xlabel('t [s]'); ylabel('horizontal error [m]');
legend('6-state |error|','4-state |error|','6-state 3\sigma_h','4-state 3\sigma_h','outage','Location','northwest');
title('Part D(c): 6-state vs 4-state EKF, both gated');
save_fig('fig_D2_4state_vs_6state.png');

k_leave = find(t_imu >= 20 & err_4_gated > sigma3h_4, 1);
fprintf('4-state error first exceeds its own 3sigma_h at t = %.1f s\n', t_imu(k_leave));
