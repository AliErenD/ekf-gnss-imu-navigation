%% Part F: Results summary (Table 1), final (gated, 6-state) filter

rmse_final = sqrt(mean(err_gated(eval_mask).^2));

heading_err_final = wrap_angle(out_Dgated.X(3,end) - psi_true(end));

ba_final = out_Dgated.X(5,end); bw_final = out_Dgated.X(6,end);
ba_final_true = ba_true(end); bw_final_true = bw_true(end);

max_err_outage_4state_gated = max(err_4_gated(outage_mask));

fprintf('\n============================================================\n');
fprintf(' TABLE 1: Results summary (final gated 6-state EKF)\n');
fprintf('============================================================\n');
fprintf('Estimated sigma_n, sigma_e (static)      : %.3f m, %.3f m\n', sigma_n, sigma_e);
fprintf('Values used in R (sqrt of diag)           : %.3f m, %.3f m\n', sqrt(R(1,1)), sqrt(R(2,2)));
fprintf('GNSS outage: start / end / duration       : %.3f s / %.3f s / %.3f s\n', t_outage_start, t_outage_end, outage_duration);
fprintf('Horizontal RMSE (t>=20s), gated 6-state   : %.3f m\n', rmse_final);
fprintf('Max horizontal error during outage, 6-st. : %.3f m\n', max(err_gated(outage_mask)));
fprintf('Max horizontal error during outage, 4-st. : %.3f m (gated) / %.3f m (ungated)\n', ...
    max_err_outage_4state_gated, max(err_4_nogate(outage_mask)));
fprintf('Final heading error                       : %.4f rad = %.3f deg\n', heading_err_final, rad2deg(heading_err_final));
fprintf('Final b_a_hat (truth)                      : %.5f (%.5f) m/s^2\n', ba_final, ba_final_true);
fprintf('Final b_omega_hat (truth)                  : %.6f (%.6f) rad/s\n', bw_final, bw_final_true);
fprintf('Rejected GNSS epochs (true outlier/false)  : %d true outliers / %d false alarms (of %d total)\n', ...
    numel(true_outlier_idx), numel(false_alarm_idx), M);
fprintf('Mean NIS / mean NEES                       : %.3f / %.3f\n', nis_mean_base, mean(nees_base(eval_mask)));
fprintf('============================================================\n');
