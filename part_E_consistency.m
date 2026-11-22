%% Part E: Consistency (NEES/NIS) and Q tuning
% Baseline filter = out_Dgated from Part D (gated, 6-state, Q as built in Part C).

state_names = {'pn [m]','pe [m]','\psi [rad]','V [m/s]','b_a [m/s^2]','b_\omega [rad/s]'};
x_true_all = [pn_true; pe_true; psi_true; V_true; ba_true; bw_true];

eval_mask = t_imu >= 20; % exclude the stationary segment, per the hint

%% (a) Error +- 3-sigma for all six states
err_all = out_Dgated.X - x_true_all;
err_all(3,:) = wrap_angle(err_all(3,:));

figure('Name','E - state errors with 3-sigma bounds', 'Position',[100 100 900 900]);
for i = 1:6
    subplot(3,2,i);
    sig3 = 3*sqrt(out_Dgated.Pdiag(i,:));
    plot(t_imu, err_all(i,:), 'b-'); hold on;
    plot(t_imu, sig3, 'r--'); plot(t_imu, -sig3, 'r--');
    xline(t_outage_start,'k:'); xline(t_outage_end,'k:');
    grid on; xlabel('t [s]'); ylabel(state_names{i});
    if i == 1
        legend('error','\pm3\sigma','','outage','','Location','southwest');
    end
end
sgtitle('Part E(a): estimation error vs \pm3\sigma bounds (gated 6-state EKF)');
save_fig('fig_E1_errors_3sigma.png');

%% (b) NEES and NIS for the baseline filter
[nees_base, nis_mean_base, nis_frac_base, ~] = compute_consistency(out_Dgated, x_true_all, eval_mask, t_gnss, t_imu);

fprintf('\n--- Part E(b): consistency of the baseline (gated, Q scale=1) filter ---\n');
fprintf('Mean NEES = %.3f  (expected 6, acceptable range ~[3,12])\n', mean(nees_base(eval_mask)));
fprintf('Mean NIS  = %.3f  (expected 2)\n', nis_mean_base);
fprintf('Fraction of NIS inside two-sided 95%% interval [0.051,7.38] = %.1f%%\n', 100*nis_frac_base);

%% (c) Q scaling sweep: 0.01x, 1x, 100x
scales = [0.01, 1, 100];
results = struct('scale',{},'rmse',{},'max_outage_err',{},'mean_nis',{},'mean_nees',{},'n_rejected',{});

for si = 1:numel(scales)
    s = scales(si);
    Qs = Q * s;
    out_s = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Qs, R, true, gamma);
    err_s = hypot(out_s.X(1,:)-pn_true, out_s.X(2,:)-pe_true);
    rmse_s = sqrt(mean(err_s(eval_mask).^2));
    max_out_s = max(err_s(outage_mask));
    [nees_s, nis_mean_s, ~, n_rej_s] = compute_consistency(out_s, x_true_all, eval_mask, t_gnss, t_imu);

    results(si).scale = s;
    results(si).rmse = rmse_s;
    results(si).max_outage_err = max_out_s;
    results(si).mean_nis = nis_mean_s;
    results(si).mean_nees = mean(nees_s(eval_mask));
    results(si).n_rejected = n_rej_s;
end

fprintf('\n--- Part E(c): Q scaling sweep ---\n');
fprintf('%8s %10s %14s %10s %10s %10s\n','Q scale','RMSE [m]','MaxOutErr[m]','MeanNIS','MeanNEES','#Rejected');
for si = 1:numel(scales)
    r = results(si);
    fprintf('%8.2f %10.3f %14.3f %10.3f %10.3f %10d\n', r.scale, r.rmse, r.max_outage_err, r.mean_nis, r.mean_nees, r.n_rejected);
end

%% --- local function ---
function [nees, nis_mean, nis_frac95, n_rej] = compute_consistency(out, x_true_all, eval_mask, t_gnss, t_imu)
    err = out.X - x_true_all;
    err(3,:) = wrap_angle(err(3,:));
    N = size(err,2);
    nees = nan(1,N);
    for k = find(eval_mask)
        Pk = out.Pfull(:,:,k);
        e = err(:,k);
        nees(k) = e' * (Pk \ e);
    end
    nis_eval_mask_gnss = t_gnss >= 20;
    nis_vals = out.d2(out.accepted & nis_eval_mask_gnss);
    nis_mean = mean(nis_vals);
    lo = chi2inv(0.025,2); hi = chi2inv(0.975,2);
    nis_frac95 = mean(nis_vals >= lo & nis_vals <= hi);
    n_rej = sum(~out.accepted);
end
