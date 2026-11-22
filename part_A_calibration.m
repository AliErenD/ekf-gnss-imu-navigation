%% Part A(b): static calibration from the stationary segment (t < 20 s)
% Stationary -> true a = 0, true omega = 0, true position constant.
% So am = ba + wa and wm = bw + ww directly during this window:
% the SAMPLE MEAN estimates the bias, the SAMPLE STD estimates the noise.

T_STATIC = 20; % [s], given by the scenario

imu_static_mask  = t_imu  < T_STATIC;
gnss_static_mask = t_gnss < T_STATIC;

am_static = u_imu(1, imu_static_mask);
wm_static = u_imu(2, imu_static_mask);
n_imu_static = sum(imu_static_mask);

%% IMU: bias and noise std estimates
ba_hat = mean(am_static);
bw_hat = mean(wm_static);
sigma_a = std(am_static);
sigma_w = std(wm_static);

fprintf('\n--- IMU static calibration (N = %d samples, %.1f s) ---\n', n_imu_static, T_STATIC);
fprintf('ba_hat  = %.5f m/s^2   (sigma_a = %.5f m/s^2)\n', ba_hat, sigma_a);
fprintf('bw_hat  = %.6f rad/s = %.4f deg/s  (sigma_w = %.6f rad/s = %.4f deg/s)\n', ...
    bw_hat, rad2deg(bw_hat), sigma_w, rad2deg(sigma_w));

%% GNSS: position std and initial position estimate
pn_static = z_gnss(1, gnss_static_mask);
pe_static = z_gnss(2, gnss_static_mask);
n_gnss_static = sum(gnss_static_mask);

sigma_n = std(pn_static);
sigma_e = std(pe_static);
pn0_hat = mean(pn_static);
pe0_hat = mean(pe_static);

fprintf('\n--- GNSS static calibration (M = %d samples, %.1f s) ---\n', n_gnss_static, T_STATIC);
fprintf('sigma_n = %.3f m, sigma_e = %.3f m\n', sigma_n, sigma_e);
fprintf('Initial position estimate: pn0 = %.3f m, pe0 = %.3f m\n', pn0_hat, pe0_hat);

%% Reliability of the std estimates: chi-square confidence interval
% For n iid Gaussian samples, (n-1)*s^2/sigma^2 ~ chi2(n-1).
% => 95% CI for sigma^2: [ (n-1)s^2 / chi2inv(0.975,n-1), (n-1)s^2 / chi2inv(0.025,n-1) ]
alpha = 0.05;

ci_sigma = @(s, n) sqrt((n-1)*s^2 ./ [chi2inv(1-alpha/2, n-1), chi2inv(alpha/2, n-1)]);

ci_a = ci_sigma(sigma_a, n_imu_static);
ci_w = ci_sigma(sigma_w, n_imu_static);
ci_n = ci_sigma(sigma_n, n_gnss_static);
ci_e = ci_sigma(sigma_e, n_gnss_static);

fprintf('\n--- 95%% confidence intervals for the std estimates ---\n');
fprintf('sigma_a in [%.5f, %.5f] m/s^2   (n=%d)\n', ci_a(1), ci_a(2), n_imu_static);
fprintf('sigma_w in [%.6f, %.6f] rad/s   (n=%d)\n', ci_w(1), ci_w(2), n_imu_static);
fprintf('sigma_n in [%.3f, %.3f] m        (n=%d)  <- WIDE: only %d GNSS samples in 20 s\n', ci_n(1), ci_n(2), n_gnss_static, n_gnss_static);
fprintf('sigma_e in [%.3f, %.3f] m        (n=%d)  <- WIDE: only %d GNSS samples in 20 s\n', ci_e(1), ci_e(2), n_gnss_static, n_gnss_static);

%% Build R with a safety margin
% The IMU estimates rest on ~2000 samples -> tight CI, use the point estimate directly.
% The GNSS estimates rest on only ~18-20 samples -> CI is wide (see print above):
% using the upper 95% CI bound instead of the point estimate protects against
% underestimating R (which would make the filter overconfident in GNSS).
SAFETY_MARGIN = true;
if SAFETY_MARGIN
    sigma_n_R = ci_n(2);
    sigma_e_R = ci_e(2);
else
    sigma_n_R = sigma_n;
    sigma_e_R = sigma_e;
end

R = diag([sigma_n_R^2, sigma_e_R^2]);
fprintf('\nR used in the filter = diag([%.3f, %.3f])  (m^2)\n', R(1,1), R(2,2));
