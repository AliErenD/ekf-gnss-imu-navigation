%% Kalman Filter-based Navigation and Sensor Fusion
% Main script. Run this file from start to end; it loads the data,
% produces every figure/number used in the report, section by section.
%
% Student: Ali Eren Dağ
% Date: 16.10.2026

clear; close all; clc;

%% Load data
data = load('navigation_data.mat');

t_imu  = data.t_imu;     % 1xN [s]
u_imu  = data.u_imu;     % 2xN [am; wm]  [m/s^2; rad/s]
t_gnss = data.t_gnss;    % 1xM [s]
z_gnss = data.z_gnss;    % 2xM [pn; pe]  [m]

% Ground truth - EVALUATION / PLOTTING ONLY. Never feed into the filter.
pn_true  = data.pn_true;
pe_true  = data.pe_true;
psi_true = data.psi_true;
V_true   = data.V_true;
ba_true  = data.ba_true;
bw_true  = data.bw_true;

N = numel(t_imu);
M = numel(t_gnss);
dt_imu = t_imu(2) - t_imu(1);   % nominal IMU sample period

fprintf('=== Part A: Data inspection ===\n');

%% ---- A(a): raw data plot, sampling rates, async proof, outage detection ----
part_A_inspection;

%% ---- A(b): static calibration (stationary segment, t < 20 s) ----
part_A_calibration;

fprintf('\n=== Part B: Dead reckoning ===\n');
part_B_deadreckoning;

fprintf('\n=== Part C: Multi-rate EKF (no gating) ===\n');
part_C_ekf;

fprintf('\n=== Part D: Robustness (gate, outage, 4-state) ===\n');
part_D_robustness;

fprintf('\n=== Part E: Consistency (NEES/NIS) and Q tuning ===\n');
part_E_consistency;

fprintf('\n=== Part G: Design checks (R, gate level, bias random walk) ===\n');
part_G_design_checks;

fprintf('\n=== Part F: Results summary ===\n');
part_F_summary;
