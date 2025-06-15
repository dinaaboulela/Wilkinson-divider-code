clc;
clear all;
close all;

%% Given Parameters
f1 = 942.5e6;       % Frequency 1 (Hz)
f2 = 1945e6;        % Frequency 2 (Hz)
er = 2.2;           % Relative permittivity (RO5880)
h = 1.57e-3;        % Substrate height (m)
c = 3e8;            % Speed of light (m/s)
theta2_f1 = 58.75;  % Electrical length for Z2 at f1 (degrees)
theta3_f1 = 117.51; % Electrical length for Z3 at f1 (degrees)
Z0 = 50;            % Port impedance (Ohm)
Z2 = 42.91;         % Impedance of T-section branch (Ohm)
Z3 = 79.102;        % Impedance of open stub (Ohm)

%% 1. Calculate Microstrip Widths (W) for Z0, Z2, Z3
[W0, er_eff0] = microstrip_width(Z0, er, h);
[W2, er_eff2] = microstrip_width(Z2, er, h);
[W3, er_eff3] = microstrip_width(Z3, er, h);

fprintf('Microstrip Widths:\n');
fprintf('  Z0 = 50 Ohm:   W0 = %.4f mm, ε_eff = %.4f\n', W0*1e3, er_eff0);
fprintf('  Z2 = 42.91 Ohm: W2 = %.4f mm, ε_eff = %.4f\n', W2*1e3, er_eff2);
fprintf('  Z3 = 79.1 Ohm:  W3 = %.4f mm, ε_eff = %.4f\n\n', W3*1e3, er_eff3);

%% 2. Calculate Physical Lengths (L) for Z2 and Z3 at f1
lambda0_f1 = c / f1;       % Wavelength at f1

% Length for Z2 (θ2 = 58.75°)
lambda_g2 = lambda0_f1 / sqrt(er_eff2);
L2 = (theta2_f1 / 360) * lambda_g2;

% Length for Z3 (θ3 = 117.51°)
lambda_g3 = lambda0_f1 / sqrt(er_eff3);
L3 = (theta3_f1 / 360) * lambda_g3;

fprintf('Branch Lengths at f1 (%.2f MHz):\n', f1/1e6);
fprintf('  L2 (Z2): %.4f mm (θ = %.2f°)\n', L2*1e3, theta2_f1);
fprintf('  L3 (Z3): %.4f mm (θ = %.2f°)\n\n', L3*1e3, theta3_f1);

%% 3. Dual-Band Input/Output Line (50 Ω) Optimization
f_center = sqrt(f1 * f2);  % Geometric mean frequency
lambda0_center = c / f_center;
lambda_g50_center = lambda0_center / sqrt(er_eff0);
L50_quarterwave = lambda_g50_center / 4;  % 90° at f_center

% Option 2: Dual-band length (90° at f1, ~180° at f2)
L50_dualband = (90/360) * (c / f1) / sqrt(er_eff0);
theta_f2 = (f2 / f1) * 90;  % Phase shift at f2

fprintf('Dual-Band 50 Ω Line Options:\n');
fprintf('  Option 1: Quarter-wave at f_center (%.2f MHz) = %.4f mm\n', f_center/1e6, L50_quarterwave*1e3);
fprintf('  Option 2: 90° at f1 (%.2f MHz) = %.4f mm (θ = %.2f° at f2)\n\n', f1/1e6, L50_dualband*1e3, theta_f2);

%% 4. Verification at Both Frequencies
% Check electrical lengths at f2
theta2_f2 = (f2 / f1) * theta2_f1;
theta3_f2 = (f2 / f1) * theta3_f1;

fprintf('Electrical Lengths at f2 (%.2f MHz):\n', f2/1e6);
fprintf('  θ2 (Z2): %.2f°\n', theta2_f2);
fprintf('  θ3 (Z3): %.2f°\n', theta3_f2);

%% Microstrip Width Calculation Function
function [W, er_eff] = microstrip_width(Z, er, h)
    % Hammerstad & Jensen approximation for W/h
    A = (Z / 60) * sqrt((er + 1)/2) + (er - 1)/(er + 1)*(0.23 + 0.11/er);
    B = 377 * pi / (2 * Z * sqrt(er));
    
    if A > 1.52
        W_over_h = (8 * exp(A)) / (exp(2*A) - 2);
    else
        W_over_h = (2/pi) * (B - 1 - log(2*B - 1) + (er - 1)/(2*er)*(log(B - 1) + 0.39 - 0.61/er));
    end
    
    W = W_over_h * h;  % Width in meters
    
    % Effective permittivity (ε_eff)
    er_eff = (er + 1)/2 + (er - 1)/2 * (1 + 12 * h/W)^(-0.5);
end
