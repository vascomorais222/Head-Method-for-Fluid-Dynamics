%% Call Variables
close all
clear;clc
load('HEAD');

%% Initialize initial conditions
x_i = x_stn(1);
x_max = x_stn(5);
L = (x_max - x_i);
N = 1000;
step = L / N;

Ue_0 = Ue(4);
delta_asterisk = delta_asterisk_stn(1);
theta = theta_stn(1);
H = delta_asterisk / theta;
H1 = calculate_H1(H);
phi = H1 * Ue_0 * theta;

%% Setting initial conditions
F = zeros(2, N+1);
r0 = [theta; phi];
H1_vals = zeros(1, N+1); % Initialize H1_vals vector
H1_vals(1) = H1; % Store initial value of H1
F(:, 1) = r0;  % Save initial values of theta and phi

xi = x_i;
velocityValues = zeros(1, N+1);
dUe_dxValues = zeros(1, N+1);
x_values_L = zeros(1, N+1);

for i = 1:(N+2)
    [u, du_dx] = calculate_ExternalVelocity(xi);
    velocityValues(i) = u;
    dUe_dxValues(i) = du_dx;
    
    x_values_L(i) = (xi - x_i)/L;
    xi = xi + step;
end

%% Calculations

i = 2;  % Start from the second element
x_l = x_i;
while x_l <= x_max
    [Y1, Y2] = f(kin_visc, theta, H1_vals(i-1), velocityValues(i-1), dUe_dxValues(i-1));  % Use i-1 to access the previous values
    
    % Update variable to continue iterating
    theta = theta + step * Y1;
    phi = phi + step * Y2;
    
    H1_vals(i) = phi / (velocityValues(i-1) * theta);  % Calculate the corresponding H1 to feed function f
    
    % Store values in the matrix
    F(:, i) = [theta; phi];  % Increment index by 1 to store in the next column
    
    x_l = x_l + step;
    i = i + 1;
end

%% Determine missing flow quantities

H_head = calculate_H(H1_vals);
delta_asterisk_stn_head = H_head.*F(1,:);
delta_stn_head = (H1_vals .* F(1, :)) + delta_asterisk_stn_head;
cf_stn_Ludwieg = calculate_Cf(H1_vals, velocityValues, F(1, :), kin_visc);

%% Plotting
figure;
sgtitle("Head's Method - Flow Quantities along x/L");

subplot(2, 3, 1);
plot(x_values_L, F(1, :), 'b');
hold on;
plot((x_stn - x_i)/L, theta_stn, 'k*');  % Add experimental values
legend("\theta by Head's Method", "\theta experimental","location","southeast");
hold off;
xlabel('x/L');
ylabel('\theta');
title('\theta');

subplot(2, 3, 2);
plot(x_values_L, H_head, 'r');
hold on;
plot((x_stn - x_i)/L, H_stn, 'k*'); % Add experimental values
legend("H by Head's Method", "H experimental","location","southeast");
hold off;
xlabel('x/L');
ylabel('H');
title('H');
legend('Location', [0.5 0.7 0.08 0.1]); % Manually set legend position

subplot(2, 3, 3);
plot(x_values_L, delta_asterisk_stn_head, 'g');
hold on;
plot((x_stn - x_i)/L, delta_asterisk_stn, 'k*');  % Add experimental values
legend("\delta^* by Head's Method", "\delta^* experimental","location","southeast");
hold off;
xlabel('x/L');
ylabel('\delta^*');
title('\delta^*');

subplot(2, 3, 4);
plot(x_values_L, delta_stn_head, 'm');
hold on;
plot((x_stn - x_i)/L, deltaBL_stn(2, 4:8), 'k*');  % Add experimental values
legend("\delta by Head's Method", "\delta experimental","location","southeast");
hold off;
xlabel('x/L');
ylabel('\delta');
title('\delta');

subplot(2, 3, 5);
plot(x_values_L, cf_stn_Ludwieg, 'c');
hold on;
plot((x_stn - x_i)/L, Cf_stn(4:8), 'k*');  % Add experimental values
legend("c_f by Ludwieg-Tillman", "c_f experimental","location","northeast");
hold off;
xlabel('x/L');
ylabel('c_f');
title('c_f');

subplot(2, 3, 6);
plot(x_values_L, velocityValues, 'k');
hold on;
plot((x_stn - x_i)/L, Ue(4:8), 'k*');  % Add experimental values
legend("U_e Quadratic Fit", "U_e experimental","location","southeast");
hold off;
xlabel('x/L');
ylabel('U_e');
title('U_e');

figure;
plot(x_values_L, H1_vals, 'm');
legend("H1 by Head's Method","location","northeast");
xlabel('x/L');
ylabel('H_1');
title('H_1');

%% Auxiliary Functions

function [Y1, Y2] = f(niu, theta, H1, Ue, dUe_dx)
    cf = calculate_Cf(H1, Ue, theta, niu);
    
    Y1 = -theta * (calculate_H(H1) + 2) / Ue * dUe_dx + cf / 2;
    Y2 = Ue * calculate_FH1(H1);
end

function H1 = calculate_H1(H)
    if H <= 1.6
        H1 = 0.8234 .* (H - 1.1).^(-1.287) + 3.3;
    else
        H1 = 1.5501 .* (H - 0.6778).^(-3.064) + 3.3;
    end
end

function H = calculate_H(H1)
    if H1 >= 5.31
        H = 0.8599 .* (H1 - 3.3).^(-1.0 / 1.287) + 1.1;
    else
        H = 1.1538 .* (H1 - 3.3).^(-1.0 / 3.064) + 0.6778;
    end
end

function FH1 = calculate_FH1(H1)
    exponent = -0.6169;
    FH1 = 0.0306 .* (H1 - 3).^exponent;
end

function cf = calculate_Cf(H1, Ue, theta, niu)
    exponent = -0.678 .* calculate_H(H1);
    arg = Ue .* theta ./ niu;
    cf = 0.248 .* 10.^exponent .* arg.^(-0.268);
end

function [Ue, dUe_dx] = calculate_ExternalVelocity(xi)
%     Ue = -8.078 * xi^4 + 26.463 * xi^3 - 30.836 * xi^2 + 15.666 * xi + 23.517;
%     dUe_dx = -32.312 * xi^3 + 79.389 * xi^2 - 61.672 * xi + 15.666;
    Ue = -0.774 * xi^2 + 2.251*xi + 25.294;
    dUe_dx = -1.548*xi + 2.251;
    
end