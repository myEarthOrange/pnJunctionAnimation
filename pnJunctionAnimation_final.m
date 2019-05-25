% Homework Animation for Semiconductor Physics.
% Developed by Dali Cheng, Dept. of EE, Tsinghua University, in May 2019.
% Tested on MATLAB R2016b.
% Conceptual illustration of p-n junction formation, thermal equilibirum
% conditions and non-equilibrium conditions. Physical accuracy may not be
% expected.
% For academic usage only. E-mail to cheng-dl16@mails.tsinghua.edu.cn.
clear all; close all; clc;
h = figure;
h.Position = [100 100 1000 500];
h.Color = 'w';

pause(10);

%% Initialize
pause_time = 0.01;
pause_interval = 4;
step = 0.008;
p_type_x_min = -4; p_type_x_max = -1; p_type_y_min = 1; p_type_y_max = 3;
n_type_x_min = 1;  n_type_x_max = 4;  n_type_y_min = 1; n_type_y_max = 3;
% p-type semiconductor
x_p = [p_type_x_max p_type_x_max p_type_x_min p_type_x_min];
y_p = [p_type_y_min p_type_y_max p_type_y_max p_type_y_min];
c_p = [0.9  1  0.9];

% n-type semiconductor
x_n = [n_type_x_min n_type_x_min n_type_x_max n_type_x_max];
y_n = [n_type_y_min n_type_y_max n_type_y_max n_type_y_min];
c_n = [0.9  1  0.9];

num_row_carrier = 20;
num_col_carrier = 40;
num_row_ion = 5;
num_col_ion = 10;

% holes
temp_x_hc = x_p(3) + (x_p(1) - x_p(3)) * linspace(0, 1, num_col_carrier);
temp_y_hc = y_p(1) + (y_p(2) - y_p(1)) * linspace(0, 1, num_row_carrier);
[x_hc, y_hc] = meshgrid(temp_x_hc, temp_y_hc);
x_hc = reshape(x_hc, [1 num_row_carrier*num_col_carrier]);
y_hc = reshape(y_hc, [1 num_row_carrier*num_col_carrier]);

% electrons
temp_x_ec = x_n(1) + (x_n(3) - x_n(1)) * linspace(0, 1, num_col_carrier);
temp_y_ec = y_n(1) + (y_n(2) - y_n(1)) * linspace(0, 1, num_row_carrier);
[x_ec, y_ec] = meshgrid(temp_x_ec, temp_y_ec);
x_ec = reshape(x_ec, [1 num_row_carrier*num_col_carrier]);
y_ec = reshape(y_ec, [1 num_row_carrier*num_col_carrier]);

% acceptors
temp_x_ac = x_p(3) + (x_p(1) - x_p(3)) * linspace(0, 1, num_col_ion + 1);
temp_x_ac = temp_x_ac + (temp_x_ac(2) - temp_x_ac(1)) / 2;
temp_y_ac = y_p(1) + (y_p(2) - y_p(1)) * linspace(0, 1, num_row_ion + 1);
temp_y_ac = temp_y_ac + (temp_y_ac(2) - temp_y_ac(1)) / 2;
[x_ac, y_ac] = meshgrid(temp_x_ac(1:end-1), temp_y_ac(1:end-1));
x_ac = reshape(x_ac, [1 num_row_ion*num_col_ion]);
y_ac = reshape(y_ac, [1 num_row_ion*num_col_ion]);

% donors
temp_x_dc = x_n(1) + (x_n(3) - x_n(1)) * linspace(0, 1, num_col_ion + 1);
temp_x_dc = temp_x_dc + (temp_x_dc(2) - temp_x_dc(1)) / 2;
temp_y_dc = y_n(1) + (y_n(2) - y_n(1)) * linspace(0, 1, num_row_ion + 1);
temp_y_dc = temp_y_dc + (temp_y_dc(2) - temp_y_dc(1)) / 2;
[x_dc, y_dc] = meshgrid(temp_x_dc(1:end-1), temp_y_dc(1:end-1));
x_dc = reshape(x_dc, [1 num_row_ion*num_col_ion]);
y_dc = reshape(y_dc, [1 num_row_ion*num_col_ion]);

% CBM and VBM of p-type
doping_percentage = 0.2;
x_Ep = linspace(p_type_x_min, p_type_x_max, 1000);
y_Ecp = -p_type_y_min * ones(size(x_Ep));
y_Evp = -p_type_y_max * ones(size(x_Ep));
y_Efp = doping_percentage * y_Ecp + (1 - doping_percentage) * y_Evp;

% CBM and VBM of n-type
x_En = linspace(n_type_x_min, n_type_x_max, 1000);
y_Ecn = -n_type_y_min * ones(size(x_En));
y_Evn = -n_type_y_max * ones(size(x_En));
y_Efn = doping_percentage * y_Evn + (1 - doping_percentage) * y_Ecn;

while x_p(1) < x_n(1)
    fill(x_p, y_p, c_p, 'EdgeColor', 'g'); hold on;
    fill(x_n, y_n, c_n, 'EdgeColor', 'g'); hold on;
    
    x_p = x_p + step;
    x_n = x_n - step;
    p_type_x_min = p_type_x_min + step;
    p_type_x_max = p_type_x_max + step;
    n_type_x_min = n_type_x_min - step;
    n_type_x_max = n_type_x_max - step;
    
    x_ac = x_ac + step;
    scatter(x_ac, y_ac, ...
        10, 'x', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    
    x_dc = x_dc - step;
    scatter(x_dc, y_dc, ...
        10, '+', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    
    x_hc = x_hc + step + step * randn(1, num_row_carrier*num_col_carrier);
    x_hc(x_hc >= x_p(1)) = 2 * x_p(1) - x_hc(x_hc >= x_p(1));
    x_hc(x_hc <= x_p(3)) = 2 * x_p(3) - x_hc(x_hc <= x_p(3));
    y_hc = y_hc + step * randn(1, num_row_carrier*num_col_carrier);
    y_hc(y_hc >= y_p(2)) = 2 * y_p(2) - y_hc(y_hc >= y_p(2));
    y_hc(y_hc <= y_p(1)) = 2 * y_p(1) - y_hc(y_hc <= y_p(1));
    scatter(x_hc, y_hc, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    
    x_ec = x_ec - step - step * randn(1, num_row_carrier*num_col_carrier);
    x_ec(x_ec >= x_n(3)) = 2 * x_n(3) - x_ec(x_ec >= x_n(3));
    x_ec(x_ec <= x_n(1)) = 2 * x_n(1) - x_ec(x_ec <= x_n(1));
    y_ec = y_ec - step * randn(1, num_row_carrier*num_col_carrier);
    y_ec(y_ec >= y_n(2)) = 2 * y_n(2) - y_ec(y_ec >= y_n(2));
    y_ec(y_ec <= y_n(1)) = 2 * y_n(1) - y_ec(y_ec <= y_n(1));
    scatter(x_ec, y_ec, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    
    x_Ep = x_Ep + step;
    plot(x_Ep, y_Ecp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Evp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Efp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    x_En = x_En - step;
    plot(x_En, y_Ecn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Evn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Efn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    axis([-5 5 -4 3]);
    axis off; hold off;
    pause(pause_time);
end
pause(pause_interval); hold on;

%% Thermal Equilibrium
recombination_tolerance = 0.025;
equilibrium_percentage = 0.35;
max_frames = 450;
energy_step_exponential_coefficient = 0.988;
energy_step_initial = (y_Efn - y_Efp) / 2 * (1 - energy_step_exponential_coefficient);
energy_step = energy_step_initial;

p_barrier = equilibrium_percentage * p_type_x_min + (1 - equilibrium_percentage) * p_type_x_max;
n_barrier = equilibrium_percentage * n_type_x_max + (1 - equilibrium_percentage) * n_type_x_min;

diffusion_length_multiple = 1.5;
diffusion_length = (n_barrier - p_barrier) / 2 * diffusion_length_multiple;
x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - diffusion_length);
x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + diffusion_length);

x_h_static = x_hc(x_hc <  p_barrier); y_h_static = y_hc(x_hc <  p_barrier);
x_h_mobile = x_hc(x_hc >= p_barrier); y_h_mobile = y_hc(x_hc >= p_barrier);
x_e_static = x_ec(x_ec >  n_barrier); y_e_static = y_ec(x_ec >  n_barrier);
x_e_mobile = x_ec(x_ec <= n_barrier); y_e_mobile = y_ec(x_ec <= n_barrier);

v_h_max = 0.04; v_e_max = -0.04;
v_h = v_h_max * (x_h_mobile - p_barrier) / (p_type_x_max - p_barrier); % one-to-one mapping between velocity and position
v_e = v_e_max * (x_e_mobile - n_barrier) / (n_type_x_min - n_barrier);

for loop = 1 : max_frames
    
    
    x_h_static = x_h_static + step * randn(1, length(x_h_static));
    x_h_static(x_h_static >= p_barrier) = 2 * p_barrier - x_h_static(x_h_static >= p_barrier);
    x_h_static(x_h_static <= p_type_x_min) = 2 * p_type_x_min - x_h_static(x_h_static <= p_type_x_min);
    y_h_static = y_h_static + step * randn(1, length(y_h_static));
    y_h_static(y_h_static >= p_type_y_max) = 2 * p_type_y_max - y_h_static(y_h_static >= p_type_y_max);
    y_h_static(y_h_static <= p_type_y_min) = 2 * p_type_y_min - y_h_static(y_h_static <= p_type_y_min);
    
    
    x_e_static = x_e_static + step * randn(1, length(x_e_static));
    x_e_static(x_e_static <= n_barrier) = 2 * n_barrier - x_e_static(x_e_static <= n_barrier);
    x_e_static(x_e_static >= n_type_x_max) = 2 * n_type_x_max - x_e_static(x_e_static >= n_type_x_max);
    y_e_static = y_e_static + step * randn(1, length(y_e_static));
    y_e_static(y_e_static >= n_type_y_max) = 2 * n_type_y_max - y_e_static(y_e_static >= n_type_y_max);
    y_e_static(y_e_static <= n_type_y_min) = 2 * n_type_y_min - y_e_static(y_e_static <= n_type_y_min);
    
    
    x_h_mobile = x_h_mobile + v_h + step * randn(1, length(x_h_mobile));
    x_h_mobile(x_h_mobile < p_barrier) = 2 * p_barrier - x_h_mobile(x_h_mobile < p_barrier);
    y_h_mobile = y_h_mobile + step * randn(1, length(y_h_mobile));
    y_h_mobile(y_h_mobile >= p_type_y_max) = 2 * p_type_y_max - y_h_mobile(y_h_mobile >= p_type_y_max);
    y_h_mobile(y_h_mobile <= p_type_y_min) = 2 * p_type_y_min - y_h_mobile(y_h_mobile <= p_type_y_min);
    
    
    x_e_mobile = x_e_mobile + v_e + step * randn(1, length(x_e_mobile));
    x_e_mobile(x_e_mobile > n_barrier) = 2 * n_barrier - x_e_mobile(x_e_mobile > n_barrier);
    y_e_mobile = y_e_mobile + step * randn(1, length(y_e_mobile));
    y_e_mobile(y_e_mobile >= n_type_y_max) = 2 * n_type_y_max - y_e_mobile(y_e_mobile >= n_type_y_max);
    y_e_mobile(y_e_mobile <= n_type_y_min) = 2 * n_type_y_min - y_e_mobile(y_e_mobile <= n_type_y_min);
    
    
    outer_loop = 1;
    while (outer_loop < length(x_h_mobile))
        inner_loop = 1;
        while(outer_loop < length(x_h_mobile) && inner_loop < length(x_e_mobile))
            if ((x_h_mobile(outer_loop) - x_e_mobile(inner_loop))^2 + ...
                    (y_h_mobile(outer_loop) - y_e_mobile(inner_loop))^2 < recombination_tolerance^2)
                x_h_mobile(outer_loop) = [];
                y_h_mobile(outer_loop) = [];
                v_h(outer_loop) = [];
                x_e_mobile(inner_loop) = [];
                y_e_mobile(inner_loop) = [];
                v_e(inner_loop) = [];
            end
            inner_loop = inner_loop + 1;
        end
        outer_loop = outer_loop + 1;
    end
    
    y_Ecp = y_Ecp + energy_step; y_Evp = y_Evp + energy_step; y_Efp = y_Efp + energy_step;
    y_Ecn = y_Ecn - energy_step; y_Evn = y_Evn - energy_step; y_Efn = y_Efn - energy_step;
    energy_step = energy_step * energy_step_exponential_coefficient;
    
    y_Ecp(x_Ep > p_barrier) = max(y_Ecp) - (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_Ep(x_Ep > p_barrier) - p_barrier) / (p_type_x_max - p_barrier)).^2;
    y_Evp = y_Ecp - (max(y_Ecp) - max(y_Evp));
    y_Ecn(x_En < n_barrier) = min(y_Ecn) + (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_En(x_En < n_barrier) - n_barrier) / (n_type_x_min - n_barrier)).^2;
    y_Evn = y_Ecn - (min(y_Ecn) - min(y_Evn));
    
    y_quasi_fermi_p_region = linspace(min(y_Efp), max(y_Efn), length(x_quasi_fermi_p_region));
    y_quasi_fermi_n_region = linspace(min(y_Efp), max(y_Efn), length(x_quasi_fermi_n_region));
    
    
    
    outer_loop = 1;
    while (outer_loop < length(x_h_mobile))
        inner_loop = 1;
        while(outer_loop < length(x_h_mobile) && inner_loop < length(x_e_static))
            if ((x_h_mobile(outer_loop) - x_e_static(inner_loop))^2 + ...
                    (y_h_mobile(outer_loop) - y_e_static(inner_loop))^2 < recombination_tolerance^2)
                x_h_mobile(outer_loop) = [];
                y_h_mobile(outer_loop) = [];
                v_h(outer_loop) = [];
            end
            inner_loop = inner_loop + 1;
        end
        outer_loop = outer_loop + 1;
    end
    
    
    
    outer_loop = 1;
    while (outer_loop < length(x_h_static))
        inner_loop = 1;
        while(outer_loop < length(x_h_static) && inner_loop < length(x_e_mobile))
            if ((x_h_static(outer_loop) - x_e_mobile(inner_loop))^2 + ...
                    (y_h_static(outer_loop) - y_e_mobile(inner_loop))^2 < recombination_tolerance^2)
                x_e_mobile(inner_loop) = [];
                y_e_mobile(inner_loop) = [];
                v_e(inner_loop) = [];
            end
            inner_loop = inner_loop + 1;
        end
        outer_loop = outer_loop + 1;
    end
    
    fill(x_p, y_p, c_p, 'EdgeColor', 'g'); hold on;
    fill(x_n, y_n, c_n, 'EdgeColor', 'g'); hold on;
    scatter(x_ac, y_ac, ...
        10, 'x', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_dc, y_dc, ...
        10, '+', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_h_static, y_h_static, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    scatter(x_e_static, y_e_static, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    scatter(x_h_mobile, y_h_mobile, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    scatter(x_e_mobile, y_e_mobile, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    plot(x_Ep, y_Ecp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Evp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Efp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_quasi_fermi_p_region, y_quasi_fermi_p_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_quasi_fermi_n_region, y_quasi_fermi_n_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_En, y_Ecn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    
    
    plot(x_En, y_Evn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Efn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    y_h_mobile(x_h_mobile > n_type_x_max) = []; v_h(x_h_mobile > n_type_x_max) = []; x_h_mobile(x_h_mobile > n_type_x_max) = [];
    y_e_mobile(x_e_mobile < p_type_x_min) = []; v_e(x_e_mobile < p_type_x_min) = []; x_e_mobile(x_e_mobile < p_type_x_min) = [];
    
    axis([-5 5 -4 3]);
    axis off; hold off;
    pause(pause_time);
end
pause(pause_interval);

%% Forward Bias

x_e = [x_e_static x_e_mobile]; y_e = [y_e_static y_e_mobile];
x_h = [x_h_static x_h_mobile]; y_h = [y_h_static y_h_mobile];
max_frames = 500;
velocity_depletion_region = 3 * step;
forward_percentage = 0.2;
p_barrier_temp = p_barrier;
p_barrier = forward_percentage * p_type_x_min + (1 - forward_percentage) * p_type_x_max;
n_barrier_temp = n_barrier;
n_barrier = forward_percentage * n_type_x_max + (1 - forward_percentage) * n_type_x_min;
recombination_tolerance = 0.02;
reflection_threshold = 0.85;

energy_step_exponential_coefficient = 0.99;
energy_step_initial = (max(y_Ecp) - min(y_Ecn)) / 2 * (1 - energy_step_exponential_coefficient) * 0.5;
energy_step = energy_step_initial;
diffusion_length_multiple = 0.8;
junction_width = n_barrier - p_barrier;
diffusion_length = junction_width / 2 * diffusion_length_multiple;
x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - junction_width / 2 - diffusion_length);
x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + junction_width / 2 + diffusion_length);
y_quasi_fermi_p_region = zeros(size(x_quasi_fermi_p_region));
y_quasi_fermi_n_region = zeros(size(x_quasi_fermi_n_region));

for loop = 1 : max_frames
    x_h(x_h <  p_barrier) = x_h(x_h <  p_barrier) + step + step * randn(1, length(x_h(x_h < p_barrier)));
    x_h(x_h >= p_barrier + step) = x_h(x_h >= p_barrier + step) + velocity_depletion_region + step * randn(1, length(x_h(x_h >= p_barrier + step)));
    temp_index = (x_h >= p_barrier & x_h < p_barrier + step & rand(1, length(x_h >= p_barrier & x_h < p_barrier + step)) > reflection_threshold);
    x_h(temp_index) = 2 * p_barrier - x_h(temp_index);
    x_h(x_h <= p_type_x_min) = 2 * p_type_x_min - x_h(x_h <= p_type_x_min);
    y_h = y_h + step * randn(1, length(y_h));
    y_h(y_h <= p_type_y_min) = 2 * p_type_y_min - y_h(y_h <= p_type_y_min);
    y_h(y_h >= p_type_y_max) = 2 * p_type_y_max - y_h(y_h >= p_type_y_max);
    
    x_e(x_e >  n_barrier) = x_e(x_e >  n_barrier) - step + step * randn(1, length(x_e(x_e > n_barrier)));
    x_e(x_e <= n_barrier - step) = x_e(x_e <= n_barrier - step) - velocity_depletion_region + step * randn(1, length(x_e(x_e <= n_barrier - step)));
    temp_index = (x_e <= n_barrier & x_e > n_barrier - step & rand(1, length(x_e <= n_barrier & x_e > n_barrier - step)) > reflection_threshold);
    x_e(temp_index) = 2 * n_barrier - x_e(temp_index);
    x_e(x_e >= n_type_x_max) = 2 * n_type_x_max - x_e(x_e >= n_type_x_max);
    y_e = y_e + step * randn(1, length(y_e));
    y_e(y_e <= n_type_y_min) = 2 * n_type_y_min - y_e(y_e <= n_type_y_min);
    y_e(y_e >= n_type_y_max) = 2 * n_type_y_max - y_e(y_e >= n_type_y_max);
    
    outer_loop = 1;
    while (outer_loop < length(x_h))
        inner_loop = 1;
        while(outer_loop < length(x_h) && inner_loop < length(x_e))
            if ((x_h(outer_loop) - x_e(inner_loop))^2 + ...
                    (y_h(outer_loop) - y_e(inner_loop))^2 < recombination_tolerance^2)
                if ((x_h(outer_loop) > n_barrier && x_e(inner_loop) > n_barrier) ||...
                        (x_h(outer_loop) < p_barrier && x_e(inner_loop) < p_barrier))
                    x_h(outer_loop) = [];
                    y_h(outer_loop) = [];
                    x_e(inner_loop) = [];
                    y_e(inner_loop) = [];
                end
            end
            inner_loop = inner_loop + 1;
        end
        outer_loop = outer_loop + 1;
    end
    
    y_h_mobile(x_h_mobile > n_type_x_max) = []; x_h_mobile(x_h_mobile > n_type_x_max) = [];
    y_e_mobile(x_e_mobile < p_type_x_min) = []; x_e_mobile(x_e_mobile < p_type_x_min) = [];
    
    if (mod(loop, ceil((p_type_x_max - p_type_x_min) / num_col_carrier / step)) == 0)
        x_h = [x_h, p_type_x_min * ones(1, num_row_carrier)];
        y_h = [y_h, linspace(p_type_y_min, p_type_y_max, num_row_carrier)];
        x_e = [x_e, n_type_x_max * ones(1, num_row_carrier)];
        y_e = [y_e, linspace(n_type_y_min, n_type_y_max, num_row_carrier)];
    end
    
    fill(x_p, y_p, c_p, 'EdgeColor', 'g'); hold on;
    fill(x_n, y_n, c_n, 'EdgeColor', 'g'); hold on;
    scatter(x_ac, y_ac, ...
        10, 'x', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_dc, y_dc, ...
        10, '+', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_h, y_h, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    scatter(x_e, y_e, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    
    if (p_barrier_temp < p_barrier && n_barrier_temp > n_barrier)
        p_barrier_temp = p_barrier_temp + step;
        n_barrier_temp = n_barrier_temp - step;
        
        y_Ecp(x_Ep > p_barrier_temp) = max(y_Ecp) - (max(y_Ecp) - min(y_Ecn)) / 2 *...
            ((x_Ep(x_Ep > p_barrier_temp) - p_barrier_temp) / (p_type_x_max - p_barrier_temp)).^2;
        y_Evp = y_Ecp - (max(y_Ecp) - max(y_Evp));
        y_Ecn(x_En < n_barrier_temp) = min(y_Ecn) + (max(y_Ecp) - min(y_Ecn)) / 2 *...
            ((x_En(x_En < n_barrier_temp) - n_barrier_temp) / (n_type_x_min - n_barrier_temp)).^2;
        y_Evn = y_Ecn - (min(y_Ecn) - min(y_Evn));
    else
        y_Ecp = y_Ecp - energy_step; y_Evp = y_Evp - energy_step; y_Efp = y_Efp - energy_step;
        y_Ecn = y_Ecn + energy_step; y_Evn = y_Evn + energy_step; y_Efn = y_Efn + energy_step;
        energy_step = energy_step * energy_step_exponential_coefficient;
        
        y_Ecp(x_Ep > p_barrier) = max(y_Ecp) - (max(y_Ecp) - min(y_Ecn)) / 2 *...
            ((x_Ep(x_Ep > p_barrier) - p_barrier) / (p_type_x_max - p_barrier)).^2;
        y_Evp = y_Ecp - (max(y_Ecp) - max(y_Evp));
        y_Ecn(x_En < n_barrier) = min(y_Ecn) + (max(y_Ecp) - min(y_Ecn)) / 2 *...
            ((x_En(x_En < n_barrier) - n_barrier) / (n_type_x_min - n_barrier)).^2;
        y_Evn = y_Ecn - (min(y_Ecn) - min(y_Evn));
        
        y_quasi_fermi_p_region(x_quasi_fermi_p_region >  p_barrier) = max(y_Efn);
        y_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier) = max(y_Efn) + (x_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier) - p_barrier) / diffusion_length * (max(y_Efn) - min(y_Efp));
        y_quasi_fermi_n_region(x_quasi_fermi_n_region < n_barrier) = min(y_Efp);
        y_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier) = min(y_Efp) + (x_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier) - n_barrier) / diffusion_length * (max(y_Efn) - min(y_Efp));
        plot(x_quasi_fermi_p_region, y_quasi_fermi_p_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
        plot(x_quasi_fermi_n_region, y_quasi_fermi_n_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
        
    end
    
    plot(x_Ep, y_Ecp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Evp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Efp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_En, y_Ecn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Evn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Efn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    axis([-5 5 -4 3]);
    axis off; hold off;
    
    pause(pause_time);
end
pause(pause_interval);

%% Reverse Bias

max_frames = 750;
velocity_depletion_region = 4 * step;
reverse_percentage = 0.4;
p_barrier_temp = p_barrier;
p_barrier = reverse_percentage * p_type_x_min + (1 - reverse_percentage) * p_type_x_max;
n_barrier_temp = n_barrier;
n_barrier = reverse_percentage * n_type_x_max + (1 - reverse_percentage) * n_type_x_min;
reflection_threshold = 0.25;

energy_step_exponential_coefficient = 0.985;
energy_step_initial = (max(y_Ecp) - min(y_Ecn)) / 2 * (1 - energy_step_exponential_coefficient) * 2;
energy_step = energy_step_initial;
diffusion_length_multiple = 0.8;
junction_width = n_barrier_temp - p_barrier_temp;
diffusion_length = junction_width / 2 * diffusion_length_multiple;
x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - junction_width / 2 - diffusion_length);
x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + junction_width / 2 + diffusion_length);
y_quasi_fermi_p_region = zeros(size(x_quasi_fermi_p_region));
y_quasi_fermi_n_region = zeros(size(x_quasi_fermi_n_region));
generation_rate = 1;
generation_cycle = 4;

minority_new_percentage = 0.05;
x_h = [x_h, n_barrier_temp + (n_type_x_max - n_barrier_temp) * rand(1, round(minority_new_percentage * length(x_h)))];
y_h = [y_h, n_type_y_min + (n_type_y_max - n_type_y_min) * rand(1, round(minority_new_percentage * length(y_h)))];
x_e = [x_e, p_type_x_min + (p_barrier_temp - p_type_x_min) * rand(1, round(minority_new_percentage * length(x_e)))];
y_e = [y_e, p_type_y_min + (p_type_y_max - p_type_y_min) * rand(1, round(minority_new_percentage * length(y_e)))];

for loop = 1 : max_frames
    static_index_temp = (x_h < p_barrier_temp) & (rand(1, length(x_h)) < reflection_threshold);
    mobile_index_temp = (x_h < p_barrier_temp) & (~static_index_temp);
    x_h(static_index_temp) = x_h(static_index_temp) + step * randn(1, sum(static_index_temp));
    x_h(mobile_index_temp) = x_h(mobile_index_temp) - step + step * randn(1, sum(mobile_index_temp));
    x_h(x_h > n_barrier_temp) = x_h(x_h > n_barrier_temp) - step + step * randn(1, length(x_h(x_h > n_barrier_temp)));
    x_h(x_h >= p_barrier_temp & x_h <= n_barrier_temp) = x_h(x_h >= p_barrier_temp & x_h <= n_barrier_temp) - velocity_depletion_region + step * randn(1, length(x_h(x_h >= p_barrier_temp & x_h <= n_barrier_temp)));
    disappear_index_temp = (x_h <= p_type_x_min);
    x_h(disappear_index_temp) = [];
    y_h = y_h + step * randn(1, length(y_h));
    y_h(y_h <= p_type_y_min) = 2 * p_type_y_min - y_h(y_h <= p_type_y_min);
    y_h(y_h >= p_type_y_max) = 2 * p_type_y_max - y_h(y_h >= p_type_y_max);
    y_h(disappear_index_temp) = [];
    
    static_index_temp = (x_e > n_barrier_temp) & (rand(1, length(x_e)) < reflection_threshold);
    mobile_index_temp = (x_e > n_barrier_temp) & (~static_index_temp);
    x_e(static_index_temp) = x_e(static_index_temp) + step * randn(1, sum(static_index_temp));
    x_e(mobile_index_temp) = x_e(mobile_index_temp) + step + step * randn(1, sum(mobile_index_temp));
    x_e(x_e < p_barrier_temp) = x_e(x_e < p_barrier_temp) + step + step * randn(1, length(x_e(x_e < p_barrier_temp)));
    x_e(x_e <= n_barrier_temp & x_e >= p_barrier_temp) = x_e(x_e <= n_barrier_temp & x_e >= p_barrier_temp) + velocity_depletion_region + step * randn(1, length(x_e(x_e <= n_barrier_temp & x_e >= p_barrier_temp)));
    disappear_index_temp = (x_e >= n_type_x_max);
    x_e(disappear_index_temp) = [];
    y_e = y_e + step * randn(1, length(y_e));
    y_e(y_e <= n_type_y_min) = 2 * n_type_y_min - y_e(y_e <= n_type_y_min);
    y_e(y_e >= n_type_y_max) = 2 * n_type_y_max - y_e(y_e >= n_type_y_max);
    y_e(disappear_index_temp) = [];
    
    fill(x_p, y_p, c_p, 'EdgeColor', 'g'); hold on;
    fill(x_n, y_n, c_n, 'EdgeColor', 'g'); hold on;
    scatter(x_ac, y_ac, ...
        10, 'x', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_dc, y_dc, ...
        10, '+', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_h, y_h, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    scatter(x_e, y_e, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    
    if (p_barrier_temp > p_barrier && n_barrier_temp < n_barrier)
        p_barrier_temp = p_barrier_temp - step;
        n_barrier_temp = n_barrier_temp + step;
    else
        reflection_threshold = 0.8;
        y_Ecp = y_Ecp + energy_step; y_Evp = y_Evp + energy_step; y_Efp = y_Efp + energy_step;
        y_Ecn = y_Ecn - energy_step; y_Evn = y_Evn - energy_step; y_Efn = y_Efn - energy_step;
        energy_step = energy_step * energy_step_exponential_coefficient; 
        if (mod(loop, generation_cycle) == 0)
            
            x_h = [x_h, n_barrier_temp + diffusion_length * rand(1, generation_rate), p_barrier_temp - diffusion_length * rand(1, generation_rate)];
            y_h = [y_h, n_type_y_min + (n_type_y_max - n_type_y_min) * rand(1, generation_rate), p_type_y_min + (p_type_y_max - p_type_y_min) * rand(1, generation_rate)];
            x_e = [x_e, p_barrier_temp - diffusion_length * rand(1, generation_rate), n_barrier_temp + diffusion_length * rand(1, generation_rate)];
            y_e = [y_e, p_type_y_min + (p_type_y_max - p_type_y_min) * rand(1, generation_rate), n_type_y_min + (n_type_y_max - n_type_y_min) * rand(1, generation_rate)];
        end
    end
    
    junction_width = n_barrier_temp - p_barrier_temp;
    diffusion_length = junction_width / 2 * diffusion_length_multiple;
    x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - junction_width / 2 - diffusion_length);
    x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + junction_width / 2 + diffusion_length);
    y_quasi_fermi_p_region = zeros(size(x_quasi_fermi_p_region));
    y_quasi_fermi_n_region = zeros(size(x_quasi_fermi_n_region));
    
    y_Ecp(x_Ep > p_barrier_temp) = max(y_Ecp) - (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_Ep(x_Ep > p_barrier_temp) - p_barrier_temp) / (p_type_x_max - p_barrier_temp)).^2;
    y_Evp = y_Ecp - (max(y_Ecp) - max(y_Evp));
    y_Ecn(x_En < n_barrier_temp) = min(y_Ecn) + (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_En(x_En < n_barrier_temp) - n_barrier_temp) / (n_type_x_min - n_barrier_temp)).^2;
    y_Evn = y_Ecn - (min(y_Ecn) - min(y_Evn));
    
    y_quasi_fermi_p_region(x_quasi_fermi_p_region >  p_barrier_temp) = max(y_Efn);
    y_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier_temp) = max(y_Efn) + (x_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier_temp) - p_barrier_temp) / diffusion_length * (max(y_Efn) - min(y_Efp));
    y_quasi_fermi_n_region(x_quasi_fermi_n_region <  n_barrier_temp) = min(y_Efp);
    y_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier_temp) = min(y_Efp) + (x_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier_temp) - n_barrier_temp) / diffusion_length * (max(y_Efn) - min(y_Efp));
    plot(x_quasi_fermi_p_region, y_quasi_fermi_p_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_quasi_fermi_n_region, y_quasi_fermi_n_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    
    plot(x_Ep, y_Ecp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Evp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Efp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_En, y_Ecn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Evn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Efn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    axis([-5 5 -4 3]);
    axis off; hold off;
    
    pause(pause_time);
end
pause(pause_interval);

%% Light Excitation
max_frames = 500;
velocity_depletion_region = 3 * step;
equilibrium_percentage = 0.3;
p_barrier_temp = p_barrier;
p_barrier = equilibrium_percentage * p_type_x_min + (1 - equilibrium_percentage) * p_type_x_max;
n_barrier_temp = n_barrier;
n_barrier = equilibrium_percentage * n_type_x_max + (1 - equilibrium_percentage) * n_type_x_min;

energy_step_exponential_coefficient = 0.99;
energy_step_initial = (max(y_Ecp) - min(y_Ecn)) / 2 * (1 - energy_step_exponential_coefficient) * 0.5;
energy_step = energy_step_initial;
diffusion_length_multiple = 0.8;
junction_width = n_barrier_temp - p_barrier_temp;
diffusion_length = junction_width / 2 * diffusion_length_multiple;
x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - junction_width / 2 - diffusion_length);
x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + junction_width / 2 + diffusion_length);
y_quasi_fermi_p_region = zeros(size(x_quasi_fermi_p_region));
y_quasi_fermi_n_region = zeros(size(x_quasi_fermi_n_region));
reflection_threshold = 0.7;
recombination_tolerance = 0.02;
generation_rate = 1;
generation_cycle = 3;

for loop = 1 : max_frames
    if (p_barrier_temp < p_barrier && n_barrier_temp > n_barrier)
        p_barrier_temp = p_barrier_temp + step;
        n_barrier_temp = n_barrier_temp - step;
        if (mod(loop, ceil((p_type_x_max - p_type_x_min) / num_col_carrier / step)) == 0)
            x_h = [x_h, p_type_x_min * ones(1, num_row_carrier)];
            y_h = [y_h, linspace(p_type_y_min, p_type_y_max, num_row_carrier)];
            x_e = [x_e, n_type_x_max * ones(1, num_row_carrier)];
            y_e = [y_e, linspace(n_type_y_min, n_type_y_max, num_row_carrier)];
        end
        x_h = x_h + step + step * randn(1, length(x_h));
        y_h = y_h + step * randn(1, length(y_h));
        y_h(y_h <= p_type_y_min) = 2 * p_type_y_min - y_h(y_h <= p_type_y_min);
        y_h(y_h >= p_type_y_max) = 2 * p_type_y_max - y_h(y_h >= p_type_y_max);
        x_e = x_e - step + step * randn(1, length(x_e));
        y_e = y_e + step * randn(1, length(y_e));
        y_e(y_e <= n_type_y_min) = 2 * n_type_y_min - y_e(y_e <= n_type_y_min);
        y_e(y_e >= n_type_y_max) = 2 * n_type_y_max - y_e(y_e >= n_type_y_max);
    else
        y_Ecp = y_Ecp - energy_step; y_Evp = y_Evp - energy_step; y_Efp = y_Efp - energy_step;
        y_Ecn = y_Ecn + energy_step; y_Evn = y_Evn + energy_step; y_Efn = y_Efn + energy_step;
        energy_step = energy_step * energy_step_exponential_coefficient; 
        
        x_h = x_h + step * randn(1, length(x_h));
        y_h = y_h + step * randn(1, length(y_h));
        x_h(x_h > p_barrier_temp & x_h < n_barrier_temp) = x_h(x_h > p_barrier_temp & x_h < n_barrier_temp) - velocity_depletion_region;
        temp_select_index = (x_h <= p_barrier_temp & rand(1, length(x_h)) > reflection_threshold);
        x_h(temp_select_index) = x_h(temp_select_index) - step;
        y_h(y_h <= p_type_y_min) = 2 * p_type_y_min - y_h(y_h <= p_type_y_min);
        y_h(y_h >= p_type_y_max) = 2 * p_type_y_max - y_h(y_h >= p_type_y_max);
        disappear_index_temp = (x_h <= p_type_x_min);
        x_h(disappear_index_temp) = [];
        y_h(disappear_index_temp) = [];
        
        x_e = x_e + step * randn(1, length(x_e));
        y_e = y_e + step * randn(1, length(y_e));
        x_e(x_e > p_barrier_temp & x_e < n_barrier_temp) = x_e(x_e > p_barrier_temp & x_e < n_barrier_temp) + velocity_depletion_region;
        temp_select_index = (x_e >= n_barrier_temp & rand(1, length(x_e)) > reflection_threshold);
        x_e(temp_select_index) = x_e(temp_select_index) + step;
        y_e(y_e <= n_type_y_min) = 2 * n_type_y_min - y_e(y_e <= n_type_y_min);
        y_e(y_e >= n_type_y_max) = 2 * n_type_y_max - y_e(y_e >= n_type_y_max);
        disappear_index_temp = (x_e >= n_type_x_max);
        x_e(disappear_index_temp) = [];
        y_e(disappear_index_temp) = [];
        
        outer_loop = 1;
        while (outer_loop < length(x_h))
            inner_loop = 1;
            while(outer_loop < length(x_h) && inner_loop < length(x_e))
                if ((x_h(outer_loop) - x_e(inner_loop))^2 + ...
                        (y_h(outer_loop) - y_e(inner_loop))^2 < recombination_tolerance^2)
                    if ((x_h(outer_loop) > n_barrier_temp && x_e(inner_loop) > n_barrier_temp) ||...
                            (x_h(outer_loop) < p_barrier_temp && x_e(inner_loop) < p_barrier_temp))
                        x_h(outer_loop) = [];
                        y_h(outer_loop) = [];
                        x_e(inner_loop) = [];
                        y_e(inner_loop) = [];
                    end
                end
                inner_loop = inner_loop + 1;
            end
            outer_loop = outer_loop + 1;
        end
        
        if (mod(loop, generation_cycle) == 0)
            x_generation = p_barrier_temp + (n_barrier_temp - p_barrier_temp) * rand(1, generation_rate);
            y_generation = p_type_y_min + (p_type_y_max - p_type_y_min) * rand(1, generation_rate);
            x_h = [x_h, x_generation];
            y_h = [y_h, y_generation];
            x_e = [x_e, x_generation];
            y_e = [y_e, y_generation];
        end
    end   
    
    junction_width = n_barrier_temp - p_barrier_temp;
    diffusion_length = junction_width / 2 * diffusion_length_multiple;
    x_quasi_fermi_p_region = x_Ep(x_Ep > p_type_x_max - junction_width / 2 - diffusion_length);
    x_quasi_fermi_n_region = x_En(x_En < n_type_x_min + junction_width / 2 + diffusion_length);
    y_quasi_fermi_p_region = zeros(size(x_quasi_fermi_p_region));
    y_quasi_fermi_n_region = zeros(size(x_quasi_fermi_n_region));
    
    y_Ecp(x_Ep > p_barrier_temp) = max(y_Ecp) - (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_Ep(x_Ep > p_barrier_temp) - p_barrier_temp) / (p_type_x_max - p_barrier_temp)).^2;
    y_Evp = y_Ecp - (max(y_Ecp) - max(y_Evp));
    y_Ecn(x_En < n_barrier_temp) = min(y_Ecn) + (max(y_Ecp) - min(y_Ecn)) / 2 *...
        ((x_En(x_En < n_barrier_temp) - n_barrier_temp) / (n_type_x_min - n_barrier_temp)).^2;
    y_Evn = y_Ecn - (min(y_Ecn) - min(y_Evn));
    
    y_quasi_fermi_p_region(x_quasi_fermi_p_region >  p_barrier_temp) = max(y_Efn);
    y_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier_temp) = max(y_Efn) + (x_quasi_fermi_p_region(x_quasi_fermi_p_region <= p_barrier_temp) - p_barrier_temp) / diffusion_length * (max(y_Efn) - min(y_Efp));
    y_quasi_fermi_n_region(x_quasi_fermi_n_region <  n_barrier_temp) = min(y_Efp);
    y_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier_temp) = min(y_Efp) + (x_quasi_fermi_n_region(x_quasi_fermi_n_region >= n_barrier_temp) - n_barrier_temp) / diffusion_length * (max(y_Efn) - min(y_Efp));
    
    fill(x_p, y_p, c_p, 'EdgeColor', 'g'); hold on;
    fill(x_n, y_n, c_n, 'EdgeColor', 'g'); hold on;
    scatter(x_ac, y_ac, ...
        10, 'x', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_dc, y_dc, ...
        10, '+', 'LineWidth', 7, 'MarkerEdgeColor', [0.2 0.2 0.2]); hold on;
    scatter(x_h, y_h, ...
        7, 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]); hold on;
    scatter(x_e, y_e, ...
        7, 'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', [0 0 1]); hold on;
    

    plot(x_quasi_fermi_p_region, y_quasi_fermi_p_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_quasi_fermi_n_region, y_quasi_fermi_n_region, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    
    plot(x_Ep, y_Ecp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Evp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_Ep, y_Efp, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    plot(x_En, y_Ecn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Evn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-'); hold on;
    plot(x_En, y_Efn, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', ':'); hold on;
    
    axis([-5 5 -4 3]);
    axis off; hold off;
    
    pause(pause_time);
end
pause(pause_interval);