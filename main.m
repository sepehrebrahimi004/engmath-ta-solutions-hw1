%MAIN  Engineering Mathematics - TA solutions, MATLAB homework 1: Fourier series.
%      Every numerically computed coefficient is checked against its closed form.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

N  = 20001;                                         % fine grid: accurate trapz up to n = 500
tp = linspace(-pi, pi, N);
trig = @(t, f, n) deal(trapz(t, f)/(2*pi), ...                     % a0, a_n, b_n over a 2*pi period
                       arrayfun(@(m) trapz(t, f.*cos(m*t))/pi, n), ...
                       arrayfun(@(m) trapz(t, f.*sin(m*t))/pi, n));

%% Q1  Trigonometric series of sinh t, sign t and the triangle pi - |t| on (-pi, pi)
f = sinh(tp);  g = sign(tp);  h = pi - abs(tp);
n = 1:30;
[a0f, af, bf] = trig(tp, f, n);  [a0g, ag, bg] = trig(tp, g, n);  [a0h, ah, bh] = trig(tp, h, n);
bf_th = 2*(-1).^(n+1) .* n * sinh(pi) ./ (pi*(n.^2 + 1));   % odd: cosine terms vanish
bg_th = 4./(pi*n) .* mod(n, 2);                             % square wave: odd harmonics
ah_th = 4./(pi*n.^2) .* mod(n, 2);                          % triangle: a0 = pi/2, odd harmonics
fprintf('Q1  max error vs closed form:  sinh %.1e   sign %.1e   triangle %.1e\n', ...
        max(abs([af, bf - bf_th])), max(abs([ag, bg - bg_th])), max(abs([a0h - pi/2, ah - ah_th, bh])));
fprintf('    decay: |b_n(sinh)| ~ 1/n,  |b_n(sign)| ~ 1/n (jumps),  |a_n(triangle)| ~ 1/n^2 (continuous)\n');
fig = figure('Visible', 'off', 'Position', [100 100 1150 600]);
names = {'sinh t', 'sign t', '\pi - |t|'};  sigs = {f, g, h};
coef  = {bf, bg, ah};  theo = {bf_th, bg_th, ah_th};  lab = {'b_n', 'b_n', 'a_n'};
for i = 1:3
    subplot(2, 3, i); plot(tp, sigs{i}, 'LineWidth', 1.3); grid on; xlim([-pi pi]); title(names{i});
    subplot(2, 3, i + 3); stem(n, coef{i}, 'filled', 'MarkerSize', 3); hold on;
    plot(n, theo{i}, 'rx'); grid on; xlabel('n'); title([lab{i} ' (stems: numerical, x: closed form)']);
end
save_figure(fig, 'q1_series_coefficients');

%% Q2  Complex coefficients of cos 2pi t + 5 sin 4pi t + cos 4pi t + 10 cos 10pi t  (T = 1)
t2 = linspace(-0.5, 0.5, N);
f2 = cos(2*pi*t2) + 5*sin(4*pi*t2) + cos(4*pi*t2) + 10*cos(10*pi*t2);
k = -10:10;
c = arrayfun(@(m) trapz(t2, f2.*exp(-1j*2*pi*m*t2)), k);
c_th = zeros(size(k));
c_th(abs(k) == 1) = 0.5;  c_th(k == 2) = 0.5 - 2.5j;  c_th(k == -2) = 0.5 + 2.5j;  c_th(abs(k) == 5) = 5;
fprintf('Q2  max |c - c_theory| = %.1e;  c_2 = %.2f%+.2fj\n', max(abs(c - c_th)), real(c(k == 2)), imag(c(k == 2)));
c(abs(c) < 1e-9) = 0;                                       % phase of a zero coefficient is noise
fig = figure('Visible', 'off', 'Position', [100 100 800 450]);
subplot(2, 1, 1); stem(k, abs(c), 'filled'); grid on; title('|c_n|');
subplot(2, 1, 2); stem(k, rad2deg(angle(c)), 'filled'); grid on; title('\anglec_n (deg)'); xlabel('n');
save_figure(fig, 'q2_complex_coefficients');

%% Q3a  Staircase -2, -1, 1, 2 on (0, 2pi): partial sums
t3 = linspace(0, 2*pi, N);
f3 = -2*(t3 < pi/2) - (t3 >= pi/2 & t3 < pi) + (t3 >= pi & t3 < 3*pi/2) + 2*(t3 >= 3*pi/2);
[a0, a, b] = trig(t3, f3, 1:100);
fig = figure('Visible', 'off', 'Position', [100 100 900 520]);
Ns = [1 5 20 100];
for i = 1:numel(Ns)
    m = 1:Ns(i);
    fr = a0 + a(m)*cos(m.'*t3) + b(m)*sin(m.'*t3);
    subplot(2, 2, i); plot(t3, f3, 'b', t3, fr, 'r', 'LineWidth', 1.1); grid on; xlim([0 2*pi]);
    title(sprintf('n = %d', Ns(i)));
end
fr = a0 + a*cos((1:100).'*t3) + b*sin((1:100).'*t3);
fprintf('Q3a overshoot at the 4-unit jump (t = 0 / 2pi) with 100 terms: %.1f%% (Gibbs: 8.9%%)\n', ...
        100*(max(fr) - 2)/4);
save_figure(fig, 'q3a_staircase_partial_sums');

%% Q3b  g(t) = cos(4t) e^{jt} = (e^{j5t} + e^{-j3t})/2, a complex curve
g3 = cos(4*t3) .* exp(1j*t3);
[a0, a, b] = trig(t3, g3, 1:5);                             % trig formulas hold for complex g too
fig = figure('Visible', 'off', 'Position', [100 100 1100 330]);
for K = 3:5
    m = 1:K;
    gr = a0 + a(m)*cos(m.'*t3) + b(m)*sin(m.'*t3);
    subplot(1, 3, K - 2); plot(real(g3), imag(g3), 'b', 'LineWidth', 2); hold on;
    plot(real(gr), imag(gr), 'r', 'LineWidth', 1.1); axis equal; grid on;
    title(sprintf('n = 1..%d  (error %.1e)', K, max(abs(gr - g3))));
end
legend('g(t)', 'partial sum', 'Location', 'southoutside');
fprintf('Q3b only harmonics 3 and 5 are present: a_3 = %.2f%+.2fj, a_5 = %.2f, b_5 = %.2fj\n', ...
        real(a(3)), imag(a(3)), real(a(5)), imag(b(5)));
save_figure(fig, 'q3b_complex_curve');

%% Q4  c'_n = c_{-n} - c_n of a square wave gives g(t) = f(-t) - f(t)
f4  = double(t3 < pi) - double(t3 >= pi);                  % +1 on (0, pi), -1 on (pi, 2pi)
K   = 500;  kk = -K:K;
c1  = arrayfun(@(m) trapz(t3, f4.*exp(-1j*m*t3)), kk)/(2*pi);
c2  = fliplr(c1) - c1;
tw  = linspace(-2*pi, 2*pi, 8001);
g4  = real(c2 * exp(1j*kk.'*tw));
fsq = @(x) 1 - 2*(mod(x, 2*pi) >= pi);
fig = figure('Visible', 'off', 'Position', [100 100 900 380]);
plot(tw, fsq(tw), 'b', tw, g4, 'r', 'LineWidth', 1.2); grid on;
legend('f(t)', 'g(t) from c''_n = c_{-n} - c_n', 'Location', 'southeast'); xlabel('t');
title('Reversing the coefficient index reverses time: g(t) = f(-t) - f(t) = -2 f(t)');
save_figure(fig, 'q4_coefficient_reversal');
fprintf('Q4  rms(g + 2f) away from the jumps = %.3f\n', ...
        sqrt(mean((g4(abs(sin(tw)) > 0.1) + 2*fsq(tw(abs(sin(tw)) > 0.1))).^2)));
