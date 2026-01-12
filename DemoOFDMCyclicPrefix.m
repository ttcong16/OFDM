clearvars; close all;

%% Parameters
n_subcarriers = 128;
n_symbols = 7;
cp_length = 12;
scs = 15e3;
fs = n_subcarriers * scs; %#ok<NASGU>
qam_order = 4;

%% OFDM Grid Generation & Modulation
qam_symbols = randi(qam_order, n_subcarriers, n_symbols) - 1;


% ofdm_grid = qammod(qam_symbols, qam_order, 'gray', 'UnitAveragePower', true);
ofdm_grid = qammod(qam_symbols, qam_order, 0, 'gray');
ref = qammod((0:qam_order-1).', qam_order, 0, 'gray');
scale = sqrt(mean(abs(ref).^2));
ofdm_grid = ofdm_grid ./ scale;
% ==========================================================================


% ofdm_signal = ofdmmod(ofdm_grid, n_subcarriers, cp_length);
ofdm_time = ifft(ofdm_grid, n_subcarriers, 1);
ofdm_signal = [ofdm_time(end-cp_length+1:end, :); ofdm_time];
ofdm_signal = ofdm_signal(:);
% ==========================================================================

figure;
plot(real(ofdm_grid), imag(ofdm_grid), '.', ...
    'MarkerSize', 30, 'Color', [1 0.341 0.341]);
axlim = (max(max(abs(ofdm_grid)))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% Multipath Propagation (ISI)
ofdm_signal_isi = ofdm_signal + 0.5 * circshift(ofdm_signal, 10);

%% OFDM Demodulation
rx = reshape(ofdm_signal_isi, n_subcarriers + cp_length, []);
rx = rx(cp_length+1:end, :);
ofdm_grid_rec = fft(rx, n_subcarriers, 1);
colors = parula(qam_order); colors = colors(randperm(qam_order), :);
figure;
hold on;
for i = 1 : length(ofdm_grid_rec(:))
    plot(real(ofdm_grid_rec(i)), imag(ofdm_grid_rec(i)), ...
        '.', 'Color', colors(qam_symbols(i) + 1, :), 'MarkerSize', 30);
end
hold off;
axlim = max(max(abs(ofdm_grid_rec))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% Symbol Error Rate Calculation
qam_symbols_rec = qamdemod(ofdm_grid_rec .* scale, qam_order, 0, 'gray');
n_errors = sum(sum(qam_symbols_rec ~= qam_symbols));
SER = n_errors / (n_subcarriers * n_symbols);
fprintf('Symbol Error Rate: %g\n', SER);
