clearvars; close all;

%% Parameters
n_subcarriers = 128;
n_symbols = 20;
scs = 15e3;
data_order = 256;
pilot_order = 4;
cp_length = 12;

fs = n_subcarriers * scs;

%% OFDM Grid Generation & Modulation
data_symbols = randi(data_order, n_subcarriers, n_symbols) - 1;

% MATLAB 2015: legacy qammod syntax
data_qam_states = qammod(data_symbols, data_order, 0, 'gray');

pilot_symbols = randi(pilot_order, n_subcarriers / 4, n_symbols / 4) - 1;
pilot_qam_states = qammod(pilot_symbols, pilot_order, 0, 'gray');

% --- Equivalent to 'UnitAveragePower', true (manual normalization) ---
data_ref  = qammod((0:data_order-1).',  data_order,  0, 'gray');
pilot_ref = qammod((0:pilot_order-1).', pilot_order, 0, 'gray');
data_scale  = sqrt(mean(abs(data_ref ).^2));
pilot_scale = sqrt(mean(abs(pilot_ref).^2));

data_qam_states  = data_qam_states  ./ data_scale;
pilot_qam_states = pilot_qam_states ./ pilot_scale;
% -------------------------------------------------------------------

pilot_locations = false(n_subcarriers, n_symbols);
pilot_locations(1 : 4 : end, 1 : 4 : end) = true;
data_locations = true(n_subcarriers, n_symbols);
data_locations(pilot_locations) = false;

ofdm_grid = zeros(n_subcarriers, n_symbols);
ofdm_grid(data_locations) = data_qam_states(data_locations);
ofdm_grid(pilot_locations) = pilot_qam_states;

%% OFDM Modulation (replace ofdmmod for MATLAB 2015)
tx_time = ifft(ofdm_grid, n_subcarriers, 1);
tx_cp   = [tx_time(end-cp_length+1:end,:); tx_time];
tx_waveform = tx_cp(:);

blue = [35 60 230]/255;

figure;
plot(real(ofdm_grid(data_locations)), imag(ofdm_grid(data_locations)), ...
    'o', 'Color', blue);
axlim = max(max(abs(ofdm_grid(data_locations)))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

figure;
plot(real(ofdm_grid(pilot_locations)), imag(ofdm_grid(pilot_locations)), ...
    'o', 'Color', blue);
axlim = max(max(abs(ofdm_grid(pilot_locations)))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% Multipath Propagation (gi? y h?t code h?)
time = (0 : length(tx_waveform) - 1) / fs;
rx_waveform = tx_waveform ...
    + 0.25 * circshift(tx_waveform, 10) .* exp(2 * pi * 100 * time.');

%% OFDM Demodulation (replace ofdmdemod for MATLAB 2015)
rx_mat = reshape(rx_waveform, n_subcarriers + cp_length, []);
rx_mat = rx_mat(cp_length + 1:end, :);
ofdm_grid_rec = fft(rx_mat, n_subcarriers, 1);

figure;
plot(real(ofdm_grid_rec(data_locations)), ...
    imag(ofdm_grid_rec(data_locations)), 'o', 'Color', blue);
axlim = max(max(abs(ofdm_grid_rec))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% Channel Estimation
pilot_qam_states_rec = ofdm_grid_rec(pilot_locations);
hest = reshape(pilot_qam_states_rec, n_subcarriers / 4, n_symbols / 4) ...
    ./ pilot_qam_states;
[X, Y] = meshgrid(1 : 4 : n_symbols, 1 : 4 : n_subcarriers);
[Xq, Yq] = meshgrid(1 : n_symbols, 1 : n_subcarriers);
hest_interp = interp2(X, Y, hest, Xq, Yq, 'spline');

ofdm_grid_rec = ofdm_grid_rec ./ hest_interp;

figure;
plot(real(ofdm_grid_rec(data_locations)), ...
    imag(ofdm_grid_rec(data_locations)), 'o', 'Color', blue);
axlim = max(max(abs(ofdm_grid_rec))) + 0.05;
ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
xlabel('In-phase'); ylabel('Quadrature');

%% Symbol Error Rate Calculation
% scale back before legacy qamdemod to match UnitAveragePower=true behavior
data_symbols_rec = qamdemod(ofdm_grid_rec(data_locations) .* data_scale, ...
    data_order, 0, 'gray');

n_errors = sum(sum(data_symbols_rec ~= data_symbols(data_locations)));
SER = n_errors / length(data_symbols_rec);
disp(['Symbol Error Rate: ' num2str(SER)]);
