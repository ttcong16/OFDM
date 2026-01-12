clearvars; close all;
color{1} = [35 60 230]/255;   % '#233ce6'
color{2} = [255 87 87]/255;   % '#ff5757'

%% Parameters
n_subcarriers = 128;
n_symbols = 20;
scs = 15e3;
data_order = 16;
pilot_order = 4;
cp_length = 20;
pilot_sc_loc = 4;
pilot_sym_loc = 4;
data{1} = 'HUST ';
data{2} = 'abcdefghijklmnopqrstuvwxyz ';

%% Data & Pilot Locations
pilot_locations{1} = false(n_subcarriers, n_symbols);
pilot_locations{1}(1 : pilot_sym_loc : end, 1 : pilot_sc_loc : end) = true;

pilot_locations{2} = false(n_subcarriers, n_symbols);
pilot_locations{2}(pilot_sc_loc / 2 + 1 : pilot_sc_loc : end, ...
    1 : pilot_sym_loc : end) = true;

data_locations = true(n_subcarriers, n_symbols);
data_locations(pilot_locations{1} | pilot_locations{2}) = false;

figure;
img = zeros(n_subcarriers, n_symbols);
img(data_locations) = 0;
img(pilot_locations{1}) = 1;
img(pilot_locations{2}) = 2;
cmap = [
    255 255 255
    35  60  230
    255 87  87
]/255; % MATLAB 2015: colormap double [0..1]
imagesc(img);
colormap(cmap);

%% Create QAM Symbols
n_ports = 2;
data_bits = cell(1, n_ports);
data_states = cell(1, n_ports);
pilot_states = cell(1, n_ports);

k_data  = log2(data_order);

% MATLAB 2015: emulate 'UnitAveragePower', true by scaling
ref_data  = qammod((0:data_order-1).',  data_order,  0, 'gray');
ref_pilot = qammod((0:pilot_order-1).', pilot_order, 0, 'gray');
scale_data  = sqrt(mean(abs(ref_data ).^2));
scale_pilot = sqrt(mean(abs(ref_pilot).^2));

for p = 1 : n_ports
    % utils.text2bits(data{p})
    u = uint8(data{p}(:));
    b = de2bi(u, 8, 'left-msb');      % N x 8
    bitstream = b.';                 % 8 x N
    bitstream = bitstream(:);        % column bits

    % utils.repeat_bits(bitstream, data_order, sum(sum(data_locations)))
    nQamSyms = sum(sum(data_locations));
    need = nQamSyms * k_data;
    rep = ceil(need / length(bitstream));
    data_bits{p} = repmat(bitstream, rep, 1);
    data_bits{p} = data_bits{p}(1:need);

    % qammod(...,'InputType','bit','UnitAveragePower',true) equivalent
    data_int = bi2de(reshape(data_bits{p}(:), k_data, []).', 'left-msb');
    data_states{p} = qammod(data_int, data_order, 0, 'gray') ./ scale_data;

    pilot_symbols = randi(pilot_order, n_subcarriers / pilot_sc_loc, ...
        n_symbols / pilot_sym_loc) - 1;

    % qammod(...,'UnitAveragePower',true) equivalent
    pilot_states{p} = qammod(pilot_symbols, pilot_order, 0, 'gray') ./ scale_pilot;
end

%% OFDM Modulation
tx_grid = cell(1, n_ports);
tx_signal = cell(1, n_ports);
for p = 1 : n_ports
    tx_grid{p} = zeros(n_subcarriers, n_symbols);
    tx_grid{p}(data_locations) = data_states{p};
    tx_grid{p}(pilot_locations{p}) = pilot_states{p};

    % ofdmmod(...) equivalent (matches their OFDMGrid convention)
    tx_time = ifft(ifftshift(tx_grid{p}, 1), n_subcarriers, 1);
    tx_cp   = [tx_time(end-cp_length+1:end,:); tx_time];
    tx_signal{p} = tx_cp(:);

    figure;
    plot(real(tx_grid{p}(data_locations)), ...
        imag(tx_grid{p}(data_locations)), 'o', 'Color', color{p});
    axlim = max(max(abs(tx_grid{p}))) + 0.05;
    ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
    xlabel('In-phase'); ylabel('Quadrature');
end

%% Channel
fs = n_subcarriers * scs;

% MATLAB 2015: build with properties that exist, then call with step()
mimoChan = comm.MIMOChannel('SampleRate', fs, ...
    'PathDelays', [0 1.5e-6], ...
    'AveragePathGains', [0 -5], ...
    'MaximumDopplerShift', 5, ...
    'RandomStream', 'mt19937ar with seed', ...
    'Seed', 12);

% Some MATLAB 2015 configs infer antenna counts from input; keep it simple.
tx_mat = [tx_signal{1}, tx_signal{2}];

rx_signal = step(mimoChan, tx_mat);   % MATLAB 2015: DO NOT use mimoChan(tx_mat)

%% OFDM Demodulation
rx_grid = cell(1, n_ports);
for p = 1 : n_ports
    % ofdmdemod(...) equivalent (inverse of the mod above)
    rx_mat = reshape(rx_signal(:, p), n_subcarriers + cp_length, []);
    rx_mat = rx_mat(cp_length+1:end, :);
    rx_grid{p} = fftshift(fft(rx_mat, n_subcarriers, 1), 1);

    figure;
    plot(real(rx_grid{p}(data_locations)), ...
        imag(rx_grid{p}(data_locations)), 'o', 'Color', color{p});
    axlim = max(max(abs(rx_grid{p}))) + 0.05;
    ylim([-axlim, axlim]); xlim([-axlim, axlim]); axis square;
    xlabel('In-phase'); ylabel('Quadrature');
end

%% Channel Estimation
hest = cell(n_ports, n_ports);
[X{1}, Y{1}] = meshgrid(1 : pilot_sym_loc : n_symbols, ...
    1 : pilot_sc_loc : n_subcarriers);
[X{2}, Y{2}] = meshgrid(1 : pilot_sym_loc : n_symbols, ...
    1 + pilot_sc_loc / 2 : pilot_sc_loc : n_subcarriers);
[Xq, Yq] = meshgrid(1 : n_symbols, 1 : n_subcarriers);

for p_tx = 1 : n_ports
    for p_rx = 1 : n_ports
        hest_pilots = reshape(rx_grid{p_rx}(pilot_locations{p_tx}), ...
            n_subcarriers / pilot_sc_loc, n_symbols / pilot_sym_loc) ...
            ./ pilot_states{p_tx};
        hest{p_tx}{p_rx} = interp2(X{p_tx}, Y{p_tx}, hest_pilots, ...
            Xq, Yq, 'spline');
    end
end

%% Equalize Data
est_grid = zeros(n_subcarriers, n_symbols, 2);
for k = 1 : n_subcarriers
    for t = 1 : n_symbols
        H = [
            hest{1}{1}(k, t), hest{2}{1}(k, t);
            hest{1}{2}(k, t), hest{2}{2}(k, t)
        ];
        y = [rx_grid{1}(k, t); rx_grid{2}(k, t)];
        est_grid(k, t, :) = pinv(H) * y;
    end
end

for p = 1 : n_ports
    figure;
    plot(real(est_grid(:, :, p)), imag(est_grid(:, :, p)), ...
        'o', 'Color', color{p});
    axlim = max(max(abs(est_grid(:, :, p)))) + 0.05;
    ylim([-axlim, axlim]); xlim([-axlim axlim]); axis square;
    xlabel('In-phase'); ylabel('Quadrature');
end

%% QAM Demodulation
for p = 1 : n_ports
    data_states_rx = est_grid(:, :, p);
    data_states_rx = data_states_rx(data_locations);

    % qamdemod(...,'OutputType','bit','UnitAveragePower',true) equivalent
    data_int_rx = qamdemod(data_states_rx .* scale_data, data_order, 0, 'gray');
    data_bits_rx = de2bi(data_int_rx(:), k_data, 'left-msb').';
    data_bits_rx = data_bits_rx(:);

    % utils.bits2text(data_bits_rx)
    n8 = floor(length(data_bits_rx)/8);
    bb = reshape(data_bits_rx(1:8*n8), 8, n8).';
    uu = bi2de(bb, 'left-msb');
    data_rx = char(uu.').';

    n_errors = sum(data_bits_rx ~= data_bits{p}(:));
    BER = n_errors / length(data_bits_rx);

    fprintf('Port %d\n', p);
    fprintf(' BER: %g\n', BER);
    fprintf(' Data: %s\n', data_rx);
end
