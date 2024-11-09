function x = sender(xQ, xI)

    %% Initiate variables    
    f_low = 80000;  % Lower cutoff frequency
    f_high = 100000; % Higher cutoff frequency
    fs_low = 20000; % Sampling frequency
    upsampling_factor = 20;
    fc = (f_low + f_high) / 2;  % Carrier frequency
    bandwidth = fs_low/2;
    margin = 0.0001;
    fs_high = fs_low * upsampling_factor;   % Upsampled sampling frequency
    Ts_low = 1 / fs_low;    % Sampling period time
    Ts_high = 1 / fs_high;  % Upsampled sampling period time
    order = 500;
    
    %% Upsample xQ
    xQ_upsampled = upsample(xQ, upsampling_factor);
    
    %% Lowpass-filtrering xQ
    [b, a] = fir1(order, 1 / upsampling_factor, "low");
    xQ_upsampled = filter(b, a, [xQ_upsampled; zeros(order / 2, 1)]);
    
    %% Remove zeroes xQ
    xQ_upsampled = xQ_upsampled((order / 2 + 1):end);
    
    %% Modulering xQ
    t = Ts_high * (0:(length(xQ_upsampled) - 1)).';
    carrier_Q = -sin(2 * pi * fc * t);
    xQ_modulated = xQ_upsampled .* carrier_Q*upsampling_factor;
    
    %% Upsample xI
    xI_upsampled = upsample(xI, upsampling_factor);
    
    %% Lowpass-filtrering xI
    [b, a] = fir1(order, 1 / upsampling_factor, "low");
    xI_upsampled = filter(b, a, [xI_upsampled; zeros(order / 2, 1)]);
    
    %% Remove zeroes xI
    xI_upsampled = xI_upsampled((order / 2 + 1):end);
    
    %% Modulate xI
    t = Ts_high * (0:(length(xQ_upsampled)-1)).';
    carrier_I = cos(2*pi*fc*t);
    xI_modulated = xI_upsampled .* carrier_I *upsampling_factor;
    
    %% Combine modulated signals
    x_no_chirp = xI_modulated + xQ_modulated;

    %% Create Chirp-signal
    chirp_lower_limit = fc-(bandwidth/2)*(margin);
    chirp_upper_limit = fc+(bandwidth/2)*(margin);
    chirp_end = 1 - Ts_high;
    chirp_time = (0:Ts_high:chirp_end)';
    chirp_signal = chirp(chirp_time, chirp_lower_limit, chirp_end, chirp_upper_limit);
    %% Add chirp to x-signal
    x = cat(1, chirp_signal, x_no_chirp);
end





