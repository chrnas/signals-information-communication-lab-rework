function [zI, zQ, A, tau] = receiver(y)
    format bank
    %% Initiate variables
    f_low = 20000;     % Lower cutoff frequency
    f_high = 40000;    % Higher cutoff frequency
    fs_low = 20000;        % Sampling frequency
    upsampling_factor = 20;
    fc = (f_low + f_high) / 2; % Carrier frequency
    bandwidth = fs_low/2;
    fs_high = fs_low * upsampling_factor;     % Upsampled sampling frequency
    Ts_low = 1 / fs_low;        % Sampling period time
    Ts_high = 1 / fs_high;  % Upsampled sampling period time
    order = 200;
    %% Create bandpass filter
    bpf_margin = 0.9;
    [b,a] = fir1(order,[fc-bandwidth*bpf_margin,fc+bandwidth*bpf_margin]/(fs_high/2),"bandpass");
    %% Apply bandpass filter
    y=filter(b,a,[y;zeros(order/2,1)]);  % Filter signal
    %% Remove zeroes from y
    y=y(order/2+1:end);

    %% Create Chirp-signal
    margin = 0.0001;
    chirp_lower_limit = fc-(bandwidth/2)*(margin);
    chirp_upper_limit = fc+(bandwidth/2)*(margin);
    chirp_end = 1 - Ts_high;
    chirp_time = (0:Ts_high:chirp_end)';
    chirp_signal = chirp(chirp_time, chirp_lower_limit, chirp_end, chirp_upper_limit);
    %% Correllation
    [corr, lag] = xcorr(y, chirp_signal);
    %% Get max peak and its index
    [max_peak, max_peak_index] = max((abs(corr)));
    %% Adjust tau based on where the max peak is by getting it from lag
    tau = (lag(max_peak_index)/fs_high);
    %% Add zeros to the end of y based on how big tau is
    y = [y;zeros(tau*fs_high,1)];
    %% Cut of the delay and chirp in the start of the signal
    y = y(length(chirp_signal) + tau*fs_high+ 1:end);
    %% Get the amplitude of the chip signal
    autocorrelation = xcorr(chirp_signal);
    [max_val_auto, max_index_auto] = max(abs(autocorrelation));
    A = (corr(max_peak_index)/autocorrelation(max_index_auto));
    y = y.*(1/(A));
    %% Demodulate the signal
    t = Ts_high*(0:(length(y)-1)).'; % get the new time vector
    % Apply the carrier for demodulation
    yI = y.*(2*cos(2*pi*fc*t));
    yQ = y.*(-2*sin(2*pi*fc*t));
    % Create the LP filter
    f_norm = (f_low/2)/(f_high/2)/2; 
    [b,a] = fir1(order, 1/upsampling_factor, "low"); 
    % Apply the filter and add zeroes to the end of the signal
    yI = filter(b, a, [yI;zeros(order/2, 1)]);
    yQ = filter(b, a, [yQ;zeros(order/2, 1)]);
    % Remove the zeroes at the start of the signal
    yI_demodulated = yI(order/2+1:end);
    yQ_demodulated = yQ(order/2+1:end);
    %% Downsample
    zI = downsample(yI_demodulated, upsampling_factor);
    zQ = downsample(yQ_demodulated, upsampling_factor);

    %% Adjust the size of the return signals
    zI = zI(1:100000);
    zQ = zQ(1:100000);

    %% Return values with correct output
    tau = tau*1000000; %Convert to mu s
    A = round(A,1);
    %% Print outut
    fprintf('%s %.1f %s\n', 'tau = ', tau, 'mu s')
    fprintf('%s %.1f\n', 'A = ', A);
end