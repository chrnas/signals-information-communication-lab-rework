function [zI, zQ, A, tau] = receiver(y)
    format bank
    %% Initiate variables
    f_low = 80000;     % Lower cutoff frequency
    f_high = 100000;    % Higher cutoff frequency
    fs_low = 20000;        % Sampling frequency
    upsampling_factor = 20;
    fc = (f_low + f_high) / 2; % Carrier frequency
    bandwidth = fs_low/2;
    fs_high = fs_low * upsampling_factor;     % Upsampled sampling frequency
    Ts_low = 1 / fs_low;        % Sampling period time
    Ts_high = 1 / fs_high;  % Upsampled sampling period time
    order = 500;
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
    
    %% Identify time delay(tau)
    delayed_sample = lag(max_peak_index);
    tau = (delayed_sample/fs_high);

    %% Identify amplitude scale factor
    [auto_corr, auto_lag] = xcorr(chirp_signal); %Autocorrelation
    max_auto = max(auto_corr);
    A = max_peak/max_auto; 
    if (max(corr) < abs(min(corr)))
        A = -A;
    end
    
    %% Normalize the signal
    y = y/A;

    %% Eliminate chirp signal and time delay
    y = y((length(chirp_signal)+1):end); %Cut off the chirp signal
    y = y((delayed_sample+1):end); %Cut off the delay                   
    %y = [y;zeros(tau*fs_high,1)];

    %% Demodulate the signal
    t = Ts_high*(0:(length(y)-1)).'; % get the new time vector
    
    % Apply the carrier for demodulation
    yI = y.*(2*cos(2*pi*fc*t));
    yQ = y.*(-2*sin(2*pi*fc*t));
    
    % Create the LP filter
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
    A = round(A,1);
    tau = tau*1000000; %Convert to u s
end