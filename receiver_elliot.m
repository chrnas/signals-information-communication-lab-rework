function [zI,zQ,A,tau] = receiver_elliot(x)
    %% Definiera konstanter
    f1 = 100e3;                                                                 %Undre gräns, från uppgift
    f2 = 120e3;                                                                 %Övre gräns, från uppgift
    fs = 20e3;                                                                  %Sampelfrekvens
    fc = (f1+f2)/2;                                                             %Bärfrekvens
    M = 20;                                                                     %Upp-/nedsamplingsfaktor
    fs_up = fs*M;                                                               %Upsamplingsfrekvens
    nyquist = fs_up/2;                                                          %Nyquistfrekvens 
    N = 500;                                                                    %Filtrets gradtal

    %% FIR-filtrering
    lowerFreq = (f1/nyquist);                                                   %Undre gräns
    higherFreq = (f2/nyquist);                                                  %Övre gräns
    [b,a] = fir1(N, [lowerFreq, higherFreq]);                                   %Passbandet
    x = filter(b, a, x);                                                        %Filtrerad signal
    x = x((N/2+1):end);                                                         %Justerad signal

    %% Identifiera fördröjningen
    %Generera chirp
    start_freq = 105e3;
    end_freq = 115e3;
    time = 0:(1/(fs_up)):1;                                                     %Skapa tidsvektor
    y = chirp(time, start_freq, 1, end_freq)';                                  %Generera chirp
    
    %Korskorrelation
    [values_cross, values_sample] = xcorr(x, y);                                %Korskorrelerar in- och utsignal
   
    %Maximal Korskorrelation
    [max_values_cross, index] = max(abs(values_cross));                         %Hitta maximala absolutvärdet
    
    %Tidsfördröjning    
    delayed_sample = values_sample(index);                                      %Hämta fördröjningsvärdet
    tau = (delayed_sample/400e3);                                               %Ta fram tau
    tau = tau * 1e6                                                             %Korriger för fördröjning
    
    %% Identifiera amplitudskalning
    %Autokorrelation
    [auto_values_cross, values_sample] = xcorr(y);                              %Autokorrelera
    
    %Skalning
    max_auto = max(auto_values_cross);                                          %Hitta maximalt värde i auto_values_cross
    A = max_values_cross/max_auto;                                              %Beräkna skalfaktorn
    if (max(values_cross) < abs(min(values_cross)))                             %Ta hänsyn till negativa värden på A
        A = -A;
    end
    x = x/A;                                                                    %Normera signalen
    A = round(A,1)                                                              %Runda av enligt instruktion, abvänds för att se resultat
        

    %% Eliminera Chirp
    x = x((length(y)+1):end);
    x = x((delayed_sample+1):end);                                              %Ta bort första sekunden av signalen (där chirp:en finns)


    %% Tidsvektor
    L= length(x);
    Ts = 1/fs;
    Ts_div = Ts/M;
    t = Ts_div*(0:L-1);
    
    
    %% I/Q-demodulation
    zI = x.*(2*cos(2*pi*fc*t)).';                                               %Demodulera enligt boken
    zQ = x.*((-2)*sin(2*pi*fc*t)).';                                            %Demodulera enligt boken
    
    
    %% FIR-filtrering
    F0 = 1/M;                                                                   %Normerad gränsfrekvens
    [b,a] = fir1(N, F0);                                                        %Skapa filter
    
    zI = filter(b, a, [zI;zeros(N/2, 1)]);                                      %Filtrera och justera signal 
    zI = zI(N/2+1:end);

    zQ = filter(b, a, [zQ;zeros(N/2, 1)]);                                      %Filtrera och justera signal 
    zQ = zQ(N/2+1:end);
   

    %% Nedsampling
    %Nedsampla
    zI = downsample(zI, M);                                                     
    zQ = downsample(zQ, M);
    
    %Skär ut element 1-100 000 (så vi kan plotta på ett bra sätt)
    zI = zI(1:1e5);
    zQ = zQ(1:1e5);
end