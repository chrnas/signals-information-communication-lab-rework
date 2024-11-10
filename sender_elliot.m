function x = sender_elliot(xI,xQ)
    %% Konstanter
    f1 = 120000;                                                                 %Undre gräns
    f2 = 140000;                                                                 %Övre gräns
    fc = (f1+f2)/2;                                                             %Bärfrekvens
    M = 20;                                                                     %Upp-/nedsamplingsfaktor
    fs = 20e3;                                                                  %Sampelfrekvens
    fs_up = fs*M;                                                               %fs uppsamplad
    N = 500;                                                                    %Filtrets gradtal

    %Läs in ljudfilerna
    %[xI,fs]=audioread("xI.wav");
    %[xQ,fs]=audioread("xQ.wav");
    
    %Uppsampla i tidsled
    xI = upsample(xI, M);                                                     %Uppsampla xI
    xQ = upsample(xQ, M);                                                     %Uppsampla xQ
    
            
    %% Filtrera i tidsled
    %Skapa filter
    f0 = 1/M;                                                                   %Normerad gränsfrekvens
    [b,a] = fir1(N, f0);                                                        %Skapa filter
    
    %Filtrera xI
    xI = filter(b, a, [xI;zeros(N/2, 1)]);                                      %Filtrera & justera signal xI
    xI = xI(N/2+1:end);                                                         %(Från lecturecommands)
    
    %Filterar xQ
    xQ = filter(b, a, [xQ;zeros(N/2, 1)]);                                      %Filtrera & justera signal xQ
    xQ = xQ(N/2+1:end);                                                         %(Från lecturecommands)


    %% Skapa tidsvektor
    L= length(xI);
    Ts = 1/fs;                                                                  %Sampleperiod
    t = (Ts/M)*[0:L-1];
    

    %% I/Q representation
    XIQ = xI.*cos(2*pi*fc*t).'-xQ.*sin(2*pi*fc*t).';                            %(Boken S.28)
    x = XIQ;

    %% Generera chirp
    time = 0:(1/fs_up):1;
    start_freq = 105e3;                                                         %Startfrekvens
    end_freq = 115e3;                                                           %Slutfrekvens
    chirped = chirp(time, start_freq, 1, end_freq)';                            %Generera chirp
    x = cat(1, chirped, x);                                                     %Slå ihop uprsprungssignal med chrip
end