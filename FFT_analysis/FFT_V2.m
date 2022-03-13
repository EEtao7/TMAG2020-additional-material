function [frequency, amplitude] = FFT_V2(signal, sample_frequency, base_frequency, number_of_circle, max_frequency)
% FFT_V2 is a better function for FFT analysis
Fs = sample_frequency;
Fb = base_frequency; % if base frequency is 1, we just use the whole window of the signals
if (Fb == 1)
    size_of_signal = length(signal);
    Y = fft(signal);
    P2 = abs(Y/size_of_signal);
    P1 = P2(1: round(size_of_signal/2)+1);
    P1(2 : end-1) = 2*P1(2 : end-1);
    if(max_frequency+1 <= length(P1))
        amplitude = P1(1: max_frequency+1);
        f = Fs*(0 : round(size_of_signal/2))/size_of_signal;
        frequency = f(1: max_frequency+1);
    else
        disp('size of signal is not enough for the max frequency')
    end
else
    size_of_signal = length(signal);
    size_of_window = number_of_circle*round(Fs/Fb);
    if(size_of_signal>=size_of_window)
        signal_of_window = signal(end-size_of_window+1 :  end); % in most situation, we will find the signal will have the steady waveform in the end
        Y = fft(signal_of_window);
        P2 = abs(Y/size_of_window);
        P1 = P2(1: round(size_of_window/2)+1);
        P1(2 : end-1) = 2*P1(2 : end-1);
        if(ceil(max_frequency/Fb)+1 <= length(P1))
            amplitude = P1(1: ceil(max_frequency/Fb)+1);
            f = Fs*(0 : round(size_of_window/2))/size_of_window;
            frequency = f(1: ceil(max_frequency/Fb)+1);
        else
            disp('size of signal is not enough for the max frequency')
        end
    else
        disp('size of signal is not enough for this base frequency')
    end
end
end