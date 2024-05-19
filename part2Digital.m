    % Parameters
    pkg load signal
    bit_rate = 1000;
    bit_duration = 1 / bit_rate;
    num_bits = 100;
    amplitude = 1;

    % generate random bits

    bits = randi([0, 1], 1, num_bits);

    % generate time vector
    t = linspace(0, num_bits * bit_duration, num_bits * 100);


    % generate NRZ signal
    signal = zeros(1, length(t));

    for i = 1:num_bits
        if bits(i) == 1
            signal((i-1)*100+1 : i*100) = amplitude;
        end
    end

    % plot
    plot(t, signal);
    xlabel('Time in seconds');
    ylabel('Amplitude');
    title('Unipolar NRZ Signal');
    ylim([-1.5, 1.5]);
    grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    carrier_frequency = 2000;

    % generate ASK signal
    signal = zeros(1, length(t));

    for i = 1:num_bits
        if bits(i) == 1
            signal((i-1)*100 + 1:i*100) = amplitude * sin(2*pi*carrier_frequency*t((i-1)*100 + 1:i*100));
        end
    end

    % plot time domain signal
    figure;
    subplot(2,1,1);
    plot(t, signal);
    xlabel('Time in seconds');
    ylabel('Amplitude');
    title('Transmitted ASK Signal');
    ylim([-1.5, 1.5]);
    grid on;

    % plot freq domain signal

    subplot(2,1,2);
    fft_signal = fft(signal);
    f = linspace(0, bit_rate, length(fft_signal));
    plot(f, abs(fft_signal));
    xlabel('Frequency ');
    ylabel('Magnitude');
    title('Spectrum of Transmitted ASK Signal');
    grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Coherent ASK Receiver

    % Parameters
    oscillator_phase = [0, 30, 60, 90]; % Oscillator phase in degrees

    shift=1;

    phase = oscillator_phase(1,shift)
    %Demodulate received signal
    received_signal = signal .* sin(2*pi*carrier_frequency*t + deg2rad(phase));

    % plot time domain signal
    figure;
    subplot(2,1,1);
    plot(t, received_signal);
    xlabel('Time in seconds');
    ylabel('Amplitude');
    title(['Received ASK signal with Phase shift', num2str(phase), ' degrees']);
     ylim([-1.5, 1.5]);
    grid on;

    % plot freq domain signal
    subplot(2,1,2);
    fft_received_signal_0 = fft(received_signal);
    plot(f, abs(fft_received_signal_0));
    xlabel('Frequency ');
    ylabel('Magnitude');
    title(['Spectrum of Received ASK Signal with Phase shift ', num2str(phase), ' degrees']);
    grid on;




    pkg load signal
    % low pass filter

    % Define parameters
    fc = 60; % Cutoff frequency
    fs = 10000; % Sampling frequency
    order = 5; % Filter order

    % Design the Butterworth low-pass filter
    [b, a] = butter(order, fc/(fs/2), "low");


    x=received_signal

    % Apply the filter
    low_pass_signal = filter(b, a, received_signal);

    % Plot the original and filtered signals
    figure;
    subplot(2,1,1);
    plot(t, received_signal)
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Spectrum of Received ASK Signal']);
     ylim([-1.5, 1.5]);
    grid on;

    subplot(2,1,2);
    plot(t, low_pass_signal)
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Low-pass Filtered Signal');
     ylim([-1.5, 1.5]);
    grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compartor implementation
    threshold=0.2
    final_bits=zeros(1,length(low_pass_signal))

    for j = 1:length(low_pass_signal)

        element = low_pass_signal(j);

       if(element>threshold)
        final_bits(j)=1;
       else
       final_bits(j)=0 ;
       end
    end
    figure;
    plot(t,final_bits);
    ylim([-1.5, 1.5]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('demodulated signal');
    grid on;





