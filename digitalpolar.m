L = 10 ;%resolution 
Fs = 8*L;
bitStreamLength = 64;
bitStream = randi([0,1],1,bitStreamLength);
data = rand(1,64)>0.5;
datacode_uni = reshape(repmat(data,L,1),1,length(data)*L);
%same for polar
datacode_poler = reshape(repmat(data,L,1),1,length(data)*L)-0.5;
%
subplot(3,1,1);plot(datacode_uni);
xlabel('bit rate');
ylabel('A');
title('unipolar');
subplot(3,1,2);plot(datacode_poler);
xlabel('bit rate');
ylabel('A');
title('polar');

N = length(datacode_uni);
unipolar_nrz_l_freq = fftshift(fft(datacode_uni)); 
magnitude_spectrum_nrz_l = abs(unipolar_nrz_l_freq);
f = (-N/2:N/2 - 1)*(Fs/N);%??

figure;
plot(f, magnitude_spectrum_nrz_l);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectral Domain of Unipolar NRZ-L Signal');

polar_nrz_l_freq = fftshift(fft(datacode_poler)); 
magnitude_spectrum_nrz_l = abs(polar_nrz_l_freq);
figure;
plot(f, magnitude_spectrum_nrz_l);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectral Domain of polar NRZ-L Signal');

 
% Manchester Encoding
manchesterEncoded = [];
for i = 1:bitStreamLength
    if bitStream(i) == 0
        manchesterEncoded = [manchesterEncoded, ones(1, 100)*1, ones(1, 100)*(-1)]; 
    else
        manchesterEncoded = [manchesterEncoded, ones(1, 100)*(-1), ones(1, 100)*1]; 
    end
end

fs = 100;

% Plot Time Domain
t = (0:length(manchesterEncoded)-1) / fs;
figure;
subplot(2,1,1);
plot(t, manchesterEncoded, 'b', 'LineWidth', 2);
title('Manchester Line Coding (Time Domain)');
xlabel('Time (seconds)');
ylabel('Amplitude');
ylim([-1.5 1.5]);
grid on;

% Plot Spectral Domain
N = length(manchesterEncoded);
f = (-N/2:N/2-1)*(fs/N);
spectrum = abs(fftshift(fft(manchesterEncoded)));
subplot(2,1,2);
plot(f, spectrum, 'b', 'LineWidth', 2);
title('Manchester Line Coding (Spectral Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
