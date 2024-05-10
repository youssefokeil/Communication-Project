# clear command window
clc;
# remove all variables from workspace
clear
# close all figures
close all
# importing signal package
pkg load signal


%%%%%% Frequencies initialization %%%%%%
# carrier frequency initialization
fc=1e4; #10KHz

%%%%%% Frequency domain %%%%%%
# df=fs/N
df=0.01; #resolution 0.01Hz
# sampling frequency, inverse of sampling tinme
fs=100; #100Hz

%%%%%% Time and frequency definition %%%%%%
# time step
ts= 1/fs;
# T is the simulation time
T=1/df;
# Number of time samples
N=ceil(T/ts);
# time vector
t= -N*ts/2:ts:N*ts/2-ts;

# frequency vector
if(rem(N,2)==0) #even
  f= -(0.5*fs): df : (0.5*fs-df);
else #odd
  f= -(0.5*fs-0.5*df): df : (0.5*fs-0.5*df);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# q1 plot function x(t)
%adding linear part
%getting indices
x=zeros(size(t));
less_neg1=(t<-1) & (t>=-2);
more_pos1=(t>1)&(t<=2);
x(less_neg1)=t(less_neg1)+2*ones(size(t(less_neg1))); %2+t graph
x(more_pos1)=-t(more_pos1)+2*ones(size(t(more_pos1))); %2-t graph
x((t>=-1) & (t<=1))=ones(size(t(t>=-1&t<=1)));


figure(1);
plot(t,x);
xlabel("Time (sec)");
ylabel("x(t)");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# q2
% manual fft
analytic_X=sinc(f).*sin(1.5*2*pi*f)./(2*pi*f);
figure(2);
stem(f,abs(analytic_X)/max(analytic_X),"--r");


# q3
# fourier transform of m(t) --> M(f)
X=fftshift(fft(x,N)*ts); %non-periodic signals
figure(2)
hold on;
stem(f,abs(X)/max(X),"-b");
xlabel("Frequency (Hz)");
ylabel("X(f)");
legend("Analytical Fourier Transform","Using FFT");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q4: Band Width
Energy_from_freq=sum(abs(X).^2)*df;
# by looking at both variables values in workspace we can verify parseval's theorem.

%%%%%% Finding Band Width %%%%%%%%%%
# we first have to find f=0
index_f0=find(f==0);
# Initializing accumulate energy
Energy_acc=0;
for index_f= index_f0:length(f)
  # we multiply by 2 to use the property of even function
  Energy_acc+=2*(abs(X(index_f)).^2)*df;
  # comparing it to 95% of total energy
  if(Energy_acc>=0.95*Energy_from_freq)
    BW=f(index_f);
    break
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Ideal LPF %%%%%%%%%%
%% q5: Ideal LPF 1Hz
lpf_1hz= abs(f)<1; % ones only if value less than 1hz
X_filtered_1hz=X.*lpf_1hz;
x_filtered_1hz=ifft(ifftshift(X_filtered_1hz)/ts);
%%%%% Plotting filtered vs original signal %%%%%%%%
figure(3);
plot(t, x_filtered_1hz, "--r");
xlabel("Time (sec)");
ylabel("x(t)");
hold on;
plot(t,x,"-b");
legend("Signal after LPF 1Hz","Original Signal");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Ideal LPF %%%%%%%%%%
%% q6: Ideal LPF 0.3Hz
lpf= abs(f)<0.3; % ones only if value less than 0.3hz
X_filtered=X.*lpf;
x_filtered=ifft(ifftshift(X_filtered)/ts);
%%%%% Plotting filtered vs original signal %%%%%%%%
figure(4);
plot(t, x_filtered, "--r");
xlabel("Time (sec)");
ylabel("x(t)");
hold on;
plot(t,x,"-b");
legend("Signal after LPF 0.3Hz","Original Signal");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q7 defining new function
close all;
m=cos(2*pi*1*t);
m((t>6)|(t<0))=0;

%%%%% plotting in time domain
figure(1);
plot(t,m)
box off;

%%%% plotting in frequency domain
M=fftshift(fft(m)*ts); % non-periodic signals
figure(2);
hold off;
plot(f,abs(M));
box off;

%%%%%%%%%%%%%% getting band width
Energy_total_freq= sum(abs(M).^2)*df;
 %%%%% Band Width %%%%%%%%
 # we first have to find f=0
index_f0=find(f==0);
# Initializing accumulate energy
Energy_acc=0;
for index_f= index_f0:length(f)
  # we multiply by 2 to use the property of even function
  Energy_acc+=2*(abs(M(index_f)).^2)*df;
  # comparing it to 95% of total energy
  if(Energy_acc>=0.95*Energy_total_freq)
    BW2=f(index_f);
    break
  end
end


%%%%%% LPF 1Hz  %%%%%%%%%%%%%%%%
M_filtered=M.*lpf_1hz;
m_filtered=ifft(ifftshift(M_filtered)/ts);
figure(3)
hold off;
plot(t,m_filtered,"--r");
hold on;
plot(t,m,"-b");
legend("Signal after LPF 1Hz","Original Signal");
box off;

%%%%%%%% LPF 0.3Hz %%%%%%%%%%%%
M_filtered=M.*lpf;
m_filtered=ifft(ifftshift(M_filtered)/ts);
figure(4)
hold off;
plot(t,m_filtered,"--r");
hold on;
plot(t,m,"-b");
legend("Signal after LPF 0.3Hz","Original Signal");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q8:
BW=3; #Each channel bandwidth is 3 Hz
# carrier frequency initialization
fc1=20; #20Hz first carrier frequency
fc2=25.5 # fc1+3+2.5 =25.5Hz second carrier frequency
c1=cos(2*pi*fc1*t);   # first carrier

%% q10 : value for c2
c2=cos(2*pi*fc2*t);    # second carrier

%% first signal
s1 = x_filtered_1hz.*c1;  %% Define signal s1 in time domain
figure(1);
plot(t,s1);
xlabel("Time (sec)");
ylabel("s1(t)");
box off;

%% second signal
s2 = m.*c2;  %% Define signal s2 in time domain
figure(2);
plot(t,s2);
xlabel("Time (sec)");
ylabel("s2(t)");
box off;


% q9 : SSB signal USB
  %% Filter
  S2=fftshift(fft(s2,N)*ts); %non-periodic signals
  H = zeros(size(f));
  H(f>fc2 & f<(fc2+BW))=1; %% +ve frequency
  H(f<-fc2 & f>-(fc2+BW))=1; %% -ve frequency
  %% Filter spectrum
  S2 = S2.*H;
  s2=real(ifft(ifftshift(S2) / ts));

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %q11 :  s(t)= s1(t) + s2(t)
  %% plot signal in time domain
  s=s1+s2;
  figure(3);
plot(t,s);
xlabel("Time (sec)");
ylabel("s(t)");
box off;

%% plot signal in frequency domain
  S=fftshift(fft(s,N)*ts);
  figure(4);
  plot(f,S);
xlabel("Frequency (HZ)");
ylabel("S(f)");
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% q12: coherent demodulator
  g1=s.*c1;
  G1 = fftshift(fft (g1)) *ts;
  %% LPF (ideal)
H = abs(f)< BW;
%% After LPF
%% retrieve first message
 x_recieved= real(ifft(ifftshift(H.*G1) /ts));
 figure(5);
plot(t,x_recieved/max(x_recieved),"--b");
hold on
plot(t,x_filtered_1hz/max(x_filtered_1hz),"-r");
xlabel("Time (sec)");
ylabel("x(t)");
legend("recieved message","input message");
box off;

%% retrieve second message
g2=s.*c2;
  G2 = fftshift(fft (g2)) *ts;
m_recieved= real(ifft(ifftshift(H.*G2) /ts));
figure(6);
plot(t,m_recieved/max(m_recieved),"--b");
hold on
plot(t,m/max(m),"-r");
xlabel("Time (sec)");
ylabel("x(t)");
legend("recieved message","input message");
box off;
