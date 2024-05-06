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



# q1 plot function x(t)
%adding linear part
%getting indices
m=zeros(size(t));
less_neg1=(t<-1) & (t>=-2);
more_pos1=(t>1)&(t<=2);
m(less_neg1)=t(less_neg1)+2*ones(size(t(less_neg1))); %2+t graph
m(more_pos1)=-t(more_pos1)+2*ones(size(t(more_pos1))); %2-t graph
m((t>=-1) & (t<=1))=ones(size(t(t>=-1&t<=1)));


figure(1);
plot(t,m);
xlabel("Time (sec)");
ylabel("x(t)");
box off;


# q2
% manual fft
%analytic_M=zeros(size(t));
%for j=1:length(t)
%    for k=1:length(t)
%      F(j)=t(k)*exp(-1i*2*pi*(j-1)*(k-1)/N);
%      analytic_M(j)=analytic_M(j)+F(j);
%    end
%end
%figure(2)
%plot(f,analytic_M,"--r");

# q3
# fourier transform of m(t) --> M(f)
M=fftshift(fft(m,N)/N); %since it's an even function scaled by 1/N
figure(2)
hold on;
stem(f,M,"-b");
xlabel("Frequency (Hz)");
ylabel("M(f)");
legend("Using FFT");
box off;

% q4: Band Width
Energy_from_freq=sum(abs(M).^2)*df;
# by looking at both variables values in workspace we can verify parseval's theorem.

%%%%%% Finding Band Width %%%%%%%%%%
# we first have to find f=0
index_f0=find(f==0);
# Initializing accumulate energy
Energy_acc=0;
for index_f= index_f0:length(f)
  # we multiply by 2 to use the property of even function
  Energy_acc+=2*(abs(M(index_f)).^2)*df;
  # comparing it to 95% of total energy
  if(Energy_acc>=0.95*Energy_from_freq)
    BW=f(index_f);
    break
  end
end


