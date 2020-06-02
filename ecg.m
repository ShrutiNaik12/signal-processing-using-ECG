clear all
clc
load('202m.mat')
fs = 360;    % sampling frequency
%% extracting first 10 mins data and Plotting it
[tm,signal,Fs,siginfo]=rdmat('202m');   % 'wfdb-app-toolbox-0-9-9' Add to path
given_data_in_mins = (length(signal)/fs)/60;
% given_data_in_secs = round(length(signal)/(given_data_in_mins*60));
% S = (given_data_in_secs * 600);
C = round(length(signal)*10/given_data_in_mins);
data= signal(1:C,:);               % first 10 mins data
figure;
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('10 mins ECG data'); % plotting the data
% figure
% plot(signal(214000:C,:));                                                          % small section of the data
%% Number of heart beats present in the data
[peakVal, locVal] = findpeaks(data,'MinPeakHeight',0.5,'MinPeakDistance',10); % detecting the highest peak
No_of_heart_beats = length(peakVal);                                          % No. of heartbeats in 10mins
S=['The data has ',num2str(No_of_heart_beats),' heart beats'];
disp(S)
%% RR-Interval
for k = 2:length(locVal)
    B = locVal(k,1)- locVal(k-1,1);
    RR_Interval(k,1) = B;
end
RR_Interval = RR_Interval(2:end,:); % RR-Interval durations
e=interp1(locVal,RR_Interval,'spline'); 
figure;
plot(RR_Interval); xlabel('time in sec'); ylabel('Amplitude(mv)'); title('RR Interval'); % Plot of RR-Interval
%% Autocorrelation for the time interval obtained
m = 100;                                        %sample shift
[X,R] = corrmtx(RR_Interval,m,'modified');      %autocorrelation of RR interval
[A,B] = ndgrid(1:m+1);                          %2D grid
figure
plot3(A,B,real(R))                              %3D plot of autocorrelation of RR interval
title('Re(R)') 
%% Histogram of RR interval
figure
histogram(RR_Interval);xlabel('samples'); ylabel('RR interval');title('Histogram of RR Interval');                          %Histogram plot
%% Power density Spectrum
t = 0:1/Fs:1-1/Fs;

N = length(RR_Interval);
xdft = fft(RR_Interval);                    %fourier transform 
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;           %magnitude square of fourier transform 
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(RR_Interval):Fs/2;       %frequency
figure
plot(freq,10*log10(psdx))                   %Power density Spectrum Plot
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
clear N
%% FIR filter to remove frequences exceeding 50 Hz

Pass_edge_freq =50; %pass edge freq
Trans_band =10; %transition band
N= round((3.3*Fs)/Trans_band);%N according to hamming window 
Mid_point= floor(N/2); % finding the midpoint
fc=(Pass_edge_freq+(Trans_band/Fs))/Fs; 
Hd(1)= 2*fc;
Wd(1)= 1; %according to Hamming window
H(1) = Hd(1) * Wd(1);
for x=2:N
    
Hd(x)=(2*fc*sin(2*pi*fc*(x-1)))/(2*fc*pi*(x-1)); 
Wd(x)=0.54 + (0.46*cos(((x-1)*2*pi)/N));

H(x)=Hd(x)*Wd(x);

end

figure;
plot(H); title('filter'); 
New_data=filter(H',1,data); %filtered data
figure;
plot(data'); hold on; plot(New_data);
legend('Input Data','Filtered Data');xlabel('samples'); ylabel('Amplitude(mv)');
%% Ploting filtered data
figure;
plot(New_data); xlabel('samples'); ylabel('Amplitude(mv)'); title('10 mins filtered ECG data'); % plotting the data
%% No. of heartbeats present in the filtered data
[peakVal, locVal] = findpeaks(New_data,'MinPeakHeight',0.5,'MinPeakDistance',10);         % detecting the highest peak
No_of_heart_beats_in_filtered = length(peakVal);                                          % No. of heartbeats in 10mins
Y=['The filtered data has ',num2str(No_of_heart_beats_in_filtered),' heart beats'];
disp(Y)
%% RR-Interval of filtered data
for k = 2:length(locVal)
    B = locVal(k,1)- locVal(k-1,1);
    RR_Interval_filtered(k,1) = B;
end
RR_Interval_filtered = RR_Interval_filtered(2:end,:); % RR-Interval durations
e=interp1(locVal,RR_Interval_filtered,'spline'); 
figure;
plot(RR_Interval_filtered); xlabel('time in sec'); ylabel('Amplitude(mv)'); title('RR Interval of filtered data'); % Plot of RR-Interval
%% cross correlation of filtered signal and downloaded data
[cor,lag] = xcorr(New_data,data);
figure
plot(lag,cor);xlabel('lag');ylabel('cross correlation');title('cross correlation of filtered signal and downloaded data');
%% Comparing filtered results with original data 
figure
subplot(2,1,1)
plot(data); xlabel('samples'); ylabel('Amplitude(mv)'); title('10 mins ECG data'); % plotting the data
subplot(2,1,2) 
plot(New_data); xlabel('samples'); ylabel('Amplitude(mv)'); title('10 mins filtered ECG data');

figure
subplot(2,1,1)
plot(RR_Interval); xlabel('time in sec'); ylabel('Amplitude(mv)'); title('RR Interval');
subplot(2,1,2)
plot(RR_Interval_filtered); xlabel('time in sec'); ylabel('Amplitude(mv)'); title('RR Interval of filtered data'); 