clear;
clc;
% Setting variables
fs = 44100; % Sampling rate in Hz
bits = 24; % Bits per sample
channels = 1; % Number of channels: 1

record_sample = audiorecorder(fs, bits, channels);

% Collect a sample of your speech with a microphone and plot the signal data:
% Record voice for 2 seconds
record_sample = audiorecorder;
disp('Start speaking.')
recordblocking(record_sample, 2);
disp('End of Recording.');

% Play back the recording
play(record_sample);

% Store data in double-precision array
recorded_sample = getaudiodata(record_sample);


% Converting audiorecorder object to wav
 audiowrite('sample01.wav',recorded_sample, fs,'BitsPerSample', bits);


[y,Fs] = audioread('sample01.wav');

%plot the input signal
t = linspace(0,length(y)/Fs,length(y));
figure(1);
plot(t,y);

% Number of samples
N = length(y)

% Generate plot of the frequency spectrum
f = fs/N .*(0:N-1); 

% Calculate each frequency component
Y=fft(y,N);
Y=abs(Y(1:N))./(N/2);

%plot the fft
figure(2);
plot(f,Y);

%import filters
Hd_m = filter_m;
Hd_f = filter_f;
%passing through the filter

filtered_wave_m = filter(Hd_m,y);
filtered_wave_f = filter(Hd_f,y);

G1 = abs(fft(filtered_wave_m,N))
G2 = abs(fft(filtered_wave_f,N));

%plotting filtered signals' fft
f1 = linspace(0,fs,16000);
figure(5);
subplot(2,1,1)
plot(f1(1:N*0.5),G1(1:N*0.5));
subplot(2,2,2)
plot(f1(1:N*0.5),G2(1:N*0.5));

%calculating rms
sum1 = 0;
for i=1:1:length(G1)
    sum1=sum1+G1(i)*G1(i);
end
disp sum1
R1 = sqrt(sum1/length(G1));

sum2 = 0;
for i=1:1:length(G2)
    sum2 = sum2 + G2(i)*G2(i)
end
R2 = sqrt(sum2/length(G2));

if(R1>R2);
    disp('THE VOICE IS OF A MALE');
else
    disp('THE VOICE IS OF A FEMALE');
    
end;
