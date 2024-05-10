function  [out1, out2]  = returnFFT(t,yy,Fs)

t = transpose(t);
yy = transpose(yy);

% Fs=1000; %Sampling frequency. Any component of the FFT cannot exceed half of the sampling frequency
T=t(end);    %The total time. The resolution of the FFT will be 1/T

% figure;
% plot(t,yy) %Plotting the input signal.

dummy = fft(yy)/length(yy); %Take fft of the signal and divide by the length.
% figure;
% plot(dummy); %Plot the first version, note that the vector has complex elements.
% figure;
% plot(abs(dummy)); %Find the magnitude of the complex values and plot the magnitude array using absolute value. 
% 
fft_yy=fftshift(dummy); %Mirror the FFT 
fft_yy=abs(fft_yy);
% figure;
% plot(fft_yy);

freq_axis = (-1*Fs/2):(1/T):(Fs/2);%Construct the frequency array, from -Fs/2 to +Fs/2 with increments of 1/T
freq_axis = freq_axis(1,1:end-1); %Throw away the last component to match the size.


% freq_axis = linspace(-1*Fs/2, Fs/2, numel(fft_yy));
closest2zero = min(abs(freq_axis)); 


ind=find(freq_axis==closest2zero); %Find the index where freq is 0;
% figure;
% plot(freq_axis(ind:end),fft_yy(ind:end)*2) %Only plot the portion that is bigger than 0 and multiply that portion by 2 to match the magnitude.
% xlabel('Frequency')
% ylabel('Amplitude')

out1 = freq_axis(ind:end);
out2 = fft_yy(ind:end).*2;