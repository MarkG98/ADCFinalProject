% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');

% read data from the file
tmp = fread(f2, 'float32');

% close the file
fclose(f2);


% since the USRP stores the data in an interleaved fashion
% with real followed by imaginary samples 
% create a vector of half the length of the received values to store the
% data. Make every other sample the real part and the remaining samples the
% imaginary part
y = zeros(length(tmp)/2,1);
y = tmp(1:2:end)+j*tmp(2:2:end);

% Store original signal
orig = y;

% Cross correlate the recieved
[r, lags] = xcorr(orig,header);
[~,I] = max(r);
delay = lags(I);
y = orig(delay:delay+1100*20);

% to visualize, plot the real and imaginary parts separately
%return;
subplot(211)
stem(real(y));
title('Real')
subplot(212)
stem(imag(y));
title('Imaginary')

%% Decode

% Plot the FFT of the recieved signal squared in order to derive
% frequency and phase offset
N = length(real(y));
frequencies_shifted = (linspace(-pi, pi-2/N*pi, N) + pi/N*mod(N,2));
figure;
a = fftshift(fft(y.^2));
plot(frequencies_shifted,abs(a));

% Calculate frequency and phase offset
[M,I] = max(abs(a));
foffset = frequencies_shifted(I)/2;
aoffset = angle(a(I))/2;

% Divide output by the exponential coefficient in order to get x * h which
% is good enough because the channel is fairly small compared to our
% message signal
times = 0:1:length(y) - 1;
expon = exp(j*(foffset*times + aoffset));
res = y.'./expon;

% Checking if phase was corrected in the correct direction
% and changing sign if necessary
if (abs(min(real(r))) > abs(max(real(r))))
    res = -real(res);
end


% Plot recieved bits
figure;
stem(res)

% Extract bits from transmitted and received messages
i = 1;
for m = 10:20:length(res)
    received_bits(i) = sign(real(res(m)));
    i = i + 1;
end