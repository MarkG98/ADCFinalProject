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
xlabel('Sample')
ylabel('Magnitude')
subplot(212)
stem(imag(y));
title('Imaginary')
xlabel('Sample')
ylabel('Magnitude')

%% Decode

% Plot the FFT of the recieved signal squared in order to derive
% frequency and phase offset
N = length(real(y));
frequencies_shifted = (linspace(-pi, pi-2/N*pi, N) + pi/N*mod(N,2));
figure;
a = fftshift(fft(y.^4));
plot(frequencies_shifted,abs(a));
xlabel('Frequency (Radians/Sample)')
ylabel('Magnitude')
title('FFT of Recieved Signal Raised to the Fourth')


% Calculate frequency and phase offset
[M,I] = max(abs(a));
foffset = frequencies_shifted(I)/4;
aoffset = (angle(a(I)) + pi)/4;

% Divide output by the exponential coefficient in order to get x * h which
% is good enough because the channel is fairly small compared to our
% message signal
times = 0:1:length(y) - 1;
expon = exp(j*(foffset*times + aoffset));
res = (y.'./expon);

% Headers
res1_h = res(1:length(header));
res2_h = res(1:length(header))*exp(j*(-pi/2));
res3_h = res(1:length(header))*exp(j*(-2*pi/2));
res4_h = res(1:length(header))*exp(j*(-3*pi/2));

% Collect header from each rotation
i = 1;
for m = 10:20:length(res1_h)
    
    h_vals(i) = sign(real(header(m))) + sign(imag(header(m)))*j;
    
    h1(i) = sign(real(res1_h(m))) + sign(imag(res1_h(m)))*j;
    h2(i) = sign(real(res2_h(m))) + sign(imag(res2_h(m)))*j;
    h3(i) = sign(real(res3_h(m))) + sign(imag(res3_h(m)))*j;
    h4(i) = sign(real(res4_h(m))) + sign(imag(res4_h(m)))*j;
    i = i + 1;
end

% Calculate header error
error = [length(find((h1 - h_vals) ~= 0)); length(find((h2 - h_vals) ~= 0)); length(find((h3 - h_vals) ~= 0)); length(find((h4 - h_vals) ~= 0))];

% Establish which rotation is correct based on the header comparison
[M,I] = min(error);
if (I == 1)
    res = res;
elseif (I == 2)
    res = res*exp(j*(-pi/2));
elseif (I == 3)
    res = res*exp(j*(-2*pi/2));
else
    res = res*exp(j*(-3*pi/2));
end

% Plot recieved bits
figure;
stem(real(res))
title('Real Bits')
xlabel('Bit Number')
ylabel('Magnitude')

figure;
stem(imag(res))
title('Imaginary Bits')
xlabel('Bit Number')
ylabel('Magnitude')

% Extract bits from transmitted and received messages
i = 1;
for m = 10:20:length(res)
    received_bits_R(i) = sign(real(res(m)));
    received_bits_I(i) = sign(imag(res(m)));
    i = i + 1;
end

received_bits = received_bits_R + received_bits_I*j;