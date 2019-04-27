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

indexes = find(abs(real(y))>0.0008);
y = y(indexes);

% to visualize, plot the real and imaginary parts separately
%return;
subplot(211)
stem(real(y));
title('Real')
subplot(212)
stem(imag(y));
title('Imaginary')

%% Our Code

N = length(y);
frequencies_shifted = (linspace(-pi, pi-2/N*pi, N) + pi/N*mod(N,2));
figure;
a = fftshift(fft(y.^2));
plot(frequencies_shifted,abs(a));

[M,I] = max(abs(a));
foffset = frequencies_shifted(I)/2;
aoffset = angle(a(I))/2;

times = 0:1:length(y) - 1;
expon = exp(j*(foffset*times + aoffset));
res = y.'./expon;

figure;
plot(real(res))

%i = 1;
%res1 = [];
%for m = 1:1:length(real(res))
    %if (abs(real(res(m))) > 3*10^-3)
        %res1(i) = real(res(m));
        %i = i + 1;
    %end
%end