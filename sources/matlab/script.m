% Loading data
load('wfmdata.mat');
wfm_array = {Wfm00,Wfm01,Wfm02,Wfm03,Wfm04,...
    Wfm05,Wfm06,Wfm07,Wfm08,Wfm09,...
    Wfm10,Wfm11,Wfm12,Wfm13,Wfm14,...
    Wfm15,Wfm16,Wfm17,Wfm18,Wfm19,...
    Wfm20,Wfm21,Wfm22,Wfm23,Wfm24,...
    Wfm25,Wfm26,Wfm27,Wfm28,Wfm29,...
    Wfm30,Wfm31,Wfm32,Wfm33,Wfm34,...
    Wfm35,Wfm36,Wfm37,Wfm38,Wfm39,...
    Wfm40,Wfm41,Wfm42,Wfm43,Wfm44,...
    Wfm45,Wfm46,Wfm47,Wfm48,Wfm49,Wfm50};

freq2dist = 3*10^8/2/(58.5*10^9*1*10^3);
asdf = @(x) 2/length(x)*(sum(x.^2)-ifft(fft(x).*conj(fft(x))));

% Random preconditioning test
rand_wfm = wfm_array{randi([1 51],1)};
% Params
Y = rand_wfm.Y(5:end);
N = length(Y);
dt = rand_wfm.Y(2);
fs = 1/dt;
T = 0:dt:dt*(N-1);

figure
plot(T*1000, Y)
xlabel('Время, мс')
ylabel('Амплитуда')

% Wavelet 
cwt(Y,fs)
title("")
xlabel("Время, мс")
ylabel("Частота, МГц")
cb = get(gca,'Colorbar');
ylabel(cb,'Амплитуда')

% Highpass
hpY = highpass(Y, 20000, fs);

figure
plot(T,hpY)
xlabel('Время, с')
ylabel('Амплитуда')

% Crop
subY = hpY(floor(N/3):2*floor(N/3)-1);
subT = T(1:floor(N/3))';
subN = length(subT);

figure
plot(subT*1000,subY);
xlabel('Время, мс')
ylabel('Амплитуда')


% Smoothing
smY = smoothdata(subY);

figure
plot(subT*1000, smY)
xlabel('Время, мс')
ylabel('Амплитуда')

% Envelope
envY = abs(hilbert(smY));

figure
plot(subT*1000,smY)
hold on
plot(subT*1000,envY)
hold off
xlabel('Время, мс')
ylabel('Амплитуда')

% Equalizer
eqY = smY./envY;

figure
plot(subT*1000, eqY)
xlabel('Время, мс')
ylabel('Амплитуда')

% Zero crossing rate
zc1 = 0;
zc2 = 0;
dzc = [];
for i = 2:subN
    if eqY(i)*eqY(i-1) < 0
        zc1 = zc2;
        zc2 = subT(i);
        dzc = [dzc, zc2-zc1];
    end
end
zcr = 0.5/median(dzc);

% Least Squares Fit
if eqY(1) > 0
    f0 = pi*(1-dzc(1)/dzc(2));
else
    f0 = -pi*dzc(1)/dzc(2);
end

ft = fittype('sin(2*pi*a*x+b)');
f = fit(subT,eqY,ft,...
    'Start',[zcr, 0],...
    'Lower',[zcr*0.9, -pi],...
    'Upper',[zcr*1.1, pi])

lsf_freq = f.a
lsf_dist = lsf_freq * freq2dist * 1000

% Fourier
fftN = 2.^nextpow2(subN);
[maxVal, maxInd] = max(abs(fft(eqY,fftN)));
fft_freq = maxInd / fftN * fs
fft_dist = fft_freq * freq2dist * 1000

% ASDF
asdfY = asdf(subY);
pasdfY = asdfY(1:floor(subN/2));
[pks, locs] = findpeaks(-pasdfY,"MinPeakWidth",5);
loc = locs(1);
testT = loc-50:loc+50;
testY = pasdfY(testT);
testFT = fittype('a*(x-x0)^2-b');
testFit = fit(testT',testY,testFT,...
    'Start',[1,pasdfY(loc),loc])

asdf_freq = 1/ (testFit.x0 * dt)
asdf_dist = asdf_freq * freq2dist * 1000

figure
plot(1:floor(subN/2),-pasdfY)

% Main loop
lsf_dists = zeros(1,51);
trf_dists = zeros(1,51);
fft_dists = zeros(1,51);
fft2_dists = zeros(1,51);
fft4_dists = zeros(1,51);
zcr_dists = zeros(1,51);
asdf_dists = zeros(1,51);
acr_dists = zeros(1,51);

for i = 1:51 
    % Params
    Y = wfm_array{i}.Y(5:end);
    N = length(Y);
    dt = wfm_array{i}.Y(2);
    fs = 1/dt;
    T = 0:dt:dt*(N-1);
    
    % Highpass
    hpY = highpass(Y, 20000, fs);
    
    % Crop
    subY = hpY(floor(N/3):floor(2*N/3)-1);
    subT = T(1:floor(N/3))';
    subN = length(subT);
    
    % Smoothing
    smY = smoothdata(subY);
    
    % Envelope
    envY = abs(hilbert(smY));
    
    % Equalizer
    eqY = smY./envY;
    
    % Zero crossing rate
    zc1 = 0;
    zc2 = 0;
    dzc = [];
    for j = 2:subN
        if eqY(j)*eqY(j-1) < 0
            zc1 = zc2;
            zc2 = subT(j);
            dzc = [dzc, zc2-zc1];
        end
    end
    zcr = 0.5/mean(dzc(2:end-1));
    zcr_dists(i) = zcr*freq2dist*1000;
    
    % Least Squares Levenberg-Marquardt Fit
    if eqY(1) > 0
        f0 = pi*(1-dzc(1)/dzc(2));
    else
        f0 = -pi*dzc(1)/dzc(2);
    end

    ft = fittype('sin(2*pi*a*x+b)');
    f = fit(subT,eqY,ft,...
        'StartPoint',[zcr, f0],...
        'Algorithm','Levenberg-Marquardt');

    lsf_freq = f.a;
    lsf_dist = lsf_freq * freq2dist * 1000;
    lsf_dists(i) = lsf_dist;
    
    trft = fittype('sin(2*pi*a*x+b)');
    trf = fit(subT,eqY,ft,...
        'StartPoint',[zcr, f0],...
        'Upper',[1.1*zcr, pi],...
        'Lower',[zcr, -pi],...
        'Algorithm','Trust-Region');
    
    trf_freq = trf.a;
    trf_dist = trf_freq * freq2dist * 1000;
    trf_dists(i) = trf_dist;

    % Fourier 
    fftN = 2.^nextpow2(subN);
    [maxVal, maxInd] = max(abs(fft(eqY,fftN)));
    fft_freq = maxInd / fftN * fs;
    fft_dist = fft_freq * freq2dist * 1000;
    fft_dists(i) = fft_dist;

    % Fourier 4x padding
    fftN = 4*2.^nextpow2(subN);
    [maxVal2, maxInd] = max(abs(fft(eqY,fftN)));
    fft2_freq = maxInd / fftN * fs;
    fft2_dist = fft2_freq * freq2dist * 1000;
    fft2_dists(i) = fft2_dist;

    % Fourier 16x padding
    fftN = 16*2.^nextpow2(subN);
    [maxVal4, maxInd] = max(abs(fft(eqY,fftN)));
    fft4_freq = maxInd / fftN * fs;
    fft4_dist = fft4_freq * freq2dist * 1000;
    fft4_dists(i) = fft4_dist;
    
    % Autocorr
    acrY = ifft(fft(eqY).*conj(fft(eqY)));
    [acrpks, locs] = findpeaks(acrY);
    acr_freq = 1 / (locs(1) * dt);
    acr_dist = acr_freq * freq2dist * 1000;
    acr_dists(i) = acr_dist;

    % ASDF
    asdfY = asdf(eqY);
    [pks, locs] = findpeaks(-asdfY);
    asdf_freq = 1 / (locs(1) * dt);
    asdf_dist = asdf_freq * freq2dist * 1000;
    asdf_dists(i) = asdf_dist;

end

% Graphs
L = 0:0.1:5;
figure
plot(L,fft_dists)
hold on
plot(L,fft2_dists,'-.')
plot(L,fft4_dists,':')
plot(L, asdf_dists)
plot(L, acr_dists)
plot(L,lsf_dists)
plot(L,trf_dists)
plot(L,zcr_dists)
hold off
legend({'FFT','FFT 4x','FFT 16x','ASDF','ACR','LMF','TRF','ZCR'},'Location','best')
grid on
% Linearity
fft_pf = fit(L',fft4_dists','poly1')

figure
hold on
plot(L, fft4_dists, 'o')
plot(L, polyval([fft_pf.p1, fft_pf.p2],L))
hold off


lsf_pf = fit(L',lsf_dists','poly1')

figure
hold on
plot(L, lsf_dists, 'o')
plot(L, polyval([lsf_pf.p1, lsf_pf.p2],L))

% Output
% Ground Truth
E = L + 74.5;

% Fourier comparison
figure
hold on
plot(L, fft_dists)
plot(L,fft2_dists)
plot(L,fft4_dists)
plot(L,E,'--')
hold off
legend({'БПФ', 'БПФ-x4','БПФ-x16','Эталон'},'Location','best')
xlabel('Смещение мишени, мм')
ylabel('Рассчитанная разность плеч, мм')
grid minor

% LSQ comparison
figure
hold on
plot(L,zcr_dists)
plot(L,lsf_dists)
plot(L,trf_dists)
plot(L,E,'--')
hold off
legend({'Частота переходов через нуль', ....
    'МНК: Левенберг-Марквардт',...
    'МНК: Доверительная область',...
    'Эталон'},'Location','best')
xlabel('Смещение мишени, мм')
ylabel('Рассчитанная разность плеч, мм')
grid minor

% ACR comparison
figure
hold on
plot(L,acr_dists,'-o')
plot(L,asdf_dists,'-x')
plot(L,E,'--')
hold off

legend({'Автокорреляция','ASDF','Эталон'},'Location','best')
xlabel('Смещение мишени, мм')
ylabel('Рассчитанная разность плеч, мм')
grid minor
