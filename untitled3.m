%% 1.Задаём ФМ сигнал
clear;
clc;
barker_code=[ +1; +1 ;+1 ;0 ;0 ;+1 ;0];
bpskModulator = comm.BPSKModulator;
bpskModulator.PhaseOffset = 0;
modData = bpskModulator(barker_code);
tau = -6:1:6;
C_tau_barker = xcorr(modData);

%% 2.Строим АКФ полученного сигнала.
figure(1)
plot(tau, C_tau_barker, 'b', 'LineWidth',2);

%% 3.Строим тело неопределённости
lvl_noise_db = -40; %%Задаем уровень шума в дб
lvl_noise = db2pow(lvl_noise_db);
noise = @() wgn(1, 7, lvl_noise_db);
t = linspace(0, 1, length(noise()));
f = 1;
T = 1/f;
dt = t(2) - t(1);
Omega = -6*pi/T:1:6*pi/T;
signal_ref=reshape(modData,[1,7]);
C = zeros(length(signal_ref), length(Omega));
for i = 1:length(signal_ref)
for j = 1:length(Omega)
ind = (i - round(length(signal_ref)/2));
C(i, j) = sum(signal_ref.*circshift(signal_ref, ind).*exp(-1j*t*Omega(j)));
end
end
mesh(abs(C))

%% Рассматриваем сечение при при нулевом сдвиге по времени.
ylim([4,4.04])

%% Рассматриваем сечение при нулевом сдвиге по частоте.
xlim([20,20.2])

%% 4.Реализовываем обнаружитель на основе коррелятора с опорным сигналом и на основе коррелятора с принимаемым сигналом.
%Создаем корреляторы

NumOfExp = 1e7;
th_1 = zeros(1, NumOfExp);
th_2 = zeros(1, NumOfExp);
for i = 1:NumOfExp
signal = 0*signal_ref + noise();
th_1(i) = sum(signal.*signal_ref); % Выход коррелятора с опорным сигналом
th_2(i) = sum(signal.*signal); % Коррелятор с приёмным сигналом
end

%% 5. Вычисляем кривую оценки вероятности ложной тревоги до вероятности 10–6. Порог для вероятности ложной тревоги
[F1, th_val1] = ecdf(th_1);
[F2, th_val2] = ecdf(th_2);
ind1 = find((F1(end) - F1 < 1e-6));
ind2 = find((F2(end) - F2 < 1e-6));
thr_Plt1 = th_val1(ind1(1));
thr_Plt2 = th_val2(ind1(2)); %% Их должно быть 6
thr_Plt3 = th_val2(ind1(3));
thr_Plt4 = th_val2(ind1(4));
thr_Plt5 = th_val2(ind1(5));
thr_Plt6 = th_val2(ind1(6));

%%
%get_power = @(sig) sum(abs(sig(:)).^2)/numel(sig);
power_ref = sum(abs(signal_ref(:)).^2)/numel(signal_ref);
SNR_db = (-40:1:0);
SNR = db2pow(SNR_db);
cnt1 = zeros(size(SNR));
cnt2 = zeros(size(SNR));
coef = lvl_noise/power_ref;
NumOfExp = NumOfExp/10;
for i = 1:numel(SNR)
for j = 1:NumOfExp
signal = sqrt(SNR(i)*coef)*signal_ref + noise();
th_1 = sum(signal.*signal_ref);
th_2 = sum(signal.*signal);
if th_1 > thr_Plt1 % Проверяем если выход коррелятора больше порогового значения
cnt1(i) = cnt1(i) + 1;
end
if th_2 > thr_Plt2
cnt2(i) = cnt2(i) + 1;
end
end
end
p_1 = cnt1/NumOfExp;
p_2 = cnt2/NumOfExp;
h_2=SNR_db+pow2db(1/(t(2)-t(1)));

figure(2)
plot(SNR,p_1) % Кривая оценки вероятности правильного обнаружения ( зависимость от SNR) для коррелятора с опорным сигнало
figure(3)
plot(h_2,p_1) % Кривая оценки вероятности правильного обнаружения ( зависимость от h^2) для коррелятора с опорным сигналом
figure(4)
plot(SNR,p_2) % Кривая оценки вероятности правильного обнаружения ( зависимость от SNR) для коррелятора с приёмным сигналом
figure(5)
plot(h_2,p_2) % Кривая оценки вероятности правильного обнаружения ( зависимость от h^2) для коррелятора с приёмным сигналом

%% Следовательно для кривой оценки вероятности ложной тревоги:
p_l1=1-p_1;
p_l2=1-p_2;
figure(6)
plot(SNR,p_l1) %Кривая оценки вероятности ложной тревоги ( зависимость от SNR) для коррелятора с опорным сигналом
figure(7)
plot(h_2,p_l1) %Кривая оценки вероятности ложной тревоги ( зависимость от h^2) для коррелятора с опорным сигналом
figure(8)
plot(SNR,p_l2) %Кривая оценки вероятности ложной тревоги ( зависимость от SNR) для коррелятора с приёмным сигналом
figure(9)
plot(h_2,p_l2) %Кривая оценки вероятности ложной тревоги ( зависимость от h^2) для коррелятора с приёмным сигналом
clear;
clc;
%% 1. Построение тела неопределенности
B = [1 1 1 -1 -1 1 -1 -1 1 -1]; % код Баркера
mseq = [1 -1 1 1 -1 1 -1 -1 -1 -1 1 -1 1 1 -1 -1 1 1 1 -1]; % m-последовательность
code = kron(B, mseq); % сложный фазоманипулируемый сигнал
f_shift = 2; % сдвиг по частоте
t_body = 0:0.1:length(code)/10-0.1; % время
f_body = -5:0.1:5; % частота
[u, v] = meshgrid(t_body, f_body);
body = zeros(length(f_body), length(t_body)); % тело неопределенности
for i = 1:length(t_body)
t = t_body(i);
code_shifted = code.*exp(1j*2*pi*f_shift*t); % сдвиг по частоте
for j = 1:length(f_body)
f = f_body(j);
signal = code_shifted.*exp(1j*2*pi*f*t_body);
body(j,i) = abs(sum(signal))^2; % вычисление тела неопределенности
end
end
figure;
surf(u,v,body);
xlabel('Время, с');
ylabel('Частота, Гц');
zlabel('Тело неопределенности');
%% 2. Построение кривой ложной тревоги
BW = 10; % пропускная способность
Fs = 4*BW; % частота дискретизации
t = 0:1/Fs:1-1/Fs; % время
s_LFM = exp(1j*pi*BW*(t-1/(2*BW)).^2); % ЛЧМ-сигнал
s_HFM = exp(1j*pi*BW*(t-1/(2*BW)).^2).*exp(1j*pi*BW*t); % ГЧМ-сигнал
s_noise = randn(1,length(t)); % белый шум
t_seq = -5:5; % t-последовательность
s_seq = (t_seq>=0); % формирование t-последовательности
PRF = 1000; % импульсно-периодическая частота
s_repeating = repmat(s_LFM,1,round(PRF/Fs)); % повторение сигнала
thresholds = 0:0.01:1; % порог
p_false_alarm = zeros(1,length(thresholds)); % вероятность ложной тревоги
for i = 1:length(thresholds)
threshold = thresholds(i); % текущий порог
detections = zeros(1,1000); % массив хранения обнаружений
for j = 1:1000
s_noisy = s_repeating + threshold*s_noise(1:length(s_repeating)); % смесь сигнала и шума
xcorr_LFM = abs(conv(s_LFM,s_noisy,'same')); % вычисление кросс-корреляции
xcorr_HFM = abs(conv(s_HFM,s_noisy,'same')); % вычисление кросс-корреляции
xcorr_seq = abs(conv(s_seq,s_noisy,'same')); % вычисление кросс-корреляции
detection = max([xcorr_LFM xcorr_HFM xcorr_seq]); % максимальное значение кросс-корреляции
detection_norm = detection/sqrt(sum(abs(s_noisy).^2)*length(s_LFM)); % нормирование
if detection_norm>threshold % обнаружение сигнала
detections(j) = 1;
end
end
p_false_alarm(i) = 1-mean(detections); % вероятность ложной тревоги
end
figure;
semilogy(thresholds, p_false_alarm, 'LineWidth', 2);
xlabel('Порог');
ylabel('Вероятность ложной тревоги');
title('Кривая ложной тревоги');

%% 3. Построение кривой правильного обнаружения
SNR_dB = -20:2:20; % отношение сигнал-шум в дБ
SNR = 10.^(SNR_dB/10); % отношение сигнал-шум
threshold = 0.3; % значение порога из кривой ложной тревоги
p_detection = zeros(1,length(SNR)); % вероятность правильного обнаружения
for i = 1:length(SNR)
snr = SNR(i); % текущее значение отношения сигнал-шум
noise = sqrt(1/snr).*1;length(s_LFM); % генерация шума
s_noisy = s_LFM + noise; % смесь сигнала и шума
xcorr = abs(conv(s_LFM,s_noisy,'same')); % вычисление кросс-корреляции
detection = max(xcorr); % максимальное значение кросс-корреляции
detection_norm = detection/sqrt(sum(abs(s_noisy).^2)*length(s_LFM)); % нормирование
if detection_norm>threshold % обнаружение сигнала
p_detection(i) = 1;
end
end
figure;
plot(SNR_dB, p_detection, 'LineWidth', 2);
xlabel('Отношение сигнал-шум, дБ');
ylabel('Вероятность правильного обнаружения');
title('Кривая правильного обнаружения');
clc;
clear;

%% Код для построения тела неопределенности:
%Задаем параметры
f = 2.1e9; % частота
Ts = 1/10e6; % интервал дискретизации
N = 63; % длина Баркер-последовательности
mseq = [1 -1 -1 1 -1 1 1 1 -1 1 -1]; % m-последовательность

%% Создаем Баркер-последовательность
barker = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

% Строим сигнал
signal = kron(barker, mseq);

% Моделирование сдвига по частоте
shift = 5000; % сдвиг на 5000 Hz
time = (0:length(signal)-1)*Ts;
freq_shift = exp(1i*2*pi*time*shift);
signal = signal.*freq_shift;

% Вычисляем тело неопределенности
Fmax = 50000; % максимальный сдвиг по частоте в Гц
dF = 100; % шаг сдвига по частоте в Гц
Faxis = -Fmax:dF:Fmax;
Tdmax = 2*N-1; % максимальное значение задержки в отсчетах
Td = 0:Tdmax;
[Td, F] = meshgrid(Td, Faxis);
tau = Td*Ts;
fd = F/f;
ossbd = zeros(length(Faxis), length(Td));
for k=1:length(Faxis)
for l=1:length(tau)
sum = 0;
for n=1:N
N = length(s1);
cross_corr = conv(s1, s2_noisy);
cross_corr = cross_corr(ceil(N/2):(end-floor(N/2)));
%sum = sum + signal(n+l-1)*exp(-1i*2*pi*(fd(k))*tau(l)*(n-1));
end
ossbd(k,l) = abs(sum)^2/N; % Дискретное преобразование Фурье
end
end

% Визуализируем тело неопределенности
figure;
surf(T, F, 10*log10(ossbd.'));
xlabel('Задержка, отсчеты');
ylabel('Частота, Гц');
zlabel('Отношение сигнал/шум, дБ');
title('Тело неопределенности');

%% Код для построения кривой ложной тревоги:
% Сигналы
Ts = 1/10e6;
t = (0:1:127)*Ts;
s1 = sin(2*pi*200e3*t);
s2 = sin(2*pi*1e6*t);

% Шум
num_iter = 100;
sigma = 0.1;
threshold = 0:0.01:1;
p_false_alarm = zeros(size(threshold));
for i=1:length(threshold)
false_alarm_count = 0;
for j=1:num_iter
noise = sigma*randn(size(s1));
cross_corr = conv(s1+noise, s2);
cross_corr = cross_corr((length(s2)+1)/2:end-(length(s2)-1)/2);
max_corr = max(abs(cross_corr));
if max_corr > threshold(i)*sigma*sqrt(length(s2))
false_alarm_count = false_alarm_count + 1;
end
end
p_false_alarm(i) = false_alarm_count/num_iter;
end

% Визуализируем кривую ложной тревоги
figure;
semilogy(threshold, p_false_alarm);
xlabel('Порог');
ylabel('Вероятность ложной тревоги');
title('Кривая ложной тревоги');

%% Код для построения кривой правильного обнаружения:
% Фиксируем порог из кривой ложной тревоги
threshold = 0.5;

% Сигналы
Ts = 1/10e6;
t = (0:1:127)*Ts;
s1 = sin(2*pi*200e3*t);
s2 = sin(2*pi*1e6*t);

% Сигнал-шум
snr_db = -10:2:20;
num_iter = 100;
p_detection = zeros(size(snr_db));
for i=1:length(snr_db)
detection_count = 0;
for j=1:num_iter
noise = wgn(size(s1), 1, snr_db(i));
s2_noisy = s2 + noise;
cross_corr = conv(s1, s2_noisy);
cross_corr = cross_corr((length(s2)+1)/2:end-(length(s2)-1)/2);
max_corr = max(abs(cross_corr));
if max_corr > threshold*sigma*sqrt(length(s2))
detection_count = detection_count + 1;
end
end
p_detection(i) = detection_count/num_iter;
end

% Визуализируем кривую правильного обнаружения
figure;
semilogy(snr_db, p_detection);
xlabel('Отношение сигнал/шум, дБ');
ylabel('Вероятность правильного обнаружения');
title('Кривая правильного обнаружения');