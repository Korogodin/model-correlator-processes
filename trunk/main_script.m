clc
clear
close all

tau_chip = 1e-3 / 511;
Tmod = 8*tau_chip;

Fd = 1*44.2e6;
Td = 1/Fd;
L = fix(Tmod / Td);
t = (1:L)*Td;

dfi  = 0e6;
f0 = 14e6 + dfi;
fp = 10e6 + dfi;
phi = pi/3;
Gdk = ((1:L)>(L/8)) - ((1:L)>(2*L/8)) + ((1:L)>(3*L/8)) - ((1:L)>(4*L/8)) + ((1:L)>(5*L/8)) - 0.5;
S = Gdk .* cos(2*pi*f0*t + phi);
pOkno = nan(1,L);


RealFilter1 = 1;
if RealFilter1 == 0
    H = 16e6;
    tau_shift = 0*tau_chip/6; % BUG: Коэффициент неправильный ниже
    Okno = ones(1,L);
    pOkno = zeros(1,L);
    if H>0
        for l = 1:L
            Okno(l) = ((fp - H/2)/Fd*L < l)&&((fp + H/2)/Fd*L > l) || ((Fd-fp - H/2)/Fd*L < l)&&((Fd-fp + H/2)/Fd*L > l);
            pOkno(l) = (l <= L/2)*(-1)*l*(Fd/L)*2*pi*tau_shift + (l > L/2)*(+1)*(L-l)*(Fd/L)*2*pi*tau_shift;
        end
    end

    j = sqrt(-1);
    y_Filt = real(ifft( ( (Okno.*exp(j*pOkno))).*fft(S) ));
elseif RealFilter1 == 1
    load Hd1.mat
    y_Filt = filter(Hd1, S);
end

hF = 1;
hF = figure(hF) + 1;
plot(1e6*t, S, 1e6*t, y_Filt)
xlabel('t, \mu')
ylabel('S, y_{Filt}')

hF = figure(hF) + 1;
plot((0:(Fd/L):(Fd*(1-1/L)))/1e6,abs(fft(y_Filt)))
xlabel('f, MHz')
ylabel('|fft(y_{Filt})|')
xlim([0 (Fd*(1-1/L))/1e6])

hF = figure(hF) + 1;
plot((0:(Fd/L):(Fd*(1-1/L)))/1e6,angle(fft(y_Filt)), (0:(Fd/L):(Fd*(1-1/L)))/1e6, pOkno)
xlabel('f, MHz')
ylabel('angle(y_{Filt}), pOkno')
xlim([0 (Fd*(1-1/L))/1e6])

q = y_Filt .* sin(2*pi*f0*t);
i = y_Filt .* cos(2*pi*f0*t);

hF = figure(hF) + 1;
plot(i, q)
xlabel('i')
ylabel('q')
maxx = max(abs(i));
if max(abs(q))>maxx
    maxx = max(abs(q));
end
xlim([-maxx*1.2 +maxx*1.2])
ylim([-maxx*1.2 +maxx*1.2])

RealFilter2 = 1;
if RealFilter2 == 0
    LPH = 5e6;
    Okno_LPH = ones(1,L);
    pOkno_LPH = zeros(1,L);
    tau_shift_LPH = 0*tau_chip/15;
    for l = 1:L
        Okno_LPH(l) = ((0 - 0/2)/Fd*L < l)&&((0 + LPH)/Fd*L >= l) || ((Fd- 0 - LPH)/Fd*L < l)&&((Fd- 0 + 0)/Fd*L >= l);
        pOkno_LPH(l) = (l <= L/2)*(-1)*l*(Fd/L)*2*pi*tau_shift_LPH + (l > L/2)*(+1)*(L-l)*(Fd/L)*2*pi*tau_shift_LPH;
    end
    i_Filt = real(ifft(Okno_LPH.*exp(j*pOkno_LPH).*fft(i)));
    q_Filt = real(ifft(Okno_LPH.*exp(j*pOkno_LPH).*fft(q)));
elseif RealFilter2 == 1
    load Hd2.mat
    i_Filt = filter(Hd2, i);
    q_Filt = filter(Hd2, q);
end

hF = figure(hF) + 1;
if RealFilter2 == 0
    plot((0:(Fd/L):(Fd*(1-1/L)))/1e6, 1.2*max(abs(fft(i)))*Okno_LPH, (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(i)), (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(q)))
else
    plot((0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(i)), (0:(Fd/L):(Fd*(1-1/L)))/1e6, abs(fft(q)))
end
xlabel('f, MHz')
ylabel('Okno_{LPH}, fft i, fft q')
xlim([0 (Fd*(1-1/L))/1e6])

hF = figure(hF) + 1;
plot(i_Filt, q_Filt)
xlabel('i_{Filt}')
ylabel('q_{Filt}')
maxx = max(abs(i_Filt));
if max(abs(q_Filt))>maxx
    maxx = max(abs(q_Filt));
end
xlim([-maxx*1.2 +maxx*1.2])
ylim([-maxx*1.2 +maxx*1.2])

hF = figure(hF) + 1;
plot(t, rad2deg(atan2(q_Filt,i_Filt))/360*187, t, -phi/(2*pi)*187 - (Gdk + 0.5 - 1) * pi/(2*pi)*187)
xlabel('t, \mu')
ylabel('atan Filt, mm')
