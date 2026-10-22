t = scope108{:,1} + 0.1404;
u = scope108{:,2};
y = scope108{:,4};

close all
clc
plot(t,[u y])


K = mean(y)/mean(u);

%Mr = (max(y)-min(y))/(max(u)-min(u));

Mr = (2.88075 - 0.37)/(2.2875-0.475);
Mr_z = Mr/K;
roots([4 0 -4 0 1/Mr_z^2]);


T = 0.0017;
F = 1/T;
z = sqrt((1-sqrt(1-1/Mr_z^2))/2)    

%z = 0.43

wr = 2*pi/T;
wn = wr/sqrt(1-2*z^2)

%ph = rad2deg(-0.00015*wr/2)

Tzero = sqrt(Mr^2 * 4 * z^2 * (1 - z^2) - 1) / wr

A = [0, 1; -wn^2, -2*z*wn];
B = [0; K*wn^2];
C = [1 Tzero];
D = 0;

Hss = ss(A,B,C,D);
[ysim,~,~] = lsim(Hss,u,t,[y(1), (y(2)-y(1))/(t(2)-t(1))]);

figure
plot(t,[y ysim])

norm(y-ysim)/norm(y-mean(y))

%bode(Hss)

%c fara zero
Te = T/10;
data_wu_id = iddata(ysim,u,Te);
data_wu_vd = iddata(y,u,Te);

arx_wu = arx(data_wu_id,[2 3 1]);
armax_wu = pem(data_wu_vd,armax(data_wu_id,[3 3 3 0]));

figure; resid(data_wu_vd,armax_wu);
figure; compare(data_wu_vd,armax_wu);

iv_wu = pem(data_wu_vd,iv4(data_wu_id,[3 3 0]));

figure; resid(data_wu_vd,iv_wu);
figure; compare(data_wu_vd,iv_wu);