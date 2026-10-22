%% incarcare date
t = scope108{:,1};
u = scope108{:,2};
y = scope108{:,3};
y0 = scope108{:,4};

plot(t,[u y y0+6])
xlabel('Time[sec]');
ylabel('Voltage[V]');
legend('u(t)','y(t)','y0(t)');
%% identificare neparametrica cu ajutorul rezonantei

K = mean(y)/mean(u)

Mr = (y(405)-y(410))/(u(403)-u(408));
Mr_z = Mr;

z = sqrt((1-sqrt(1-1/Mr_z^2))/2)

T = 0.0018;
F = 1/T; 


wr = 2*pi/T;
wn = wr/sqrt(1-2*z^2)

A = [0, 1; -wn^2, -2*z*wn];
B = [0; K*wn^2];
C = [1 0];
D = 0;

Hss = ss(A,B,C,D);
[ysim,~,~] = lsim(Hss,u,t,[y(1), (y(2)-y(1))/(t(2)-t(1))]);

figure
plot(t,[y ysim])

norm(y-ysim)/norm(y-mean(y))
%% bode + panta

w1 = pi/(t(78)-t(47));
M1 = (y(47)-y(78))/(u(46)-u(76));
phi1 = -w1*(t(47)-t(46));

w2 = pi/(t(99)-t(78));
M2 = (y(99)-y(78))/(u(98)-u(76));
phi2 = -w2*(t(99)-t(98));

w3 = pi/(t(410)-t(405));
M3 = (y(405)-y(410))/(u(403)-u(408));
phi3 = -w3*(t(405)-t(403));

w4 = pi/(t(414)-t(410));
M4 = (y(414)-y(410))/(u(412)-u(408));
phi4 = -w4*(t(414)-t(412));

w5 = pi/(t(613)-t(610));
M5 = (y(610)-y(613))/(u(608)-u(611));
phi5 = -w5*(t(610)-t(608));

w6 = pi/(t(616)-t(613));
M6 = (y(616)-y(613))/(u(614)-u(611));
phi6 = -w6*(t(616)-t(614));

w7 = pi/(t(872)-t(870));
M7 = (y(870)-y(872))/(u(868)-u(870));
phi7 = -w7*(t(870)-t(868));

w8 = pi/(t(874)-t(872));
M8 = (y(874)-y(872))/(u(872)-u(870));
phi8 = -w8*(t(874)-t(872));

phi7_m = mean([phi5,phi7]);
phi8_m = mean([phi6,phi8]);

figure
semilogx([w1 w2 w3 w4 w5 w6 w7 w8],20*log10([M1 M2 M3 M4 M5 M6 M7 M8]),'*r')
hold on
bode(Hss)

semilogx([w1 w2 w3 w4 w5 w6 w7 w8],rad2deg([phi1 phi2 phi3 phi4 phi5 phi6 phi7_m phi8_m]),'*g')

panta = 20*log10(M8/M5)*10*w5/w8

%% identificare parametrica cu ARMAX si IV + minimizarea erorii de predictie cu PEM

Te = t(2)-t(1)
data_validare = iddata(y,u,Te);

armax_wu = pem(data_validare,armax(data_validare,[3 3 3 0]));

figure; resid(data_validare,armax_wu);
figure; compare(data_validare,armax_wu);

iv_wu = pem(data_validare,iv4(data_validare,[2 3 0]));

figure; resid(data_validare,iv_wu);
figure; compare(data_validare,iv_wu);
%% pentru sistem fara zero cu ARMAX

[num,den] = tfdata(armax_wu,'v');

Hd = tf(num,den,t(2)-t(1))

[A,B,C,D] = tf2ss(num,den);

A_t = A';
B_t = B';
C_t = C';

x0 = @(n)(pinv(B_t(1))*(y(n)-D*u(n)));
x1 = @(n)(x0(n+1)-A_t(1,1)*x0(n)-C_t(1)*u(n));
x2 = @(n)(x1(n+1)-A_t(2,1)*x0(n)-C_t(2)*u(n));

y2c = dlsim(A_t,C_t,B_t,D,u,[x0(1),x1(1),x2(1)]);

plot(t,[y y2c])
%hold on
%plot(t(1),y2c(1),'*g',t(1),y(1),'ob')

norm(y - y2c)/norm(y-mean(y))

%% pentru sistem fara zero cu IV

[num,den] = tfdata(iv_wu,'v')

Hd = tf(num,den,t(2)-t(1))

[A,B,C,D] = tf2ss(num,den)

A_t = A';
B_t = B';
C_t = C';

x0 = @(n)(pinv(B_t(1))*(y(n)-D*u(n)));
x1 = @(n)(x0(n+1)-A_t(1,1)*x0(n)-C_t(1)*u(n));

y2c = dlsim(A_t,C_t,B_t,D,u,[x0(1),x1(1)]);

plot(t,[y y2c])

norm(y - y2c)/norm(y-mean(y))

%%
Te = t(2) - t(1)
data_validare_0 = iddata(y0,u,Te);

armax_wu_0 = pem(data_validare_0,armax(data_validare_0,[3 3 3 0]));

figure; resid(data_validare_0,armax_wu_0);
figure; compare(data_validare_0,armax_wu_0);

iv_wu_0 = pem(data_validare_0,iv4(data_validare_0,[3 3 0]));

figure; resid(data_validare_0,iv_wu_0);
figure; compare(data_validare_0,iv_wu_0);

%% sistem cu zero cu ARMAX

[num,den] = tfdata(armax_wu_0,'v')

Hd = tf(num,den,t(2)-t(1))

[A,B,C,D] = tf2ss(num,den)

A_t = A';
B_t = B';
C_t = C';

x0 = @(n)(pinv(B_t(1))*(y0(n)-D*u(n)));
x1 = @(n)(x0(n+1)-A_t(1,1)*x0(n)-C_t(1)*u(n));
x2 = @(n)(x1(n+1)-A_t(2,1)*x0(n)-C_t(2)*u(n));

y2c0 = dlsim(A_t,C_t,B_t,D,u,[x0(1),x1(1),x2(1)]);

plot(t,[y0 y2c0])
% hold on
% plot(t(1),y2c0(1),'*g',t(1),y0(1),'ob')

norm(y0 - y2c0)/norm(y0-mean(y0))
%% sistem cu zero cu IV
[num,den] = tfdata(iv_wu_0,'v')

Hd = tf(num,den,t(2)-t(1))

[A,B,C,D] = tf2ss(num,den)

A_t = A';
B_t = B';
C_t = C';

x0 = @(n)(pinv(B_t(1))*(y0(n)-D*u(n)));
x1 = @(n)(x0(n+1)-A_t(1,1)*x0(n)-C_t(1)*u(n));
x2 = @(n)(x1(n+1)-A_t(2,1)*x0(n)-C_t(2)*u(n));

y2c = dlsim(A_t,C_t,B_t,D,u,[x0(1),x1(1),x2(1)]);

plot(t,[y0 y2c])
% hold on
% plot(t(1),y2c(1),'*g',t(1),y(1),'ob')

norm(y0 - y2c)/norm(y0-mean(y0))