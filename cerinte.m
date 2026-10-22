%raspuns in frecventa nyquist sau bode + identificarea parametrica(autocorelatie+intercorelatie)

%cateva puncte z frecv inalte+joase+rezonanta de la intrare si iesire(maxime si minime)
%facem diferenta intre ele t2-t1 t3-t2 pt frecv joase
%w1 = pi/(t2 - t1) w2 = pi/(t3 - t2) %jumate de perioada deci nu e 2*pi
%aflam amplificarea M1,2 pt t12 si t23
%repetam procedeul si pt rezonanta si frecv inalte
%subplot(semilog(w1..w6),20*log(M1...M6))
%determinam faza cu valoarea medie de la modul
%identificarea parametrica arx, armax + minimizarea erorii de predictie
%iddata t = t(2)-t(1)
%Mss = ss(iddata,1:7)
%(num,den) = ss2tf(Mss.A,Mss.B,Mss.C,Mss.D)
%resid(ysim,Mss) autocorelatie
%compare(ysim,Mss)
%[A,F] = bode(num,den,1)

%%

%n4sid pt valori singulare
%scoatem ss
%validare+corectie
%[num_z,den_z] = ss2tf(Mss_A,Mss_B,Mss_C(2,:),Mss_D(2,:));
%y2c = lsim(num_z,den_z,u)
%plot(t,[y2,y2c]) pt discret
%[num_c,den_c] = d2c(num_z,den_z,dt,'zoh');(dt-perioada de esant)
%Hc = tf(num_c,den_c)
%k = den_c(1,3)/num_c(1,3)
%wn = sqrt(den_c(1,3));
%z = den_c(1,2)/wn/2
%Ts = 1/z/wn
%tf(k*[Ts 1]*wn^2,[1 2*z*wn wn^2]);
%y2c_1 = lsim(Hc,u,t)
%plot(t,[y2 y2c y2c_1]);
%fie in continuu sau discret