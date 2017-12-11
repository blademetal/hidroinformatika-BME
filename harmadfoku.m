clc; clf;
clear all;


% Befoly�s adatai --------------------------------------------------------------
ido = [0; 3; 6; 9; 12; 15; 18; 21; 24; 27; 30; 33; 36; 48; 60];
arhullam = [0.40, 0.40, 0.44, 1.10, 1.65, 5.20, 2.40, 2.16, 1.60, 0.90, 0.96, 0.60, 0.32, 0.32, 0.48];

% A t�rfogatok �s a hozz�juk tartoz� magss�gok ---------------------------------
v = [0;	150;	3220;	11466;	32330;	86593;	101376];
z = [220; 221;	222;	223;	224;	225;	226];

% A be�rkez� �rhull�m �br�zol�sa ----------------------------------------------- 
figure(1);
plot(ido, arhullam, '-b');
title('A be�rkez� �rhull�m', 'fontsize', 15);
xlabel('Ido (h)', 'fontsize', 15);
ylabel('Vizhozam (m3/s)', 'fontsize', 15);
legendText = legend ("Az �rhull�m v�zhozama");
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca, 'XTick', [0:3:ido(length(ido))]);


% Harmadfok� polinomi�lis regresszi� k�plete y = a*x3 + b*x2 + c*x + d ---------
p3 = polyfit(z, v, 3);

function result = getValue(val, pol)
  result = pol(1)*val**3+pol(2)*val**2+pol(3)*val+pol(4);
  return
end

% A 227 mBf-hez tartoz� t�rfogat -----------------------------------------------
V = getValue(227, p3);
z(8) = 227;
v(8) = V;

% A 226.7 mBf-hez tartoz� t�rfogat ---------------------------------------------
V_interpolalt = interp1(z, v, 226.7)

% Az t�rfogatadatok megjelen�t�se ----------------------------------------------
figure(2);
plot(z,v,'b-*');
title({'V�zszintekhez rendelt t�rfogatok','harmadfok� polinomi�lis regresszi�val'}, 'fontsize', 15);
xlabel('V�zszintek (m)', 'fontsize', 15);
ylabel('T�rfogatok (m3)', 'fontsize', 15);
legendText1 = legend ({"T�roz� v�zszinj�nek alakul�sa", "A 226.7 mBf-hez tartoz� t�roz�t�rfogat"});
legend(legendText1, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'YTick',0:20000:v(length(v)), 'fontsize', 10);
set(gca,'ylim', [0, 150000]);
hold on;
plot(226.7, V_interpolalt, 'r*');



% Kifolyas----------------------------------------------------------------------
mu = 0.6; % az �teresz vizhozamt�nyez�je
Aa = 0.02; % az �teresz szelv�nyter�lete
g = 9.81; % neh�zs�gi gyorsul�s

Q0 = []; % a taroz�b�l kifoly� v�zhozam
H = []; % az �teresz szintje folott mert vizoszlop magassag
h = 0;
 
% A kifoly�s sz�m�t�sa ---------------------------------------------------------
H_max = 0;
for i=221.3:0.01:228
  H_max++;
end

for i=1:H_max
  Q0(i) = mu*Aa*sqrt(2*g*h);
  H(i) = h+221.300;
  h+=0.01;
end

V_kif_alap = 150+((221.2-221)*(3220-150))/(222-221);

function result = getVolume(h, z, v)
  for i=2:size(v)(1)
    if ((220+h) > z(i-1)) && (z(i) >= (220+h))
      x = 220+h;
      x1 = z(i-1);
      x2 = z(i);
      y1 = v(i-1);
      y2 = v(i);
      result = y1+(((x-x1)*(y2-y1))/(x2-x1));
      return
    end
  end
end

% A kifoly�s �br�zol�sa --------------------------------------------------------
figure(3);
plot(Q0, H, '-o');
title('Kifoly� v�zhozam a magass�g f�ggv�ny�ben', 'fontsize', 15);
xlabel('A kifolyo vizhozam (m3/s)', 'fontsize', 15);
ylabel('Magass�gv�ltoz�s (m)', 'fontsize', 15);
legendText = legend ("T�roz� v�zszinj�nek alakul�sa");
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:0.01:Q0(length(Q0)), 'fontsize', 10);
set(gca,'xlim', [0, Q0(length(Q0))+0.01]);
set(gca,'YTick',221.3:0.5:H(length(H)), 'fontsize', 10);
set(gca,'ylim', [221.3, 228]);

% Magass�gok sz�m�t�sa ---------------------------------------------------------

%h_np1 = h_n + (1/(V/h220)) * (Qi - mu * Aa * sqrt(2*g*h_n)) * dt;
mag = [0];
for i=1:size(ido)(1)
   V_tar = getVolume(mag(i)+6.7, z, v)
   mag(i+1) = mag(i)+(1/(V_tar/(mag(i)+6.7)))*(arhullam(i)-(mu * Aa * sqrt(2*g*mag(i))))*ido(i);
end

kezd_mag = 226.700;
for i=1:length(mag)
  mag(i)=mag(i)+kezd_mag;
end

% Az �rhull�m hat�s�ra bek�vetkez� v�zszint �br�zol�sa -------------------------
figure(4);
plot(ido', mag(2:length(mag)), 'b-*');
title('V�zszintv�ltoz�s a t�roz�ban a be�rkez� �rhull�m hat�s�ra', 'fontsize', 15);
xlabel('Idosor (h)', 'fontsize', 15);
ylabel('Magass�g (m)', 'fontsize', 15);
legendText = legend ("T�roz� v�zszinj�nek alakul�sa");
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:3:ido(length(ido)), 'fontsize', 10);
set(gca,'xlim', [0, ido(length(ido))]);
set(gca,'YTick',226.7:0.002:mag(length(mag)), 'fontsize', 10);
set(gca,'ylim', [226.7, mag(length(mag))]);
print -djpg valtozas_harmad.jpg


% Az �rhull�m �s az t�roz� v�zszintj�nek egy�ttes �br�zol�sa -------------------

figure(5);
subplot(1, 2, 1, 'align');
plot(ido, arhullam, '-b');
title({'Az �rhull�m �s az t�roz�','v�zszintj�nek egy�ttes �br�zol�sa'}, 'fontsize', 12);
xlabel('Ido (h)', 'fontsize', 12);
ylabel('Vizhozam (m3/s)', 'fontsize', 12);
legendText = legend ("Az �rhull�m v�zhozama");
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 10);
set(gca, 'XTick', [0:3:ido(length(ido))]);

subplot(2, 2, 2, 'align');
plot(ido', mag(2:length(mag)), 'b-*');
title('V�zszintv�ltoz�s a t�roz�ban a be�rkez� �rhull�m hat�s�ra', 'fontsize', 12);
xlabel('Idosor (h)', 'fontsize', 12);
ylabel('Magass�g (m)', 'fontsize', 12);
legendText = legend ("T�roz� v�zszinj�nek alakul�sa");
legend(legendText, 'location', 'southeast');
set(legendText, 'fontsize', 10);
set(gca,'XTick',0:3:ido(length(ido)), 'fontsize', 8);
set(gca,'xlim', [0, ido(length(ido))]);
set(gca,'YTick',226.7:0.002:mag(length(mag)), 'fontsize', 8);
set(gca,'ylim', [226.7, mag(length(mag))]);

subplot(2, 2, 4)
plot(Q0, H, '-0');
title('Kifoly� v�zhozam a magass�g f�ggv�ny�ben', 'fontsize', 12);
xlabel('A kifolyo vizhozam (m3/s)', 'fontsize', 12);
ylabel('Magass�gv�ltoz�s (m)', 'fontsize', 12);
legendText = legend ("T�roz� v�zszinj�nek alakul�sa");
legend(legendText, 'location', 'southeast');
set(legendText, 'fontsize', 10);
set(gca,'XTick',0:0.01:Q0(length(Q0)), 'fontsize', 8);
set(gca,'xlim', [0, Q0(length(Q0))+0.01]);
set(gca,'YTick',221.3:0.5:H(length(H)), 'fontsize', 8);
set(gca,'ylim', [221.3, 228]);








  