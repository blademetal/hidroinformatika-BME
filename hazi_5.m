clc; clf;
clear all;


% A kiindul�si adatok megad�sa -------------------------------------------------
Qbe = 0.3; % m3/s


x_krit = [0, 3500, 6000, 8500, 11000];
zmeder = [87.82, 85.61, 84.45, 84.03, 83.21];

x_koz = [1200, 2700, 3400, 4700, 5700, 6400, 7000, 8100, 9400, 11000];
x_koz_min = [0, 1200, 2700, 3400, 4700, 5700, 6400, 7000, 8100, 9400];
x_koz_vegig = [0, 1200, 2700, 3400, 4700, 5700, 6400, 7000, 8100, 9400, 11000];
zfelsz = [87.68, 86.77, 86.40, 85.75, 85.28, 85.15, 85.13, 84.87, 84.46, 84.00];

hidr_h = [0.03, 0.06, 0.18, 0.33, 0.36, 0.48, 0.51, 0.72, 0.87, 0.96, 1];
hidr_P = [0.385, 0.583, 0.989, 1.359, 1.835, 3.923, 4.301, 4.896, 5.319, 5.574, 5.702];
hidr_B = [0.381, 0.566, 0.889, 1.107, 1.576, 3.651, 4.02, 4.44, 4.739, 4.919, 5];
hidr_A = [0.009, 0.023, 0.111, 0.26, 0.3, 0.613, 0.73, 1.618, 2.306, 2.74, 2.957];


% 1. l�p�s ---------------------------------------------------------------------
% A mederes�sek meg�llap�t�sa --------------------------------------------------
S_medereses = [];
for i=length(zmeder):-1:2
  S_medereses(i-1) = (zmeder(i-1)-zmeder(i))/(x_krit(i)-x_krit(i-1));
end

% V�zszintek a mederben --------------------------------------------------------
h_vizszintek = [];
S_keresztmetszetekben = [];
for i=length(x_koz):-1:1
  for j=length(x_krit):-1:1
    if(x_koz(i) == x_krit(j))
      h_vizszintek(i) = zfelsz(i)-zmeder(j);
      S_keresztmetszetekben(i) = S_medereses(j-1);
    elseif((x_koz(i) < x_krit(j)) && (x_koz(i) > x_krit(j-1)))
      dh = (x_krit(j)-x_koz(i))*S_medereses(j-1);
      z_helyi = zmeder(j)+dh;
      h_vizszintek(i) = zfelsz(i)-z_helyi;
      S_keresztmetszetekben(i) = S_medereses(j-1);
    endif
  end
end

% 2. l�p�s ---------------------------------------------------------------------
% V�zt�k�r sz�less�ge ----------------------------------------------------------
B_lin = interp1(hidr_h, hidr_B, h_vizszintek);

% Nedeves�tett szelv�nyker�let -------------------------------------------------
P_lin = interp1(hidr_h, hidr_P, h_vizszintek);

% Nedves�tett szelv�nyter�let --------------------------------------------------
A_lin = interp1(hidr_h, hidr_A, h_vizszintek);

% Hidraulikus sug�r ------------------------------------------------------------
R = [];
for i=1:length(A_lin)
  R(i) = A_lin(i)/P_lin(i);
end

% Szelv�ny k�z�psebess�gek -----------------------------------------------------
v = [];
for i=1:length(A_lin)
  v(i) = Qbe/A_lin(i);
end

% 3. l�p�s ---------------------------------------------------------------------
% Gravit�ci�s gyorsul�s --------------------------------------------------------
g = 9.81;

% Froude-sz�m ------------------------------------------------------------------
Fr = [];
for i=1:length(v)
  Fr(i) = v(i)*sqrt(B_lin(i)/(g*A_lin(i)));
end
% A jobb olvashat�s�g �rdek�ben �rdemes k�l�nbontani a k�pleteket, m�g ha 
% dr�g�bb is az algoritmus futtat�sa

% A Manning-f�le m�rt simas�gi egy�tthat� --------------------------------------
k_mert = [];
for i=1:length(v)
  k_mert(i) = v(i)/(power(R(i),(2/3)) * power(S_keresztmetszetekben(i),0.5));
end


% Az energiavonal es�se (m�rt) -------------------------------------------------
S_energia = [];
for i=1:length(v)
  S_energia(i) = power(v(i),2)/(power(k_mert(i),2) * power(R(i),(4/3)));
end

% 4. l�p�s ---------------------------------------------------------------------

% A keresztmetszetek k�z�tti t�vols�gok meghat�roz�sa --------------------------
deltaX = [];
deltaX(1) = x_koz(1);
for i=2:length(x_koz)
  deltaX(i) = x_koz(i)-x_koz(i-1);
end

% A prediktor l�p�sben meghat�rozott hi-1* v�zszintek --------------------------
h_vizszintek_csillag = [];
for i=length(h_vizszintek):-1:1
  h_vizszintek_csillag(i) = h_vizszintek(i)-deltaX(i)*((S_keresztmetszetekben(i)-S_energia(i))/(1-power(Fr(i),2)));
end


% ITT ELK�L�N�L A FELADAT ------------------------------------------------------
% Miut�n megkaptuk a prediktor l�p�sben a hi-1* v�zszintek, ezekb�l visszasz�moljuk
% a i-1. pontban l�v� param�tereket


% V�zt�k�r sz�less�ge ----------------------------------------------------------
B_lin2 = interp1(hidr_h, hidr_B, h_vizszintek_csillag);

% Nedeves�tett szelv�nyker�let -------------------------------------------------
A_lin2 = interp1(hidr_h, hidr_A, h_vizszintek_csillag);

% Nedves�tett szelv�nyter�let --------------------------------------------------
P_lin2 = interp1(hidr_h, hidr_P, h_vizszintek_csillag);

% Hidraulikus sug�r ------------------------------------------------------------
R2 = [];
for i=1:length(A_lin2)
  R2(i) = A_lin2(i)/P_lin2(i);
end

% Szelv�ny k�z�psebess�g -------------------------------------------------------
v2 = [];
for i=1:length(A_lin2)
  v2(i) = Qbe/A_lin2(i);
end

% Froude-sz�m ------------------------------------------------------------------
Fr2 = [];
for i=1:length(v2)
  Fr2(i) = v2(i)*sqrt(B_lin2(i)/(g*A_lin2(i)));
end

  
% Iter�ci�k sorozat�val megkeresni a megfelel� medersimas�gi t�nyez�t, mindezt �gy,
% hogy a m�rt �s a modellezett felsz�ng�rb�k k�l�nbs�ge a n�gyzetes k�z�phib�t tekintve minim�lis legyen

k_lista = [];
RMSE_vegso = [];
count = 1;

for k=19:0.01:25   
    S_energia2 = [];
    for i=1:length(v2)
      S_energia2(i) = power(v2(i),2)/(power(k,2) * power(R2(i),(4/3)));
    end

    H_i_1 = [];
    for i=length(h_vizszintek):-1:1
      H_i_1(i) = h_vizszintek(i)-deltaX(i)*0.5*(((S_keresztmetszetekben(i)-S_energia(i))/(1-Fr(i))) + 
      ((S_keresztmetszetekben(i)-S_energia2(i))/(1-Fr2(i))));
    end
    
    RMSE_h = 0;
    for j=1:length(h_vizszintek_csillag)
      RMSE_h += power(h_vizszintek_csillag(j)-H_i_1(j),2);
    end
    RMSE_vegso(count) = sqrt(RMSE_h/length(h_vizszintek_csillag));
    k_lista(count) = k; 
    count++;  
end

% Megtal�lni a legide�lisabb k medersimas�gi t�nyez�t --------------------------
RMSE_model_ind = ismember(RMSE_vegso,min(RMSE_vegso));
RMSE_model_best = find(RMSE_model_ind);
k_model_best = k_lista(RMSE_model_best);

Legkisebb_negyezetes_kozehiba = min(RMSE_vegso);
Legkisebb_negyezetes_kozehiba
Kalibralt_medersimasagi_ertek = k_model_best;
Kalibralt_medersimasagi_ertek

% Medersimas�gi egy�tthat� kalibr�l�s�nak �br�zol�sa grafikonon, ahol a miinimum
% pontba h�zott v�zszintes �rint� jelzi a legmegfelel�bb, a rendelkez�sre �ll�
% adatokb�l kikalibr�lhat� medersimas�gi egy�tthat�t
figure(1);
scatter(k_lista, RMSE_vegso);
hold on;
plot([min(k_mert), max(k_mert)], [min(RMSE_vegso), min(RMSE_vegso)]);
hold on;
title({'A n�gyzetes k�z�phiba alakul�s','a simas�gi egy�tthat� f�ggv�ny�ben'}, 'fontsize', 15);
xlabel('Simas�gi egy�tthat� (k) ((m1/3)/s)', 'fontsize', 15);
ylabel('N�gyzetes k�z�phiba (-)', 'fontsize', 15);
legendText = legend ({'N�gyzetes k�z�phiba sz�zados l�pt�kkel m�rve','Legalacsonyabb k�z�phiba'});
legend(legendText, 'location', 'southeast');
set(legendText, 'fontsize', 12);
set(gca,'XTick',min(k_lista):0.5:max(k_lista), 'fontsize', 10);
set(gca,'xlim', [min(k_lista), max(k_lista)]);


% Miut�n siker�lt meghat�rozni a megfelel� simas�gi egy�tthat�t
% kisz�m�that�v� v�lik a felsz�ng�rbe i-1.-ik tagja is, minden keresztemetszetben

% Az energiavonal es�se a m�r kalibr�lt modelben -------------------------------
S_energia2 = [];
for i=1:length(v2)
  S_energia2(i) = power(v2(i),2)/(power(k_model_best,2) * power(R2(i),(4/3)));
end

% A feladat megold�s�nak korrektor l�p�se, itt a jav�tott hi-1.-ik v�z�ll�sokat
% kapjuk meg
H_i_1 = [];
for i=length(h_vizszintek):-1:1
  H_i_1(i) = h_vizszintek(i)-deltaX(i)*0.5*(((S_keresztmetszetekben(i)-S_energia(i))/(1-Fr(i))) + 
  ((S_keresztmetszetekben(i)-S_energia2(i))/(1-Fr2(i))));
end

% A teljes csatorn�ra kalibr�lt felsz�ng�rbe magass�gi pontjai -----------------
hullam_javitott_teljes = H_i_1;
hullam_javitott_teljes(length(hullam_javitott_teljes)+1) = h_vizszintek(length(h_vizszintek));


% A felsz�ng�rb�k egym�sra halmozott �br�zol�sa --------------------------------
figure(2);
plot(x_koz_min, H_i_1, 'or-');
hold on,
plot(x_koz, h_vizszintek, '*b-');
hold on;
plot(x_koz_vegig, hullam_javitott_teljes, '*g-');
title('A felsz�ng�rbe alakul�sa', 'fontsize', 15);
xlabel('A befoly�si keresztmetszett�l m�rt t�vols�g (fm)', 'fontsize', 15);
ylabel('Az adott keresztszelv�nyben m�rt v�zszint (m)', 'fontsize', 15);
legendText = legend ({'Az hi-1 v�zszintek','A hi v�zszintek', 'A teljes hosszon modellezett felsz�ng�rbe'});
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:1000:11000, 'fontsize', 10);
set(gca,'xlim', [0, 11000]);
set(gca,'YTick',0:0.2:1.5, 'fontsize', 8);
set(gca,'ylim', [0, 1.5]);


% A felsz�ng�rb�k egym�st�l vertik�lisan elk�l�n�tett �br�i --------------------
figure(3);
subplot(3, 1, 1, 'align');
plot(x_koz_min, H_i_1, 'or-');
title({'A felsz�ng�rbe alakul�sa', '�sszetett �bra'}, 'fontsize', 15);
legendText = legend ('Az hi-1 v�zszintek');
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:1000:11000, 'fontsize', 10);
set(gca,'xlim', [0, 11000]);
set(gca,'YTick',0:0.2:1.5, 'fontsize', 8);
set(gca,'ylim', [0, 1.5]);

subplot(3, 1, 2, 'align');
plot(x_koz, h_vizszintek, '*b-');
ylabel('Az adott keresztszelv�nyben m�rt v�zszint (m)', 'fontsize', 15);
legendText = legend ('A hi v�zszintek');
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:1000:11000, 'fontsize', 10);
set(gca,'xlim', [0, 11000]);
set(gca,'YTick',0:0.2:1.5, 'fontsize', 8);
set(gca,'ylim', [0, 1.5]);

subplot(3, 1, 3)
plot(x_koz_vegig, hullam_javitott_teljes, '*g-');
xlabel('A befoly�si keresztmetszett�l m�rt t�vols�g (fm)', 'fontsize', 15);
legendText = legend ('A teljes hosszon modellezett felsz�ng�rbe');
legend(legendText, 'location', 'northwest');
set(legendText, 'fontsize', 12);
set(gca,'XTick',0:1000:11000, 'fontsize', 10);
set(gca,'xlim', [0, 11000]);
set(gca,'YTick',0:0.2:1.5, 'fontsize', 8);
set(gca,'ylim', [0, 1.5]);











