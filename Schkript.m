
close all; clear; clc; format long

dat = readmatrix("aufgabe5.dat");
u = dat(:,3);
t = dat(:,2);

clearvars dat

%% Schritt 1: Sicherstellen, dass u(t) eine gerade Anzahl an Werten enthält

if mod(length(u),2)==0
else
u = u(1:end-1); % Kürzen der Zeitreihe um einen Wert
end

figure(1)
hold on
plot(t,u,'b.-')
title("GNSS Up-Reihe")
xlabel("Zeit [t]")
ylabel("u[t]")
hold off

%% Schritt 2: Bestimmen der Kenngrößen

I = length(u); % Länge der Zeitreihe
dt = t(2)-t(1);
T = I*dt; % Abtastweite im Frequenzbereich
fg = 1/(2*T); % Grenzfrequenz oder Nyquistfrequenz

%% Schritt 3: Erste Parameterschätzung

A = [ones(I,1) (0:dt:T-dt)' (cos(2*pi*(1:I)'/365.25)) (sin(2*pi*(1:I)'/365.25)) (cos(2*pi*(1:I)'/(365.25/2))) (sin(2*pi*(1:I)'/(365.25/2)))]; % Erstellen der Designmatrix

l = u;

[ xdach, v, sx ] = ausgleichung( A, l );

% [offset, trend, j1, j2, hj1, hj2] = xdach; % Zerlegen des Vektors b in die gesuchten Parameter

sigma2 = 1/(I-6)*sum((u'-A*xdach).^2); % Ermitteln der Varianz des Rauschens
s0 = sqrt(sigma2); % Standardabweichung der Gewichtseinheit

% Berechne Standardabweichungen
s = sqrt((b-A*x)'*(b-A*x)/(length(t)-p));
sigma = diag(inv(A'*A))*s^2;

% Ausgabe der Parameter
fprintf('Offset: %f \n', x(1));
fprintf('Trend: %f mm/a \n', x(2)*365.25);
fprintf('Kosinus-Term des jährlichen Signals: %f \n', x(3));
fprintf('Sinusterm des jährlichen Signals: %f \n', x(4));
fprintf('Kosinus-Term des halbjährlichen Signals: %f \n', x(5));
fprintf('Sinusterm des halbjährlichen Signals: %f \n', x(6));

% Ausgabe der Standardabweichungen der Parameter
fprintf('Standardabweichungen der Parameter: \n');
disp(sqrt(sigma));

% Ausgabe der Standardabweichung der Gewichtseinheit
fprintf('Standardabweichung der Gewichtseinheit s0 = qs02 : %f \n', s);

%% Schritt 4: Darstellen der Zeitreihe der Residuen

e = u - A*b; % Berechnen der Residuen

figure;
plot(1:I,e); % Diagramm (t, e(t))
xlabel('t'); ylabel('e(t)');
title('Diagramm der Residuen');

%% Schritt 5: Schätzen der Energiedichte P(f)

E = abs(fft(e,I)).^2; % Berechnen des Spektrums E

E1 = E(1:end/2+1); % Berechnen der einseitigen Energiedichte
E1(2:end-1) = 2*E1(2:end-1);

P = E1/I; % Berechnen der Energiedichte

f = 0:1/T:fg; % Erstellen der Frequenzachse

figure;
plot(f,P); % Diagramm (f, P(f))
xlabel('f'); ylabel('P(f)');
title('Diagramm der Energiedichte');


I=length(u);
%1
T=I;
df=1/T;
fg=1/(2*df);
%2
A=[ones(I,1),[1:I]'];
x=A\u';
s0=sqrt(sum((u'-A*x).^2)/(I-2));
%3
e=u'-A*x;
%4
E=fft(e);
P=abs(E).^2/(I^2);
P(1)=2*P(1);
P=P(1:round(I/2));
%5
f=[1:length(P)]'*df;
X=[log(f),log(P)];
b=X\ones(length(f),1);
k=b(1);
%6
c0=1;
U=eye(I);
h1=1;
for i=2:I
    h=(h1/(i-1))*(i-(k/2)-2);
    for j=1:(I-i+1)
        U(j,j+i-1)=h;
    end
    h1=h;
end
CPL=c0*U'*U;
%7 
A=[ones(I,1),[1:I]'];
Q=A'*CPL\A;
R=A'*CPL\u';
x2=Q\R;
s2=sqrt(diag(inv(Q)));
%8
C=corrcov(Q);
figure;
imagesc(C);
figure;
plot(u);
hold on
plot(A*x2);
