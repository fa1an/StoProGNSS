function [Amplituden, Frequenzen] = amplitudenspektrum(Zeitpunkte, Zeitreihe)
  % Trend abziehen
  ZR_detrend = detrend(Zeitreihe,'linear');
  
  % Zeige Zeitreihen
  %figure
    %plot(Zeitpunkte, Zeitreihe)
    %hold on
    %plot(Zeitpunkte, ZR_detrend)
    %hold off

  % Abtastfrequenz: 365.25 mal pro Jahr
  Abtastfrequenz = 365.25;
  % Frequenzen fuer die Amplituden berechnet werden
  Nyquistfrequenz = Abtastfrequenz/2;
  Frequenzintervall = Abtastfrequenz/length(Zeitpunkte);
  Frequenzen = 0:Frequenzintervall:Nyquistfrequenz;

  % Diskrete Fourier-Transformation mittels Fast-Fourier-Transformation
  ZR_fft = fft(ZR_detrend);

  % Amplitudenspektrum ergibt sich aus dem Betrag der komplexen Fourierkoeffizienten
  amp_spek = abs(ZR_fft);

  % Normalisierung der Amplituden
  Amplituden = 2 * amp_spek/length(Zeitpunkte);
  Amplituden(1,1) = Amplituden(1,1)/2;
  
  % Visualisierung
  % Einschraenken des Frequenzbereichs auf "alle 4 Jahre" bis "dreimal pro Jahr"
  % periodisches Signal, das alle 4 Jahre vorkommt: f = 0.25
  % periodisches Signal, das dreimal pro Jahr vorkommt: f = 3
%   ind_plt = find(Frequenzen>=0.25 & Frequenzen<=3);
%   Frequenzen_plt = Frequenzen(ind_plt);
%   Amplituden_plt = Amplituden(ind_plt);
    
    f = 
  figure
    plot(Frequenzen, )
    xlabel('Frequenz in pro Jahr')
    ylabel('Amplitude')

end
