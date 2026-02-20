function [ bilineare ] = bilinATC( curva, vincoloK, PLOT )
% Calcola un'approssimazione bilineare di una qualsiasi curva con rigidezza
% iniziale fissata. (Metodo ATC 1996 "generalizzato").
%
%
% L'algoritmo calcola un'approssimazione della curva imponendo che questa 
% abbia una rigidezza iniziale fissata. Si fornisce lo spostamento di un
% punto (vincoloK) che viene usato per definire la rigidezza iniziale.
% L'approssimazione bi-lineare ? ottenuta eguagliando l'energia sottesa
% dalla curva bilineare a quella della curva originale.
%
%
% Il senso di definire la rigidezza come "secante che passa per un punto
% dato" ? quello di estendere i criteri ATC (basati su rigidezza iniziale
% tangente) ai casi in cui nell'analisi pushocer si considera il cracking nei  
% modelli costitutivi delle sezioni.
%
%
% INPUT
%
% curva     vettore Nx2 di punti che definisce la curva da bilineare
%
% vincoloK  scalare che rappresenta l'ascissa di un punto usato per
%           definire la rigidezza iniziale della curva bilineare.
%
% 
% OUTPUT
%
% 'bilineare' ? un vettore 3x2 che contiene il punto (0,0), il punto
% Y(xginocchio, yginocchio) e il punto U(xultimo, yultimo)

% ESEMPIO
% a=[0 0; 1 3; 2 5; 3 5.5];
% vincolo=a(2,1);
% bilin=bilinATC(a,vincolo,'plot');


% parametro che risolve i problemi di monotonia nell'interpolazione
control = 10;


% punto ultimo della curva
xultimo = curva(end,1);
yultimo = curva(end,2);

%area sottesa dalla curva di partenza
E = trapz(curva(:,1),curva(:,2)); 

% rigidezza iniziale (tangente o secante a seconda di come scelgo vincoloK)
try
    Kiniz = (interp1(curva(:,1),curva(:,2), vincoloK)) / vincoloK;
catch
    Kiniz = (interp1(curva(1:control:end,1),curva(1:control:end,2), vincoloK)) / vincoloK; % se il precedente rigo da errore vuol dire che la curva non ? monotona. Allora salto dei punti per renderla monotona
end
   
% spostamento al ginocchio (basato su ugual energia)
xginocchio = (2*E - yultimo*xultimo) / (Kiniz*xultimo - yultimo);

% forza al ginocchio
yginocchio = Kiniz * xginocchio;

% curva bi-lineare
bilineare = [   0         	0;
                xginocchio 	yginocchio;
                xultimo  	yultimo];
    
% plot di controllo
if strcmp(PLOT,'plot')
    figure
    plot(curva(:,1),curva(:,2))
    hold on
    plot(bilineare(:,1),bilineare(:,2))
end

end

