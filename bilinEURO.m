function [ bilineare ] = bilinEURO( curva, DISPmechanism, PLOT )
% bilinEURO Calcola un'approssimazione bilineare secondo EuroCodice 8.
%
%   Una qualsiasi curva è bilinearizzata fissando lo spostamento (e
%   relativo taglio alla base) relativo alla formazione del meccanismo
%   plastico. (Metodo EuroCodice 8).
%
%
% INPUT
%
% curva          vettore Nx2 di punti che definisce la curva da bilineare
%
% DISPmechanism  spostamento in corrispondenza della formazione completa
%                del meccanismo plastico
%
% 
% OUTPUT
%
% 'bilineare' ? un vettore 3x2 che contiene il punto (0,0), il punto
% Y(xginocchio, yginocchio) e il punto U(xultimo, yultimo)

% % ESEMPIO
% a = [0 0; 1 3; 2 5; 3 5.5];
% vincolo = 2.5;
% bilin=bilinEURO(a,vincolo,'plot');

%% Function

% aggiungi DISPmechanism alla curva, eliminando eventuali doppioni (se DISPmechanism era già definito nella curva)
temp = unique([curva(:,1); DISPmechanism]);
temp(:,2) = interp1(curva(:,1) , curva(:,2), temp(:,1) );

% aggiorna la definizione della curva originale
curva = temp;


%area sottesa dalla curva di partenza (fino a DISPmechanism)
E = trapz(curva(curva(:,1)<=DISPmechanism,1),curva(curva(:,1)<=DISPmechanism,2));

% snervamento della bilineare
forceY = interp1(curva(:,1) , curva(:,2), DISPmechanism);
deltaY = 2 * (DISPmechanism - E / forceY);


% curva bi-lineare
bilineare = [   0               0;
                deltaY          forceY;
                curva(end,1)  	forceY];
    
%% plot di controllo

if strcmp(PLOT,'plot')
    figure
    plot(curva(:,1),curva(:,2))
    hold on
    plot(bilineare(:,1),bilineare(:,2))
end

end

