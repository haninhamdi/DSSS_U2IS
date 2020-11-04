%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cette fonction implémente l'algorithme de synchronisation de    %
% phase qui permet d'estimer le déphasage dans le mode Data       % 
% aided                                                           %
% Y =[y_1,...,y_N] : Les symboles reçus (complexes)               %
% S =[s_1,...,s_N] : Les symboles envoyés (complexes)             %
% mu : Le pas adaptatif                                           %
% iterations : le nombre des aller-retour                         %
% On commence par beta_F(1)=0                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [estimated_phase] = DA(Y,S,mu,iterations)
    % Calculer la dimension de Y 
    Nb = length(Y);
    % Initialisation des deux vecteurs beta_F et beta_B à 0
    beta_F = zeros(1,Nb);
    beta_B = zeros(1,Nb);
    % Calculer beta_B et beta_F pour plusieurs allers-retours 
    for k = 1:iterations 
        % Calculer beta_F
        beta_F(1) = beta_B(1);
        for k1 = 2:Nb
            a = Y(k1)*conj(S(k1))*exp(-1i*beta_F(k1-1));
            beta_F(k1) = beta_F(k1-1) + mu*imag(a) ; 
        end
        % Calcul de beta_B 
        beta_B(Nb) = beta_F(Nb);
        for k2 = 2:Nb
            b = Y(Nb-k2+1)*conj(S(Nb-k2+1))*exp(-1i*beta_B(Nb-k2+2));
            beta_B(Nb-k2+1) = beta_B(Nb-k2+2) + mu*imag(b);
        end
    end
    % Calcul de l'estimateur du vecteur déphasage 
    estimated_phase = zeros(1,Nb);
    estimated_phase(1) = beta_B(1);
    estimated_phase(Nb) = beta_F(Nb);
    for k = 2:Nb-1
        estimated_phase(k) = (beta_F(k)+beta_B(k))/2;
    end
    % Limiter estimated_phase à -pi/4 et pi/4
    for i = 2:Nb
        if estimated_phase(i)>pi/4 
            %estimated_phase(i) = pi/4;
            estimated_phase(i) = estimated_phase(i)-estimated_phase(i-1)/4 ;
        elseif estimated_phase(i)<-pi/4
            %estimated_phase(i) = -pi/4;
            estimated_phase(i) = estimated_phase(i)+estimated_phase(i-1)/4;
        end
    end
end
