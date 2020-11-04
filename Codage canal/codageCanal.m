%%%%%%%%%%%%%%%%%%%% PRE_U2IS : 2ème et 3ème semaine %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Objectif : Synchronisation de phase + codage canal %%%%%%%%%%%%
%************************************************************************%                                                           %
% mu : le pas adaptatif utilisé dans l'algorithme d'estimation de phase :%
%      il vaut mieux trouver un algorithme qui le fixe. On utilise       %
%      différents mu pour les deux modes NDA et DA.                      %
% iteration : nombre des allers retours dans l'algorithme d'estimation   %
% Seuil : un seuil utilisé pour le décodage des signaux reçus            %
% Nb : le nombre des bits envoyés                                        %
% ksi : le facteur qui reprèsente l'effet Doppler                        %
% RSB : le rapport signal sur bruit. Varie de 0 à 20 .                   %
% b : les bits réellement générés                                        %
% S : les bits codés avec une modulation BPSK : 1 et -1                  %
% S_complexe : l'enveloppe complexe de S                                 %
% n : bruit blanc additif gaussien : BBAG                                %
% dephasage : un déphasage ajouté au cours de la transmission, ce qu'on  %
%             cherche à estimer.                                         %
% Y_complexe : l'enveloppe complexe du signal reçu.                      %
% beta_estime_DA/NDA : le déphasage estimé dans les deux modes DA et NDA %
% Y_DA_coded, Y_NDA_coded : les bits décodés après une synchronisation de% 
%                           phase dans les deux modes DA et NDA.         %
% Y_coded : les bits décodés à partir du signal reçu                     %
% Yr_coded : les bits décodés dans le cas où aucun déphasage n'est       %
%            introduit                                                   %
% TEB_... : le taux d'erreur binaire dans le cas où : aucun déphasage    %
%           n'est introduit, aucune sychronisation n'est faite,          %
%           synchronisation de phase dans le mode DA et NDA.             %
% MSE_DA/NDA : l'erreur quadratique moyenne des deux estimateurs de phase%
%              pour les deux modes DA et NDA.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;
% Définir les paramètres de simulation
trellis = poly2trellis(3,{'1 + x^2','1 + x + x^2'}); 
Nb = 100000;
iterations = 3;
ksi = 0.01;
Seuil = 0;
mu_DA = 0.3 ;
sigma_phase = sqrt(0.005);  % A changer

for RSB = 0 : 40
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Générer les bits
    b = randi([0,1], 1, Nb);
    % Codage canal 
    dataEnc = convenc(b,trellis);
    % Coder les bits : modulation BPSK  
    for k1 = 1 : length(b)
        if b(k1) == 0
            S(k1) = -1;
        else
            S(k1) = 1;
        end
    end
    for k1 = 1 : length(dataEnc)
        if dataEnc(k1) == 0
            Senc(k1) = -1;
        else
            Senc(k1) = 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Canal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Introduire un BBAG : bruit blanc additif gaussien 
    N0 = 10^(-RSB/10);
    sigma_bruit = sqrt(N0);
    n = sigma_bruit * randn(1,length(S));
    nEnc = sigma_bruit * randn(1,length(Senc));
    % Introduire un déphasage causé par le canal, la phase suit une marche
    % Brownienne, appartient à [-pi/4, pi/4] et affectée pour l'effet
    % Doppler (ksi)
    
    bruit_phase = sigma_phase * randn(1,length(S)-1)*pi/180;
    dephasage = zeros(1,length(S));
    for k1 = 2 : length(dephasage)
        c = dephasage(k1-1) +ksi + bruit_phase(k1-1);
        if (c>=-pi/4)&&(c<=pi/4) dephasage(k1)=c;
        elseif c>pi/4 dephasage(k1) = c - dephasage(k1-1)/4;
        elseif c<-pi/4 dephasage(k1) = c + dephasage(k1-1)/4;
        end
        %elseif c>pi/4 dephasage(k1) = pi/4;
        %elseif c<-pi/4 dephasage(k1) = -pi/4;
        %end 
    end
    bruit_phaseEnc = sigma_phase * randn(1,length(Senc)-1)*pi/180;
    dephasageEnc = zeros(1,length(Senc));
    for k1 = 2 : length(dephasageEnc)
        c = dephasageEnc(k1-1) +ksi + bruit_phaseEnc(k1-1);
        if (c>=-pi/4)&&(c<=pi/4) dephasageEnc(k1)=c;
        elseif c>pi/4 dephasageEnc(k1) = c - dephasageEnc(k1-1)/4;
        elseif c<-pi/4 dephasageEnc(k1) = c + dephasageEnc(k1-1)/4;
        end
        %elseif c>pi/4 dephasage(k1) = pi/4;
        %elseif c<-pi/4 dephasage(k1) = -pi/4;
        %end 
    end
    
    erreurDePhase = exp(1i*dephasage);
    erreurDePhaseEnc = exp(1i*dephasageEnc);
    % Construire le signal reçu : y = s*exp(j*beta)+n
    Y_complexe = S.*erreurDePhase + n ;
    Y_complexeEnc = Senc.*erreurDePhaseEnc + nEnc ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Estimation de la phase : mode DA
    %mu_DA = 2*((ksi/sigma_bruit)^(2/3)); % A changer  
    beta_estime_DA = DA(Y_complexe,S,mu_DA,iterations);
    beta_estime_DAenc = DA(Y_complexeEnc,Senc,mu_DA,iterations);
    % Estimation de la phase : mode NDA
    % Ces conditions sont à changer pour optimiser l'algorithme
    beta_estime_NDA = NDA(Y_complexe,mu_DA,iterations,N0);
    beta_estime_NDAenc = NDA(Y_complexeEnc,mu_DA,iterations,N0);
    %%%%%%%%%%%%% Test de la convergence de l'estimateur
    %X = 0.02:0.02:0.02*length(dephasage);
    %figure 
    %plot(X,dephasage,'r-')
    %hold on 
    %plot(X,beta_estime_DA,'b-')
    %hold on 
    %plot(X,beta_estime_NDA,'g--')
    %grid on 
    %legend("déphasage réel","déphasage estimé en mode DA","déphasage estimé en mode NDA")
    %title('Le déphasage')
    %ylabel("déphasage (rad)")
    %xlabel("temps (\mus)")
    
    % Synchronisation à partir de beta estimé 
    % Cas DA 
    Y_DA = Y_complexe.*exp(-1i*beta_estime_DA/2.3);
    Y_DAenc = Y_complexeEnc.*exp(-1i*beta_estime_DAenc/2.4);
    % Cas NDA 
    Y_NDA = Y_complexe.*exp(-1i*beta_estime_NDA/2.4);
    Y_NDAenc = Y_complexeEnc.*exp(-1i*beta_estime_NDAenc/2.4);
    % Décodage du signal reçu, avec synchronisation 
    % Cas DA 
    Y_I_DA = real(Y_DA);
    for k1 = 1 : length(Y_I_DA)
        if Y_I_DA(k1) > Seuil
            Y_DA_coded(k1) = 1;
        else 
            Y_DA_coded(k1) = 0;
        end
    end
   Y_I_DAenc = real(Y_DAenc);
    for k1 = 1 : length(Y_I_DAenc)
        if Y_I_DAenc(k1) > Seuil
            Y_DA_codedEnc(k1) = 1;
        else 
            Y_DA_codedEnc(k1) = 0;
        end
   end
   % Cas NDA 
   Y_I_NDA = real(Y_NDA);
   for i = 1 : length(Y_I_NDA)
       if Y_I_NDA(i) > Seuil
            Y_NDA_coded(i) = 1;
        else 
            Y_NDA_coded(i) = 0;
        end
   end
   Y_I_NDAenc = real(Y_NDAenc);
   for i = 1 : length(Y_I_NDAenc)
       if Y_I_NDAenc(i) > Seuil
            Y_NDA_codedEnc(i) = 1;
        else 
            Y_NDA_codedEnc(i) = 0;
        end
   end
   % Calcul de Y_coded dans le cas où la synchronisation est absente
   r_I = real(Y_complexe);
   for k1 = 1 : length(r_I)
       if r_I(k1) > Seuil
            Y_coded(k1) = 1;
        else 
            Y_coded(k1) = 0;
        end
   end
   r_Ienc = real(Y_complexeEnc);
   for k1 = 1 : length(r_Ienc)
       if r_Ienc(k1) > Seuil
            Y_codedEnc(k1) = 1;
        else 
            Y_codedEnc(k1) = 0;
        end
   end
   % calcul de Y_coded dans le cas idéal : pas de dépahasage
   Yr_complexe = S + n;
   Y_I = real(Yr_complexe);
   for k1 = 1 : length(Y_I)
       if Y_I(k1) > Seuil
            Yr_coded(k1) = 1;
        else 
            Yr_coded(k1) = 0;
        end
   end
   Yr_complexeEnc = Senc + nEnc;
   Y_Ienc = real(Yr_complexeEnc);
   for k1 = 1 : length(Y_Ienc)
       if Y_Ienc(k1) > Seuil
            Yr_codedEnc(k1) = 1;
        else 
            Yr_codedEnc(k1) = 0;
        end
   end
   % Résultats du TEB 
   Ne_Sans_Deph = sum(abs(Yr_coded-b));
   Ne_Sans_Synch=sum(abs(Y_coded-b));
   Ne_DA=sum(abs(Y_DA_coded-b));
   Ne_NDA=sum(abs(Y_NDA_coded-b));
   
   tbdepth = 34;
    
   Yr_codedEnc = vitdec(Yr_codedEnc,trellis,tbdepth,'trunc','hard');
   Y_codedEnc = vitdec(Y_codedEnc,trellis,tbdepth,'trunc','hard');
   Y_DA_codedEnc = vitdec(Y_DA_codedEnc,trellis,tbdepth,'trunc','hard');
   Y_NDA_codedEnc = vitdec(Y_NDA_codedEnc,trellis,tbdepth,'trunc','hard');
   Ne_Sans_DephEnc = sum(abs(Yr_codedEnc-b));
   Ne_Sans_SynchEnc=sum(abs(Y_codedEnc-b));
   Ne_DAenc=sum(abs(Y_DA_codedEnc-b));
   Ne_NDAenc=sum(abs(Y_NDA_codedEnc-b));
   
   TEB_Sans_Deph(RSB+1) = Ne_Sans_Deph/Nb;
   TEB_Sans_Synch(RSB+1)=Ne_Sans_Synch/Nb;
   TEB_DA(RSB+1)=Ne_DA/Nb;
   TEB_NDA(RSB+1)=Ne_NDA/Nb;
   
   TEB_Sans_DephEnc(RSB+1) = Ne_Sans_DephEnc/Nb;
   TEB_Sans_SynchEnc(RSB+1)=Ne_Sans_SynchEnc/Nb;
   TEB_DAenc(RSB+1)=Ne_DAenc/Nb;
   TEB_NDAenc(RSB+1)=Ne_NDAenc/Nb;
end
% Représentation graphique des résultats, dans une échelle logarithmique
% TEB
figure 
semilogy(0:40,TEB_Sans_Synch,'- b')
hold on
semilogy(0:40,TEB_NDA,'- r')
hold on
semilogy(0:40,TEB_DA,'- g')
hold on
semilogy(0:40,TEB_Sans_Deph,'- k')
hold on 
semilogy(0:40,TEB_Sans_SynchEnc,'-o b')
hold on
semilogy(0:40,TEB_NDAenc,'-o r')
hold on
semilogy(0:40,TEB_DAenc,'-o g')
hold on
semilogy(0:40,TEB_Sans_DephEnc,'-o k')
xlabel('RSB (dB)')
ylabel('TEB')
grid on
legend('SANS SYNCH','SYNCH mode NDA','SYNCH mode DA','SANS DEPH','SANS SYNCH avec CC','SYNCH mode NDA avec CC','SYNCH mode DA avec CC','SANS DEPH avec CC')



