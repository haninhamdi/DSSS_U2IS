% Cette fonction inverse l'étalement de spectre 
function [Data] = InverseEtalementDeSpectre(spreadData,spreadSequence,facteur,N)
Data = zeros(N,1);
j=1;
    for i = 1:N
      msg_demod(j:j+facteur-1) = spreadData(j:j+facteur-1)'.*spreadSequence(1:facteur);
      j = facteur*i+1; 
    end
    j=1;
    for i = 1:N
      sum=0;
      for k = j:j+facteur-1
        sum=sum+msg_demod(k);
      end
      
      if (sum >0)
        Data(i)=1;
      else
        Data(i)=0;
      end  
       j = facteur*i+1;    
    end
end
