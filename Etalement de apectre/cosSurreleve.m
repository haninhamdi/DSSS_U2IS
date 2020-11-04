function [filtre] = cosSurreleve(alpha,t,T)
if t == 0
    filtre = (4*alpha+pi*(1-alpha))/(pi*sqrt(T));
else
    a = cos(pi*t*(1+alpha)/T);
    b = sin(pi*t*(1-alpha)/T);
    c = 4*alpha*t/T;
    d = 4*alpha/(pi*sqrt(T));
    e = 1 - (4*alpha*t/T)^2;
    filtre = (a+b/c)*d/e;
end
end
