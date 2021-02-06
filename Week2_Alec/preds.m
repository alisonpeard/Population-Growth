function [y,u] = preds(k,u0)
    
    u = 1+(k*log(k./(1+k-u0)));
    y = -k*log(k./(1+k-u0));

end