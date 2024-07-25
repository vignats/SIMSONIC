function S = semblance(sig)

S = sum(sum(sig,2).^2)/sum(sum(sig.^2,2));
