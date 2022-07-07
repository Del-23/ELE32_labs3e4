
probs = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
ebn0 = qfuncinv(probs).^2/2;
ps = qfunc(sqrt(2*ebn0));
loglog(ebn0, ps)
xlabel("Relação sinal-ruído")
ylabel("Probabilidade de erro de bit na saída")