lgd = cell(size(probs_biterr, 1)+1, 1);
lgd{1} = "Sem codificação";
loglog(probs, probs)
hold on
for i=1:size(probs_biterr, 1)
    loglog(probs, probs_biterr(i, :))
    lgd{i+1} = strcat("m=", num2str(ms(i)));
end
legend(lgd)
xlabel("Probabilidade de erro de bit no canal")
ylabel("Probabilidade de erro de bit na saída")
set(gca, 'xdir', 'reverse')