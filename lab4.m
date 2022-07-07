N = 1e5;
ms = [3 4 6];
Gs = [
    13, 15, 17;
    25, 33, 37;
    117, 127, 155
];
P = 1;
Rb = 1;
Eb = P / Rb;

s = randsample([-1 1], 1e4, true);
u = s == 1;
probs = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
ebn0 = qfuncinv(probs).^2/2;
probs_biterr = zeros(size(ms, 2), size(probs, 2));
tic
for i=1:length(ms)
    trellis = poly2trellis(ms(i) + 1, Gs(i, :));
    v = convenc(u, trellis);
    s_t = zeros(size(v));
    for k=1:length(v)
        if v(k) == 0
            s_t(k) = v(k) - 1;
        else
            s_t(k) = v(k);
        end
    end
    for j=1:length(ebn0)
        N0 = Eb/ebn0(j);
        stderror = sqrt(N0/2);
        err = stderror * randn(size(v));
        z = s_t + err;
        [X, M1, M2] = viterbi_decoder2(z, trellis);
        x = zeros(size(X));
        x(1) = find(trellis.nextStates(1,:) == X(1)) - 1;
        for k=2:length(X)
            x(k) = find(trellis.nextStates(X(k-1)+1,:) == X(k)) - 1;
        end
        probs_biterr(i, j) = double(1 - sum(x == u, 'all')/numel(x))
    end
end
timeElapsed = toc;