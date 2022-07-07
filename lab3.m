N = 1e5;
ms = [3 4 6];
Gs = [
    13, 15, 17;
    25, 33, 37;
    117, 127, 155
];

u = randi([0 1], 1, 1e4);
probs = [0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
probs_biterr = zeros(size(ms, 2), size(probs, 2));
for i=1:length(ms)
    trellis = poly2trellis(ms(i) + 1, Gs(i, :));
    v = convenc(u, trellis);
    for j=1:length(probs)
        err = rand(size(v)) < probs(j);
        z = mod(v + err, 2);
        [X, M1, M2] = viterbi_decoder(z, trellis);
        x = zeros(size(X));
        x(1) = find(trellis.nextStates(1,:) == X(1)) - 1;
        for k=2:length(X)
            x(k) = find(trellis.nextStates(X(k-1)+1,:) == X(k)) - 1;
        end
        probs_biterr(i, j) = double(1 - sum(x == u, 'all')/numel(x));
    end
end