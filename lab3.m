N = 1e5;
ms = [3 4 6];
syms x;
Gs = [
    [1 0 1 1], [1 1 0 1], [1 1 1 1];
    [1 0 1 0 1], [1 1 0 1 1], [1 1 1 1 1];
    [1 0 0 1 1 1 1], [1 0 1 0 1 1 1], [1 1 0 1 1 0 1]
    ];

u = randi([1 0], 1, 1e2);
for i=1:length(ms)
    v1 = mod(conv(Gs(i, 1), u), 2);
    v2 = mod(conv(Gs(i, 2), u), 2);
    v3 = mod(conv(Gs(i, 3), u), 2);
    % precisa achar v = v1(D^n) + v2(D^n)*D + v3(D^n)*D^2
    for j=1:length(prob)
        err = rand(size(v)) < prob(j);
        z = mod(v + err, 2);
    end
end