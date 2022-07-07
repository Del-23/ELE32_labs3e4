function [X, M1, M2] = viterbi_decoder2(from_channel, trellis)
    K = trellis.numStates;
    N = trellis.numOutputSymbols;
    T = length(from_channel) / log2(N);
    M1 = Inf(K, T);  % dist hamming
    M2 = zeros(K, T) - 1;  % simbolo
    from_zero = trellis.nextStates(1,:) + 1;
    for i=1:length(from_zero)
        bin_output = custom_oct2poly(trellis.outputs(1, i));
        actual_value = zeros(size(bin_output)) - 1;
        actual_value(bin_output == 1) = 1; 
        obs = from_channel(1:3);
        dist = pdist([obs; actual_value], 'squaredeuclidean');
        M1(from_zero(i), 1) = dist;
        M2(from_zero(i), 1) = 0;
    end
    for j=2:T
        for s=1:K
            obs = from_channel(((j-1)*3 + 1):(3*j));
            [rows_states, cols_states] = find(trellis.nextStates == s - 1);
            alpha = rows_states(1);
            beta = rows_states(2);
            bin_output_alpha = custom_oct2poly(trellis.outputs(alpha, cols_states(1)));
            bin_output_beta = custom_oct2poly(trellis.outputs(beta, cols_states(2)));

            actual_value_alpha = zeros(size(bin_output_alpha)) - 1;
            actual_value_alpha(bin_output_alpha == 1) = 1;
            actual_value_beta = zeros(size(bin_output_beta)) - 1;
            actual_value_beta(bin_output_beta == 1) = 1;

            dist_alpha = pdist([obs; actual_value_alpha], 'squaredeuclidean');
            dist_beta = pdist([obs; actual_value_beta], 'squaredeuclidean');

            if M2(alpha, j-1) ~= -1 && M2(beta, j-1) ~= -1
                dists = [M1(alpha, j-1) + dist_alpha; M1(beta, j-1) + dist_beta];
                [dist, dists_idx] = min(dists);
                M1(s, j) = dist;
                if dists_idx == 1
                    M2(s, j) = alpha - 1;
                else
                    M2(s, j) = beta - 1;
                end
            elseif M2(alpha, j-1) ~= -1
                M1(s, j) = M1(alpha, j-1) + dist_alpha;
                M2(s, j) = alpha - 1;
            elseif M2(beta, j-1) ~= -1
                M1(s, j) = M1(beta, j-1) + dist_beta;
                M2(s, j) = beta - 1;
            end
        end
    end
    Z = zeros(1, T);
    X = zeros(1, T);
    [~, Z(T)] = min(M1(:, T));
    X(T) = Z(T) - 1;
    for j=T:-1:2
        Z(j-1) = M2(Z(j), j) + 1;
        X(j-1) = Z(j-1) - 1;
    end
end

function bin = custom_oct2poly(oct)
    dec = base2dec(num2str(oct), 8);
    bin = int2bit(dec, max(1, ceil(log2(dec+1))), false)';
    bin = [zeros(1, 3 - length(bin)), bin];
end
