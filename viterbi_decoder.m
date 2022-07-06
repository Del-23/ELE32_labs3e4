function [X, M1, M2] = viterbi_decoder(coded, trellis)
    K = trellis.numStates;
    N = trellis.numOutputSymbols;
    T = length(coded) / log2(N);
    M1 = Inf(K, T);  % dist hamming
    M2 = zeros(K, T) - 1;  % simbolo
    oct_array = bin2oct(coded);
    from_zero = trellis.nextStates(1,:) + 1;
    for i=1:length(from_zero)
        bin_output = custom_oct2poly(trellis.outputs(1, i));
        bin_obs = custom_oct2poly(oct_array(1));
        dist = pdist([bin_obs; bin_output], 'hamming') * 3;
        M1(from_zero(i), 1) = dist;
        M2(from_zero(i), 1) = 0;
    end
    for j=2:T
        for s=1:K
            bin_obs = custom_oct2poly(oct_array(j));
            [rows_states, cols_states] = find(trellis.nextStates == s - 1);
            alpha = rows_states(1);
            beta = rows_states(2);
            bin_output_alpha = custom_oct2poly(trellis.outputs(alpha, cols_states(1)));
            bin_output_beta = custom_oct2poly(trellis.outputs(beta, cols_states(2)));
            dist_alpha = pdist([bin_obs; bin_output_alpha], 'hamming') * 3;
            dist_beta = pdist([bin_obs; bin_output_beta], 'hamming') * 3;
%             bin_output_alpha = custom_oct2poly(trellis.outputs(alpha, cols_states(1)));
%             bin_output_beta = custom_oct2poly(trellis.outputs(beta, cols_states(2)));
%             dist_alpha = pdist([bin_obs; bin_output_alpha], 'hamming') * 3;
%             dist_beta = pdist([bin_obs; bin_output_beta], 'hamming') * 3;
%             dists = [M1(s, j-1) + dist_alpha; M1(s, j-1) + dist_beta];
%             [dist, dists_idx] = min(dists);
%             M1(s, j) = dist;
%             if dists_idx == 1
%                 M2(s, j) = alpha - 1;
%             else
%                 M2(s, j) = beta - 1;
%             end
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

function oct_array = bin2oct(bit_array)
    n = length(bit_array);
    oct_array = zeros(1, n/3);
    for i=1:3:n
        b = bit_array(i:i+2);
        oct_array((i+2)/3) = b(1)*2^2 + b(2)*2^1 + b(3);
    end
end

function bin = custom_oct2poly(oct)
    dec = base2dec(num2str(oct), 8);
    bin = int2bit(dec, max(1, ceil(log2(dec+1))), false)';
    bin = [zeros(1, 3 - length(bin)), bin];
end