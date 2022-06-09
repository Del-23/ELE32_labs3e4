function x = channel(data, prob)
    data_size = size(data);
    rand_vals = rand(data_size);
    err = double(rand_vals < prob);
    x = cast(xor(data, err), 'like', data);
end