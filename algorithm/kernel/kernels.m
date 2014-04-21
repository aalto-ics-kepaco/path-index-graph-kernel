function [ ] = kernels(input_file_features, output_file_kernels)

    % normalization function
    normalize = @(K) K ./ sqrt(diag(K)*diag(K)');

    disp('Reading path count matrices.');

    features = load(input_file_features);
    Phi = spconvert(features)
    
    disp('Computing kernels.');

    K = full(normalize(Phi * Phi'));
    dlmwrite(output_file_kernels, K, ' ');

end