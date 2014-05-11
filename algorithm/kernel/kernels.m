function [K ] = kernels(input_file_features, output_file_kernels)

    % normalization function
    normalize = @(K) K ./ sqrt(diag(K)*diag(K)');

    disp('Reading path count matrix.');

    features = load(input_file_features);
    Phi = spconvert(features);
    
    disp('Computing kernels using cosine distance.');

    K = full(1 - squareform(pdist(Phi,'cosine')));

    dlmwrite(output_file_kernels, K, ' ');

end