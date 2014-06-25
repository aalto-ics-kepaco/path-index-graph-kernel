function [ K ] = kernels(input_file_features, output_file_kernels)

    % normalization function
    normalize = @(K) K ./ sqrt(diag(K)*diag(K)');

    disp('Reading path count matrix.');

    features = load(input_file_features);
    Phi = spconvert(features);
    
    tic
    
    % euclidean distance
    disp('Computing kernels using euclidean distance.');
    Q=repmat(dot(Phi,Phi,2),1,size(Phi,1));
    D=sqrt(Q+Q'-2*Phi*Phi');
    
    scaledD = (D-min(D(:))) ./ (max(D(:)-min(D(:))));
    
    K = full(1- scaledD);
    
    toc

    dlmwrite(output_file_kernels, K, ' ');

end