function [] = randomsample(input_targets, input_kernels, output_targets, output_kernels, size)

    shuffle = linspace(1, size, size);
    shuffle = shuffle(randperm(size));
    inputtargets = load(input_targets);
    inputkernels = load(input_kernels);
    targets_sample = inputtargets(shuffle, :);
    kernels_sample = inputkernels(shuffle, shuffle);

    dlmwrite(output_targets, targets_sample, ' ');
    dlmwrite(output_kernels, kernels_sample, ' ');
    
end