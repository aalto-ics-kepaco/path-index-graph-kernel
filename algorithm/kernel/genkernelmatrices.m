function [] = genkernelmatrices(mtldir, outputdir, eccodes)
%
% generates following kernel matrices out of the suffix tree sparse data:
%
% - l=50, lambda=0.90, forward dp
% -                    backward
% -                    f tensor b
% -                    max(ff,fb,bf,bb)
%  and all of these also as lin+quad (K + K.^2)
%  and all of these also as ind/basic
%
% 
% what we want is 
%  l={15,50}, lambda={0.90}, dp={bin,dp}, order={lin,quad}, type={f,fb,tensor,max}, paths={all,core}
%  results in 2*2*2*4*2=64 kernels
%
%
%
%


% normalization function
normalize = @(K) K ./ sqrt(diag(K)*diag(K)');

% test function
validkernel = @(K) max(K(:)) < 1.001 & min(eig(K)) > -0.0001;



% read ec-codes
ecs = dlmread(eccodes);
% get non-zero ecs
inds = find(ecs(:,1) > 0);
% take forward and backward reactions separately
forwards = inds(1:2:length(inds));
backwards = inds(2:2:length(inds));
% keep only those rows
%ecs = ecs(inds,:);
ecs = ecs(forwards,:);

NNN = 17430;
NN = 15566; % graphs
N = 7783;
D = 20665838; % features
minlen = 1;
maxlen = 50;


disp('reading path count matrices..');


% cell array containing different partial feature sparse matrices
CA = cell(50,3);
CAc = cell(50,3);
tic;
% read data length-wise
for i=minlen:maxlen
  temp = load(strcat(mtldir, 'result-kegg_', num2str(i), '.mtl'));
  CA{i,1} = sparse(NN, D);
  CA{i,1} = spconvert(temp)';
  CA{i,2} = CA{i,1}(forwards,:);
  CA{i,3} = CA{i,1}(backwards,:);
  
  temp = load(strcat(mtldir, 'result-kegg-core_', num2str(i), '.mtl'));
  CAc{i,1} = sparse(NN,D);
  CAc{i,1} = spconvert(temp)';
  CAc{i,2} = CAc{i,1}(forwards,:);
  CAc{i,3} = CAc{i,1}(backwards,:);
end
toc
disp('constructing kernels..');


% huge loop-de-loop constructing all kinds of kernels, 64 in total

lambda = 0.90;
for l=[15,50]
  for dp={'ind','float'}
    for paths={'allpaths','corepaths'}
     
     if strcmp(paths,'allpaths')
      CAcurr = CA;
     else
      CAcurr = CAc;
     end
     
     [dp ' ' num2str(l) ' ' paths]
     
     tic;
     
     % initialize kernel matrix with small constant terms to prevent zeros on diagonal
     Kff = 0.001 * ones(N,N);
     Kfb = 0.001 * ones(N,N);
     Kbf = 0.001 * ones(N,N);
     Kbb = 0.001 * ones(N,N);
     
     for k=1:l
      if strcmp(dp,'ind')
       Kff = Kff + (lambda^k * (CAcurr{k,2}>0)*(CAcurr{k,2}>0)');
       Kfb = Kfb + (lambda^k * (CAcurr{k,2}>0)*(CAcurr{k,3}>0)');
       Kbf = Kbf + (lambda^k * (CAcurr{k,3}>0)*(CAcurr{k,2}>0)');
       Kbb = Kbb + (lambda^k * (CAcurr{k,3}>0)*(CAcurr{k,3}>0)');
      else
       Kff = Kff + (lambda^k * CAcurr{k,2}*CAcurr{k,2}');
       Kfb = Kfb + (lambda^k * CAcurr{k,2}*CAcurr{k,3}');
       Kbf = Kbf + (lambda^k * CAcurr{k,3}*CAcurr{k,2}');
       Kbb = Kbb + (lambda^k * CAcurr{k,3}*CAcurr{k,3}');
      end
     end
     
     toc
     
     % normalization of crossdirections: 
     % Kfb = Kfb ./ sqrt(diag(Kbb)*diag(Kff)'); 
     % Kbf = Kfb ./ sqrt(diag(Kff)*diag(Kbb)');
     
%     % 'f'
%     K = normalize(Kff);
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_f_1.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
%     K = normalize(K.^2);
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_f_2.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
     
%     % 'fb'
%     K = Kfb ./ sqrt(diag(Kbb)*diag(Kff)');
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_fb_1.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
%     K = K.^2;
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_fb_2.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
     
     % 'tensor'
%     K = normalize(Kff.*Kbb + Kfb.*Kbf + Kbf.*Kfb + Kbb.*Kff);
%%     K = normalize(normalize(Kff) .* normalize(Kbb));
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_tensor_1.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
%     K = normalize(K.^2);
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_tensor_2.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
     
     % 'max'
%     Kcross = Kfb ./ sqrt(diag(Kbb)*diag(Kff)');
%     K = normalize(max(max(Kff, Kbb), Kfb));
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_max_1.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
%     K = normalize(K.^2);
%     filename = strcat(outputdir, sprintf('%s_%d_%d_%s_max_2.kernel', paths{1},l,round(lambda*100), dp{1}));
%     dlmwrite(filename, K);
     
    end
 end
end


