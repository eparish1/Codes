function [uhat] = dealias(uhat,nx)
  kcut = ceil(nx/3);
  uhat(1:kcut,1) = 0; %set first N/3 modes = 0
  uhat(end-kcut+2:end) = 0; % set corresponding highest N/3 = 0. The +2 is since we have N/2-1 modes
end


