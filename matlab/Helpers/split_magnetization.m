function [meas, Mz] = split_magnetization(M)
meas = permute(M(1,:,:,:,:,:) + 1i*M(2,:,:,:,:,:), [2, 3, 4, 5, 6, 1]);
Mz = permute(M(3,:,:,:,:,:), [2, 3, 4, 5, 6, 1]);
end