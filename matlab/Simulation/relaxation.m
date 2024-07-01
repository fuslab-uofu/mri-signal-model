function [Rel, Rec] = relaxation(t, T1, T2, M0)
%% [Rel, Rec] = relaxation(t, T1, T2, M0)
% Computes relaxation and recovery due to T1 and T2 given t and M0
% 
% ~ Output~
% Note that, to reduce memory requirements, outputs are formatted as column
% vectors--NOT as a matrix and a column vector. Rather than using matrix 
% multiplication, relaxed magnetization can be computed as:
%   M.*Rel + Rec
%
%% 2023-05-04 Samuel Adams-Tew

e1 = exp(-t./T1);
e2 = exp(-t./T2);

% Relaxation
Rel = [e2; e2; e1];

% Recovery
Rec = zeros(size(Rel));
Rec(3, :, :, :, :, :) = M0.*(1 - e1);

end