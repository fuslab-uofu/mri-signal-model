function plot_sequence(t, B1, grad, sampleComb)
%% plot_sequence(t, B1, grad, sampleComb)
% Plots the given sequence waveforms in the current figure
%
% ~ Inputs ~
% * t (1,T): time base
% * B1 (1,T): RF waveform
% * grad (3,T): gradient waveforms for x y and z
% * sampleComb (1,T) (optional): time points where samples are saved
%
% 2023-05-22 Samuel Adams-Tew

lims = [min(t), max(t)];

if nargin == 3
    nplot = 4;
elseif nargin == 4 && ~any(sampleComb)
    nplot = 4;
elseif nargin == 4 && any(sampleComb)
    nplot = 5;
else
    error('Incorrect number of input arguments to plot_sequence')
end
nplot = nplot - sum(~any(grad, 2));

plotnum = 1;

subplot(nplot, 1, 1)
plot(t, abs(B1)); grid on; xlim(lims)
ylabel('RF')
plotnum = plotnum + 1;

if any(grad(1, :))
    subplot(nplot, 1, plotnum)
    plot(t, grad(1, :)); grid on; xlim(lims)
    ylabel('G_x')
    plotnum = plotnum + 1;
end
if any(grad(2, :))
    subplot(nplot, 1, plotnum)
    plot(t, grad(2, :)); grid on; xlim(lims)
    ylabel('G_y')
    plotnum = plotnum + 1;
end
if any(grad(3, :))
    subplot(nplot, 1, plotnum)
    plot(t, grad(3, :)); grid on; xlim(lims)
    ylabel('G_z')
    plotnum = plotnum + 1;
end

if nargin == 4 && any(sampleComb)
    subplot(nplot, 1, plotnum)
    stem(t(sampleComb == 1), sampleComb(sampleComb == 1), 'MarkerSize', 1); grid on; xlim(lims)
    ylabel('ADC')
end

end