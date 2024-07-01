function plot_magnetization(xaxis, meas, Mz)
%% plot_sequence(t, B1, grad, sampleComb)
% Plots the magnetization data in the current figure
%
% ~ Inputs ~
% * xaxis (1,X): units along horizontal axis
% * meas (1,X): transverse magnetization (complex valued)
% * Mz (1,X): longitudinal magnetization
%
% 2023-05-22 Samuel Adams-Tew

lims = [min(xaxis), max(xaxis)];

subplot(3, 1, 1)
plot(squeeze(xaxis), squeeze(abs(meas)));
grid on; xlim(lims); ylim([0, 1])
ylabel('Magnitude')
subplot(3, 1, 2)
plot(squeeze(xaxis), squeeze(angle(meas)));
ylabel('Phase')
grid on;
xlim(lims); ylim([-pi, pi]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'})
subplot(3, 1, 3)
plot(squeeze(xaxis), squeeze(Mz));
ylabel('M_z')
grid on; xlim(lims); ylim([-1, 1])

end