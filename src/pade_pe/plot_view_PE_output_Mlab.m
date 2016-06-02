% function plot_view_PE_output_Mlab

return

%%
pth1 = '/home/doru/infra/Cpp_code/fromJelle/pade-pe';
fn1 = 'starter_field.dat'

a = textread(fullfile(pth1, fn1));

figure
plot(a (:,2), a(:,1), 'r')
hold on; grid on;