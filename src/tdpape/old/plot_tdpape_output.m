function plot_tdpape_output



% cd /home/doru/infra/Cpp_code/DoruV/cpp_tests/FourierSynth_PE1


fn = 'mywavef.dat';

a = textread(fn);

%%
figure
plot(a(:,1), a(:,2))
hold on; zoom on;
xlabel('Time [sec]')


%% stack waveforms at multiple receivers
fn = 'wfgrid3.dat';
a = textread(fn);

%%
[b, I, J] = unique(a(:,1));

%%
figure
hold on; box on; zoom on;
offset = 0;
for j=1:length(b)
    K = find(J==j);
    y = a(K,3)/max(abs(a(K,3)));
    plot(a(K,2), y-offset, 'r' )
    offset = j;
end

xlabel('Time [sec]')
title(sprintf(fn))











return