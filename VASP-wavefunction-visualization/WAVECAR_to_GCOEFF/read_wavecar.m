clc;clear all;

fid = fopen('C:\Users\almul\Desktop\ORNL\nanocluster\108SiP\WAVECAR'); %Give path to WAVECAR to be read.
a = 394492; %Maximum sized plane-wave basis. To find correct value, use WaveTrans90.exe in directory with WAVECAR, and just see the max plane-wave basis it prints while running.
wavecar_data = fread(fid,[2*a Inf],'real*4')';
fclose('all')
% I know this looks dumb, but fopen() has to be used twice in a row like this.
fid = fopen('C:\Users\almul\Desktop\ORNL\nanocluster\108SiP\WAVECAR'); %Give path to WAVECAR to be read
wavecar_data_indices = fread(fid,[a Inf],'real*8')';
fclose('all')

%%% COEFFICIENTS STORED AS 2-TUPLES
for i=1:2
    wavecar_data(i,1:a) = wavecar_data_indices(i,:);
end
for i=0
    wavecar_data(i*9+3,1:1036) = wavecar_data_indices(i*9+3,1:1036) ;
end

% If you try setting manual row size, matlab picks
% the minimum number of columns required to
% hold all data.
%a=size(wavecar_data)
%a(1)*a(2) 



