% % read audio file example
% [signal,fs] = audioread('song_viet.wav');
% signal = signal(:,1);
% %plot(psd(spectrum.periodogram,signal,'Fs',fs,'NFFT',length(signal)));
% figure;
% spectrogram(signal,128,120,128,fs,'yaxis');

% read audio file
[signal,fs] = audioread('meadowlark.wav');
signal = signal(:,1);
figure;
[s,f,t] = spectrogram(signal,128,120,256,fs,'yaxis');
spectrogram(signal,128,120,256,fs,'yaxis');

% build supervisor
supervisor = abs(s);
supervisor = log(supervisor(end:-1:1,:));%/max(supervisor,[],'all'));

% normalize to range 0-1
supervisor = (supervisor - min(supervisor,[],'all')) ;
supervisor = supervisor/max(supervisor,[],'all');
%supervisor = supervisor.^(1.2);
save('supervisor.mat','supervisor')

figure;
imagesc(supervisor)
colorbar

