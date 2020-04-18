

clear
load movie_data.mat

N = length(EPSILON)
j = 1
for j=1:N
    clf
    plot(DT,squeeze(abs(EIGS(j,:,:))));
    hold on
    plot(DT,ones(size(DT)),'--')
    yy = 0:.01:2;
    % plot(dt_thresh(j)*ones(size(yy)),yy,'--')
    % TITLE = strcat('all eigenvalues, epsilon = ', num2str(EPSILON(j)), ', dt_{thresh} = ',num2str(dt_thresh(j)));
    TITLE = strcat('all eigenvalues, epsilon = ', num2str(EPSILON(j)));
    title(TITLE)
    axis([0,.04,0,1.2])
    xlabel('dt')
    figure(1)
    pause(1)
end
