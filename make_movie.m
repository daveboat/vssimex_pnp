clear
close all;

load movie_data.mat

vidObj = VideoWriter('stability_roots.avi');
vidObj.FrameRate = 5;
open(vidObj);

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
    TITLE = strcat('all eigenvalues, epsilon =   ', num2str(EPSILON(j),'%4.4f'));
    title(TITLE)
    axis([0,.04,0,1.2])
    xlabel('dt')
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    figure(1)
    % pause(1)
end

close(vidObj);








% % clear; 
% close all;
% 
% vidObj = VideoWriter('peaks.avi');
% open(vidObj);
% 
% % Creates a 2D Mesh to plot surface
% x=linspace(0,1,100);
% [X,Y] = meshgrid(x,x);
% 
% N=100; % Number of frames
% for i = 1:N
%     % Example of plot 
%     Z = sin(2*pi*(X-i/N)).*sin(2*pi*(Y-i/N));
%     surf(X,Y,Z)
%     
%            % Write each frame to the file.
%        currFrame = getframe(gcf);
%        writeVideo(vidObj,currFrame);
%     
%     % Store the frame   
%    %  M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie. 
% end 
% 
% % Output the movie as an avi file
% close(vidObj);