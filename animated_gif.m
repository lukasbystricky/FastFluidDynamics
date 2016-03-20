function animated_gif(u, v, fname, dt)

vmag = sqrt(u.^2 + v.^2);

[~,~,nsteps] = size(u);

h = figure();

for i = 1:10:nsteps
   surf(u(:,:,i));
   view(0,90)
   axis off
  
   shading interp
   axis square
   set(gcf, 'renderer', 'painters');
   
   title(sprintf('t = %6.3f', dt*(i-1)));
   drawnow
   
   frame = getframe(h);
   im = frame2im(frame);
   
   [imind,cm] = rgb2ind(im,256);
   
   if i == 1;
       
       imwrite(imind,cm,fname,'gif', 'Loopcount',inf, 'DelayTime', 0);
       
   else
       
       imwrite(imind,cm,fname,'gif','WriteMode','append', 'DelayTime',0);       
   end
end
