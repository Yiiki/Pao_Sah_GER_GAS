% main_ode_script_matR.m [required]

delay=0.5;

figure(1)
filename = './images/iv-evo.gif';

for zz=1:mag_Rlen

idsmat_edg=idsmat;

for z=1:cub_rows
    [zi,zj]=vec2squ(vdslen,vgslen,z);
    idsmat_edg(zi,zj)=idscub(z,zz);
end

plot(vdslis,idsmat_edg)
xlim([0,5])
ylim([0,2e-5])
stag=sprintf('%dth-R=%.3e',zz,mag_Rlis(zz));
title(stag)

pause(0.5)


      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if zz == 1
          imwrite(imind,cm,filename,'gif', 'DelayTime', delay, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', delay);
      end
end