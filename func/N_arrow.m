function handles = N_arrow(x1,x2,y1,y2,linewidth,color,arrowpos,ind,direction)

% This function version is by by Michael Lindner, Raul Vicente, Michael Wibral
% Frankfurt 2009 Version 1.0 and it can be found in TRENTOOL TEarrow.m
%I modified the function for the need of DDIFTOOL toolbox

% =============================================
% calculate the arrow head coordinates
% =============================================
den         = x2 - x1 + eps;
teta        = atan((y2-y1)/den)+pi*(x2<x1)-pi/2;
cs          = cos(teta);
ss          = sin(teta);
R           = [cs -ss;ss cs];
%linelength  = sqrt((y2-y1)^2+(x2-x1)^2);
headlength  = .05;%min(linelength*alpha,maxlength);
headwidth   = .025;%min(linelength*beta,maxlength);
if arrowpos == 1
    x0          = x2*cs + y2*ss;
    y0          = -x2*ss + y2*cs;
elseif arrowpos == 2
    x0          = 0.5*(x1+x2)*cs + 0.5*(y1+y2)*ss;
    y0          = -0.5*(x1+x2)*ss + 0.5*(y1+y2)*cs;
end
coords      = R*[x0 x0+headwidth/2 x0-headwidth/2; y0 y0-headlength y0-headlength];

% =============================================
% plot arrow  (= line + patch of a triangle)
% =============================================
%xf,yf,0.5*(xt+xf),0.5*(yt+yf)
%putative sync is not use now (unreliable estimate)
if 1 && strcmp(direction,'putative sync')
    p=brewermap(100,'YlOrRd');


    h1          = plot([x1,x2],[y1,y2],'k','Linewidth',linewidth,'Color',p(ind,:));
    cb = colorbar; 
    cb.Label.String = 'percentage of reconstructability';
    caxis([0 100]) % sets colorbar limits 
    %set(cb,'XTicks',[0,10])
    handles = h1;
    handles =h1;
elseif strcmp(direction,'bidirectional')
    p=brewermap(100,'YlOrRd');
    h1          = plot([x1,x2],[y1,y2],'k','Linewidth',linewidth,'Color','k');

    arrow([x1,y1],[x2,y2],'EdgeColor',p(ind(1),:),'FaceColor',p(ind(1),:));
    arrow([x2,y2],[x1,y1],'EdgeColor',p(ind(1),:),'FaceColor',p(ind(2),:));
    cb = colorbar; 
    cb.Label.String = 'reconstructability in % (inverse NRMSE)';
    caxis([0 100]) % sets colorbar limits 
    %set(cb,'XTicks',[0,10])
    handles = h1;
    
    handles = h1;
%if else direct information, so arrow in one direction    
else
    p=brewermap(100,'YlOrRd');
%     p=colormap('parula');
    h1          = plot([x1,x2],[y1,y2],'k','Linewidth',linewidth,'Color',p(ind,:));

    arrow([x1,y1],[x2,y2],'EdgeColor',p(ind,:),'FaceColor',p(ind,:));
    colormap(p)
    cb = colorbar; 
    cb.Label.String = 'reconstructability in % (inverse NRMSE)';
    caxis([0 100]) % sets colorbar limits 
    
    handles = h1;
end

end
