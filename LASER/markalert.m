%subroutine of dLTplot_synthetic_highermodes.m
%mark the alerting parts:
rtol=0.5;
for i=iplot_meso(1:end)
        faces=[1 2 3 4];
        ymin=min(ylim);        
        ymax=max(ylim);
       % ymax=ymax-(ymax-ymin)/2;
        xmin=(dist_axis(i));
        xmax=(dist_axis(i+1));
        vertices=[xmin ymin; xmax ymin; xmax ymax; xmin ymax];
        if alphafactor_var(i)>rtol^2
        patch('Faces',faces,'Vertices',vertices,...
            'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none')
        end
end
axis tight
