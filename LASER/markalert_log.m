%Subroutine to be inserted into the plotting algorithm
%Modified from:
%Synthetic\20191110\run6\dLTplot_synthyetic_hifhermodes.m
%Works to: mark the alerting parts. This code applies to log-log plots
for i=1:(length(dist_axis)-1)
        faces=[1 2 3 4];
        ymin=min(ylim);        
        ymax=max(ylim);
        
        ymin_log=log(ymin);
        ymax_log=log(ymax);
        
        ymin=exp(ymin_log+(ymax_log-ymin_log)/2);
        
        xmin=(dist_axis(i));
        xmax=(dist_axis(i+1));
        vertices=[xmin ymin; xmax ymin; xmax ymax; xmin ymax];
        patch('Faces',faces,'Vertices',vertices,...
            'FaceColor','blue','FaceAlpha',0.7*alphafactor_ratio(i)^3,'EdgeColor','none')
end
hold on
for i=1:(length(dist_axis)-1)
        faces=[1 2 3 4];
        ymin=min(ylim);        
        ymax=max(ylim);
          
        ymin_log=log(ymin);
        ymax_log=log(ymax);
        
        ymax=exp(ymax_log-(ymax_log-ymin_log)/2);

        xmin=(dist_axis(i));
        xmax=(dist_axis(i+1));
        vertices=[xmin ymin; xmax ymin; xmax ymax; xmin ymax];
        patch('Faces',faces,'Vertices',vertices,...
            'FaceColor','red','FaceAlpha',0.5*alphafactor_var(i),'EdgeColor','none')
end
axis tight