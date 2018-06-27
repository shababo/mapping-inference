function [fig_handle]=visualize_3D_GP(fig_handle,this_sample,locations)

threshold = 0.05;
alpha_scale = ( (this_sample-min(this_sample)) ./range(this_sample)) ;

alpha_scale(alpha_scale<threshold)=threshold;

for i=1:length(this_sample)
    if alpha_scale(i)>threshold
        scatter3(locations(i,1),locations(i,2),locations(i,3),...
            'MarkerFaceAlpha',alpha_scale(i),...
            'MarkerFaceColor','k','MarkerEdgeColor','k')
        hold on;
    end
end
xlabel('x');ylabel('y');zlabel('z');






