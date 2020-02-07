function original_obstacle_layers = GenerateOriginalObstacleLayers()
global xyt_graph_search_ obstacles_
original_obstacle_layers = cell(1, xyt_graph_search_.num_nodes_t);
for ii = 1 : xyt_graph_search_.num_nodes_t
    x_obs = [];
    y_obs = [];
    for jj = 1 : size(obstacles_,2)
        elem = obstacles_{1,jj};
        for kk = 1 : (length(elem{1,1}.x) - 1)
            x = linspace(elem{1,ii}.x(kk), elem{1,ii}.x(kk+1), 20);
            y = linspace(elem{1,ii}.y(kk), elem{1,ii}.y(kk+1), 20);
            x_obs = [x_obs, x];
            y_obs = [y_obs, y];
        end
    end
    temp.x = x_obs; temp.y = y_obs;
    original_obstacle_layers{1,ii} = temp;
end
end