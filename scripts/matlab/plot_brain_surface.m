function plot_brain_surface(c_v, c_f, normals, vertex_val, title_text, view_angles, max_value, true_pos)
    % Plots a 3D brain surface with mapped values
    
    % --- correct format
    if size(c_v,2) ~= 3, c_v = c_v.'; end
    if size(c_f,2) ~= 3, c_f = c_f.'; end
    
    c_v = double(c_v);
    c_f = double(c_f);
    vertex_val = vertex_val(:); % Ensure column vector
    
    p1 = patch('Vertices', c_v, ...
          'Faces', c_f, ...
          'FaceVertexCData', vertex_val, ... 
          'FaceColor', 'interp', ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    
    % Assign pre-calculated normals
    p1.VertexNormals = normals;
    
    hold on;
    axis equal;
    axis vis3d off;
    view(view_angles(1), view_angles(2));
    set(gca, 'Projection', 'perspective');
    
    % --- colors
    colormap(gca, turbo);
    clim([0, max_value]);
    
    title(title_text, 'FontSize', 12);
    
    if nargin >= 8 && ~isempty(true_pos)
        scatter3(true_pos(1), true_pos(2), true_pos(3), 200, 'k', 'filled', 'p', 'MarkerEdgeColor', 'w');
    end
    
    % --- lighting
    camlight headlight;
    lighting gouraud;
end