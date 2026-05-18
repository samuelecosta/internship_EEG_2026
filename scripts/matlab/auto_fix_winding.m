function [faces_out, is_flipped] = auto_fix_winding(verts, faces)
    % AUTO_FIX_WINDING checks normal direction and fixes
    
    center_of_mass = mean(verts, 1);

    v1 = verts(faces(:, 1), :);
    v2 = verts(faces(:, 2), :);
    v3 = verts(faces(:, 3), :);
    
    face_centers = (v1 + v2 + v3) / 3;
    outward_vectors = face_centers - center_of_mass;
    
    face_normals = cross(v2 - v1, v3 - v1, 2);
    
    dot_products = sum(face_normals .* outward_vectors, 2);
    
    if sum(dot_products < 0) > sum(dot_products > 0)
        faces_out = faces(:, [1, 3, 2]);
        is_flipped = true;
    else
        faces_out = faces;
        is_flipped = false;
    end
end