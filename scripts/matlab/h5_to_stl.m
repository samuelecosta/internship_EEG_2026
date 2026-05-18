function h5_to_stl(surf_struct, outputDir, flip_normals)
% H5_TO_STL Exports 3D surface structures to .stl files.
%
%   h5_to_stl(surf_struct) exports each surface in the structure
%   to a folder named 'exported_stl' in the current working directory.
%
%   h5_to_stl(surf_struct, outputDir) exports the surfaces to the 
%   specified directory 'outputDir'.
%
%   INPUTS:
%   surf_struct: A structure where each field represents a 3D entity and 
%                contains at least 'vertices' and 'triangles'.
%   outputDir:   (Optional) Destination folder path.
% flip_normals: (Optional) Changes direction of normals

if nargin < 1 || ~isstruct(surf_struct)
    error('Not a valid structure.');
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(pwd, 'exported_stl');
end
if nargin < 3 || isempty(flip_normals)
    flip_normals = false;
end
        
if ~exist(outputDir, 'dir')
    [status, msg] = mkdir(outputDir);
    if ~status, error('Impossible to create the folder %s: %s', outputDir, msg); end
end
        
surf_names = fieldnames(surf_struct);
if isempty(surf_names)
    warning('Structure is empty.'); return; 
end

for i = 1:length(surf_names)
    entity_name = surf_names{i};
    temp_surf = surf_struct.(entity_name);
    file_path = fullfile(outputDir, ['Surface_', entity_name, '.stl']); 
            
    try
        if ~isfield(temp_surf, 'vertices'), error('Vertices missing.'); end
        if ~isfield(temp_surf, 'triangles'), error('Triangles missing.'); end
        
        tr_check = double(temp_surf.triangles) + 1;
        v1 = temp_surf.vertices(:, tr_check(1,:));
        v2 = temp_surf.vertices(:, tr_check(2,:));
        v3 = temp_surf.vertices(:, tr_check(3,:));
        face_normals_calc = cross(v2 - v1, v3 - v1, 1);

        v1_normals = temp_surf.vertex_normals(:, tr_check(1,:));

        alignment = dot(face_normals_calc, v1_normals, 1);
        
        original_is_inward = mean(alignment) < 0; 

        verts = double(temp_surf.vertices);
        faces = double(temp_surf.triangles) + 1; 
        
        if size(verts, 2) ~= 3, verts = verts'; end
        if size(faces, 2) ~= 3, faces = faces'; end

        if original_is_inward
            faces = faces(:, [1, 3, 2]);
        end

        if flip_normals
            faces = faces(:, [1, 3, 2]);
        end
        
        TR = triangulation(faces, verts);
        stlwrite(TR, file_path);
        
        status_msg = ' (Outward)';
        if flip_normals
            status_msg = ' (Inward)';
        end
        fprintf('Entity exported correctly: Surface_%s.stl%s\n', entity_name, status_msg);
                
    catch ME
        warning('Failed to export entity: %s\n%s', entity_name, ME.message);
    end
end
end