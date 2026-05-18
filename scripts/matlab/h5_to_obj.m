function h5_to_obj(surf_struct, outputDir, flip_normals)
% EXPORT_SURF_TO_OBJ Exports 3D surface structures to .obj files.
%
%   h5_to_obj(surf_struct) exports each surface in the structure
%   to a folder named 'exported_obj' in the current working directory.
%
%   h5_to_obj(surf_struct, outputDir) exports the surfaces to the 
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
    outputDir = fullfile(pwd, 'exported_obj');
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
    file_path = fullfile(outputDir, ['Surface_', entity_name, '.obj']); 
        
    try
        fid = fopen(file_path, 'w');
        if fid == -1, error('Could not open file.'); end
            
        if isfield(temp_surf, 'vertices')
            fprintf(fid, 'v %f %f %f\n', temp_surf.vertices);
        else
            error('Vertices missing.');
        end

       tr_check = double(temp_surf.triangles) + 1;
        v1 = temp_surf.vertices(:, tr_check(1,:));
        v2 = temp_surf.vertices(:, tr_check(2,:));
        v3 = temp_surf.vertices(:, tr_check(3,:));
        face_normals_calc = cross(v2 - v1, v3 - v1, 1);

        v1_normals = temp_surf.vertex_normals(:, tr_check(1,:));

        alignment = dot(face_normals_calc, v1_normals, 1);
        
        original_is_inward = mean(alignment) < 0;

        if isfield(temp_surf, 'vertex_normals')
            normals_to_write = temp_surf.vertex_normals;
            if flip_normals
                normals_to_write = -normals_to_write; 
            end
            fprintf(fid, 'vn %f %f %f\n', normals_to_write);
            has_normals = true;
        else
            has_normals = false;
        end
            
        if isfield(temp_surf, 'triangles')
            F = double(temp_surf.triangles + 1);

            if original_is_inward
                F = F([1, 3, 2], :);
            end
            
            if flip_normals
                F = F([1, 3, 2], :);
            end
                
            if has_normals
                face_data = [F(1,:); F(1,:); F(2,:); F(2,:); F(3,:); F(3,:)];
                fprintf(fid, 'f %d//%d %d//%d %d//%d\n', face_data);
            else
                fprintf(fid, 'f %d %d %d\n', F);
            end
        else
            error('Triangles missing.');
        end
        fclose(fid);
        
        status_msg = ' (Outward)';
        if flip_normals
            status_msg = ' (Inward)';
        end
        fprintf('Entity exported correctly: %s.obj%s\n', entity_name, status_msg);
            
    catch ME
        if exist('fid', 'var') && fid ~= -1, fclose(fid); end
        warning('Failed to export entity: %s\n%s', entity_name, ME.message);
    end
end
end