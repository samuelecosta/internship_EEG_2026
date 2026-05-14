function [points, cells, normals, dim, nNodes, matId] = readMshFile_Samuele(fileName)
%READMSHFILE Read a Gmsh ASCII 2.2 .msh surface mesh made of triangles
%
% Output:
%   points  : Nx3 node coordinates
%   cells   : Mx3 triangle connectivity
%   normals : Mx3 unit normals
%   dim     : mesh space dimension (deduced from coordinates, here 3)
%   nNodes  : number of nodes per element (for triangles = 3)
%   matId   : Mx1 material/physical region id (first tag if present)
%
% Compatible with files like:
%   $MeshFormat
%   2.2 0 8
%   $EndMeshFormat
%   $Nodes
%   ...
%   $EndNodes
%   $Elements
%   ...
%   $EndElements

    fileID = fopen(fileName, 'r');
    if fileID == -1
        error('Cannot open file: %s', fileName);
    end

    cleanupObj = onCleanup(@() fclose(fileID));

    points = [];
    cells = [];
    matId = [];
    normals = [];
    dim = 3;
    nNodes = 3;

    %% -------------------------
    % Find $Nodes section
    %% -------------------------
    line = '';
    while ischar(line)
        line = strtrim(fgetl(fileID));
        if strcmp(line, '$Nodes')
            break;
        end
    end

    if ~strcmp(line, '$Nodes')
        error('Section $Nodes not found.');
    end

    nPoint = str2double(strtrim(fgetl(fileID)));
    if isnan(nPoint) || nPoint <= 0
        error('Invalid number of nodes.');
    end

    % Gmsh node line: node-id x y z
    points = zeros(nPoint, 3);
    nodeIds = zeros(nPoint, 1);

    for i = 1:nPoint
        vals = sscanf(fgetl(fileID), '%f')';
        if numel(vals) < 4
            error('Invalid node line at node %d.', i);
        end
        nodeIds(i) = vals(1);
        points(i, :) = vals(2:4);
    end

    line = strtrim(fgetl(fileID));
    if ~strcmp(line, '$EndNodes')
        error('Expected $EndNodes, found: %s', line);
    end

    % If node IDs are not 1:N, create remapping
    maxNodeId = max(nodeIds);
    if isequal(nodeIds(:), (1:nPoint)')
        id2idx = [];
    else
        id2idx = zeros(maxNodeId, 1);
        id2idx(nodeIds) = 1:nPoint;
    end

    %% -------------------------
    % Find $Elements section
    %% -------------------------
    line = '';
    while ischar(line)
        line = strtrim(fgetl(fileID));
        if strcmp(line, '$Elements')
            break;
        end
    end

    if ~strcmp(line, '$Elements')
        error('Section $Elements not found.');
    end

    nElemTot = str2double(strtrim(fgetl(fileID)));
    if isnan(nElemTot) || nElemTot <= 0
        error('Invalid number of elements.');
    end

    % We only keep triangular surface elements: elm-type = 2
    tempCells = zeros(nElemTot, 3);
    tempMatId = zeros(nElemTot, 1);
    nTri = 0;

    for i = 1:nElemTot
        vals = sscanf(fgetl(fileID), '%f')';

        if numel(vals) < 4
            error('Invalid element line at element %d.', i);
        end

        elmType = vals(2);
        nTags   = vals(3);

        % Only triangles
        if elmType ~= 2
            continue;
        end

        % Format:
        % [elm-number elm-type number-of-tags tags... node1 node2 node3]
        if numel(vals) < 3 + nTags + 3
            error('Invalid triangle definition at element %d.', i);
        end

        tags = vals(4 : 3+nTags);
        conn = vals(4+nTags : 6+nTags);

        nTri = nTri + 1;
        tempCells(nTri, :) = conn(:)';

        if ~isempty(tags)
            % First tag in Gmsh is usually the physical entity
            tempMatId(nTri) = tags(1);
        else
            tempMatId(nTri) = 0;
        end
    end

    line = strtrim(fgetl(fileID));
    if ~strcmp(line, '$EndElements')
        error('Expected $EndElements, found: %s', line);
    end

    cells = tempCells(1:nTri, :);
    matId = tempMatId(1:nTri);

    % Remap node IDs if necessary
    if ~isempty(id2idx)
        cells = id2idx(cells);
    end

    %% -------------------------
    % Normals
    %% -------------------------
    v1 = points(cells(:,2), :) - points(cells(:,1), :);
    v2 = points(cells(:,3), :) - points(cells(:,1), :);
    normals = cross(v1, v2, 2);

    nrm = vecnorm(normals, 2, 2);
    zeroMask = (nrm == 0);
    nrm(zeroMask) = 1;  % avoid division by zero
    normals = normals ./ nrm;

    %% -------------------------
    % Output metadata
    %% -------------------------
    dim = size(points, 2);   % should be 3
    nNodes = size(cells, 2); % for triangles = 3
end