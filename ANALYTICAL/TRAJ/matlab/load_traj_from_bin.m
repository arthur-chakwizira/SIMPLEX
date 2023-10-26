function [Npart, T, Nt, X,Y,Z] =  load_traj_from_bin(r_fn)
if isfile(r_fn)
    fileID = fopen(r_fn);
    Npart = fread(fileID, [1,1], 'int32');
    T = fread(fileID, [1,1], 'single');
    Nt = fread(fileID, [1,1], 'int32');
    tmp_x = fread(fileID, [Npart*Nt,1], 'single');
    tmp_y = fread(fileID, [Npart*Nt,1], 'single');
    tmp_z = fread(fileID, [Npart*Nt,1], 'single');
    fclose(fileID);
else
    error("File " + r_fn + " not found.")
end

X = zeros(Npart, Nt);
Y = zeros(size(X));
Z = zeros(size(X));

for c = 1:Npart
    idx = ((c-1)*Nt+1):c*Nt;
    X(c,:) = tmp_x(idx);
    Y(c,:) = tmp_y(idx);
    Z(c,:) = tmp_z(idx);
end


% %crop?
% Npart = 100000;
% idx = 1:Npart;
% X = X(idx, :);
% Y = Y(idx, :);
% Z = Z(idx, :);

end