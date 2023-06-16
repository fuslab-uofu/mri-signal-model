function Rz = zrot(theta)
% Right handed rotation around z axis
O = zeros(size(theta));
Rz = [cos(theta), -sin(theta),   O
      sin(theta),  cos(theta),   O
               O,           O, 1+O];
end