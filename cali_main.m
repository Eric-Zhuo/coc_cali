%% Main simulate_data

clc
clear
close all

%%

% Basic parameters %
fu = 720; 
fv = 600;
r = fu / fv;
s  = 0.1;
u0 = 640;
v0 = 512;
K = [fu ,  s, u0; % Intrinsic parameters, K, 3*3, 5 degrees of freedom%
     0  , fv, v0;
     0  ,  0,  1];
    
l = 1; % Spherical projection parameter, (0, +¡Þ)%

d_0 = 0.1; % Distance between mirror contour and sphere center, [0,1)%
FoV = rad2deg(2 * acos(- d_0));

% FoV = 210;
% d_0 = - cos(deg2rad(FoV) / 2);


k1 = 0; 
k2 = 0;
k = [k1; k2]; % Radial distortion coefficient, the order is 2
if norm(k) == 0 
    radial_distortion_order = 0; 
else
    radial_distortion_order = size(k, 1);
end


%% Generate simulate point and line image

% Generate chessboard point
row = 10;
col = 7;
side = 30; % Unit:mm
Image_num = 2; % Number of image


Rotation = zeros(3, 3, Image_num); 
Translation = zeros(3, Image_num);

Rotation(:, :, 1) = [-0.285719303307216,0.229691442167011,0.930379665036204;-0.378598294904411,-0.918929585006904,0.110597237282455;0.880356638421675,-0.320640389212398,0.349516709175127];
Translation(:, 1) = [1.205137167200613e+02; 31.531943555993170; 72.442303511119550];

Rotation(:, :, 2) = [-0.613512623064142,-0.778155002749227,0.134450931708587;0.783017781820781,-0.577373370357345,0.231350696032377;-0.102398313918838,0.247214042680734,0.963535055100762];
Translation(:, 2) = [-112.012830124373; 52.599969211831; 165.549759288735];


% Generate mirror contour
[mirror_contour, mirror_contour_point] = generate_mirror_contour(d_0, K, l, k); 

% Generate simulation chessboard point and line image
[chessboard_cornerpoint, chessboard_linepoint] = generate_chessboard(row, col, side);

% chessboard_linepoint
[chessboard_cornerpoint_image, chessboard_linepoint_image, flag] = generate_cos_chessboard_image(chessboard_cornerpoint, chessboard_linepoint, Image_num, K, l, Rotation, Translation, k, d_0);
if flag == 0
    error('Not all world points are included in the mirror contour. ')
end


chessboard_line_image = fit_chessboard_line_image_fromcornerpoint(row, col, Image_num, chessboard_cornerpoint_image);


%% Calibration


[est_K, est_K_A, est_K_M, est_l, est_k, est_Rotation, est_Translation, est_d_0, mean_reprojection_error, key_line_index] = coc_cali(chessboard_cornerpoint_image, chessboard_cornerpoint, chessboard_line_image, Image_num,mirror_contour, radial_distortion_order, []) % unknown l 
%[est_K, est_K_A, est_K_M, est_l, est_k, est_Rotation, est_Translation, est_d_0, mean_reprojection_error, key_line_index] = coc_cali(chessboard_cornerpoint_image, chessboard_cornerpoint, chessboard_line_image, Image_num,mirror_contour, radial_distortion_order, l) % known l 


%% Plot image

Color = [51, 153, 0;
         51, 153, 0] / 255;

plot_est_coc_image(chessboard_cornerpoint_image, chessboard_cornerpoint, mirror_contour_point, est_K, est_l, est_Rotation, est_Translation, est_k, est_d_0, Color)

%% Helper function

% Generate Random Central Catadioptric Mirror Contour
% Inputs:
%    d0 - 1-by-1 - Distance between mirror contour and sphere center
%    K  - 3-by-3 - Intrinsic parameters
%    l  - 1-by-1 - Mirror parameter
%    k  - order-by-1 - Distortion coefficient, the order is 2 or 3 usually

%
% Outputs:
%    dis_mirror_contour       - 3-by-3 - Mirror contour
%    dis_mirror_contour_point - 3-by-n - Mirror contour points

function [dis_mirror_contour, dis_mirror_contour_point] = generate_mirror_contour(d0, K, l, k)

  order = size(k);
  
  if d0 >= 1 && d0 < 0
        error('The range of d0 should be 0 <= d0 < 1.')
  end
  
  if l <= 0
        error('The range of l should be l >= 0.')
  end
  
  % mirror_contour
  d = d0;
  Q = [(d - l)^2  0              0;
       0             (d - l)^2   0;
       0             0              d^2 - 1];
  
  r = (1 - d^2) / (l - d)^2;
  dis = 1; 
  for i = 1 : order
      dis = dis + k(i) * r^i;
  end
  D = diag([1/dis, 1/dis, 1]);
  C = inv(K') * D * Q * D *inv(K);
  C = C / C(3,3);
  dis_mirror_contour = C;
  
  % mirror_contour_point
  n = 200; 
  dis_mirror_contour_point = zeros(3, n);

  r0 = sqrt(1 - d0^2);
  theta = linspace(0, 2 * pi, n);
  x0 = r0 * cos(theta);
  y0 = r0 * sin(theta);
  z0 = - r0 * d0 * ones(1, length(theta));
  m_inf = [x0; y0; z0 + l];
  
  for i = 1 : n
        m_inf(:, i) = m_inf(:, i) / m_inf(3, i);
        
        r_m = m_inf(1, i)^2 + m_inf(2, i)^2;
        k_m = 1;
        for j = 1 : size(k) 
            k_m = k_m + k(j) * r_m^(j);
        end        
        d_m_inf = [k_m; k_m; 1] .* m_inf(:, i); 
        
        m_d = K * d_m_inf; 
        dis_mirror_contour_point(:, i) = m_d / m_d(3); 
  end

end


% Generate Chessboard Point on World
% Inputs:
%    row        - 1-by-1    - Row parameter of checkerboard
%    col        - 1-by-1    - Column parameter of checkerboard
%    side       - 1-by-1    - Length parameter of checkerboard, unit mm

% Outputs:
%    chessboard_point     - 4-by-(row * col)       - Chessboard points
%    chessboard_linepoint - 4-by-N0-by-(row + col) - Chessboard linepoints

function [chessboard_cornerpoint, chessboard_linepoint] = generate_chessboard(row, col, side)
    arguments
       row         (1, 1)     double       {mustBePositive, mustBeInteger}
       col         (1, 1)     double       {mustBePositive, mustBeInteger}
       side        (1, 1)     double       {mustBePositive}
    end

    % Method needs at least 4 world points to calibration
    if row * col < 4
        error('Method needs at least 4 world points to calibration.')
    end

    % Cornerpoint
    M = zeros(3, row * col); 
    for i = 1 : col
        for m = 1 : row
            M_col = m + (i - 1) * row; 
            M(1, M_col) = (i - 1) * side;
            M(2, M_col) = (m - 1) * side;
        end
    end
    chessboard_cornerpoint = [M; ones(1, row * col)];

    % Linepoint
    N0 = 50; 
    LP = ones(4, N0, row + col); 

    for i = 1 : col
        LP(1, :, i) = (i - 1) * side * ones(1, N0);
        LP(2, :, i) = linspace(0, (row - 1) * side, N0);
        LP(3, :, i) = zeros(1, N0);
    end

    for m = 1 : row
        LP(1, :, col + m) = linspace(0, (col - 1) * side, N0);
        LP(2, :, col + m) = (m - 1) * side * ones(1, N0);
        LP(3, :, col + m) = zeros(1, N0);
    end
    chessboard_linepoint = LP;

end

