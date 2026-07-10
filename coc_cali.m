% Central Omnidirectional Camera Calibration Using World Point, Line and Mirror Contour
% Inputs:
%    world_point_image       - 3-by-(row * col)-by-Image_num      - World point image
%    world_point             - 4-by-(row * col)                   - World or coplane point
%    chessboard_line_image   - 3-by-3-by-(row + col)-by-Image_num - Chessboard line image
%    Image_num               - 1-by-1                             - Image number
%    mirror_contour          - 3-by-3                             - Mirror contour
%    radial_distortion_order - 1-by-1                             - Radial distortion coefficient, 0 | 1 | 2 | 3
%    l                       - 1-by-1                             - Spherical projection parameter, default is null
%
% Outputs:
%    est_K                    - 3-by-3      - Estimate intrinsic parameters
%    est_K_A                  - 3-by-3      - Estimate intrinsic parameters of the affine projection
%    est_K_M                  - 3-by-3      - Estimate intrinsic parameters of the metric projection
%    est_l                    - 3-by-3      - Estimate spherical projection parameter
%    est_k                    - order-by-1  - Estimate radial distortion coefficient
%    est_Rotation             - 3-by-3-by   - Estimate rotation
%    est_Translation          - 3-by-1      - Estimate translation
%    est_d_0                  - 1-by-1      - Estimate distance between mirror contour and sphere center
%    mean_reprojection_error  - 1-by-1      - Mean reprojection error
%    key_line_index           - 1-by-2      - The index of valid line image

function [est_K, est_K_A, est_K_M, est_l, est_k, est_Rotation, est_Translation, est_d_0, mean_reprojection_error, key_line_index] = coc_cali(world_point_image, world_point, line_image, Image_num, mirror_contour, radial_distortion_order, l) % use world point images to calibration

    arguments
       world_point_image        (:, :, :)    double 
       world_point              (:, :)       double
       line_image                            double
       Image_num                (1, 1)       double   {mustBePositive, mustBeInteger}
       mirror_contour           (3, 3)       double    
       radial_distortion_order  (1, 1)                {mustBeMember(radial_distortion_order, [0, 1, 2, 3])}    % 径向畸变系数只能为0 2 3，默认为0
       l                                              {mustBeScalarOrEmpty, mustBeGreaterThanOrEqual(l, 0)} % 镜面参数 >= 0
    end

    validate_line_images(line_image, Image_num);
    
    %% Recovery of intrinsic parameters and spherical projection parameter
    if norm(line_image, 'fro') == 0
        error('Input must be reasonable line images. ')
    end

    est_K_A = find_K_A(mirror_contour); 

    [est_K_M, est_l, key_line_index, flag] = find_fl(line_image, est_K_A, l); 

    if flag == 0
        error('You may input wrong line image, please select the appropriate line image. ')
    end

    est_K = est_K_A * est_K_M;

    if isempty(world_point_image) && isempty(world_point)
            
            est_k = [0; 0];
            est_d_0 = find_d_0(mirror_contour, est_K, est_l, est_k);

            est_Rotation = [];
            est_Translation = [];
            mean_reprojection_error = [];

            return

    elseif isempty(world_point_image) || isempty(world_point)
            error('world_point_image and world_point must be both empty or non empty at the same time.');
    else
        validate_world_points(world_point_image, world_point, Image_num);

        %% Recovery of the extrinsic parameter
        [est_Rotation, est_Translation, flag] = Deng(world_point_image, world_point, est_K, est_l); 
        if sum(flag) == 0
            error('The image cannot be calibrated by this algorithm, or you may input wrong line image. ')
        end
    
        %% Distortion
        if isempty(l)
            est_k = [0; 0];
        else 
            est_l = l;
            world_point_inf_image = world_to_inf(world_point, est_Rotation, est_Translation, est_l);
            [est_k, est_K, est_K_M] = fit_k(world_point_inf_image, world_point_image, est_K_A, est_K_M, radial_distortion_order); % 这里est_K_M可能求不到值
        end
    
        est_d_0 = find_d_0(mirror_contour, est_K, est_l, est_k);
        [est_point_image, flag] = generate_cos_world_point_image(world_point, Image_num, est_K, est_l, est_Rotation, est_Translation, est_k, []);
        mean_reprojection_error = find_mre(world_point_image, est_point_image);
    
        if flag == 0
            error('The image cannot be calibrated by this algorithm, or you may input wrong mirror contour.')
        end

    end

    
end



%% Helper function

% Recover intrinsic parameters of the affine projection from mirror contour
function est_K_A = find_K_A(mirror_contour)
    C = mirror_contour / mirror_contour(3, 3);
    a = C(1, 1);
    b = C(2, 1);
    c = C(2, 2);
    d = C(3, 1);
    e = C(3, 2);
    f = C(3, 3);

    r = sqrt(- b^2 / a^2 + c / a);
    s_a = - b / a;
    u0 = (b * e - c * d) / (a * c - b^2);
    v0 = (b * d - a * e) / (a * c - b^2);
    
    est_K_A = [r  s_a  u0;
               0  1    v0;
               0  0    1];

end

% Recover spherical projection parameter and intrinsic parameters of the metric projection from line image
function [est_K_M, est_l, key_line_index, flag] = find_fl(line_image, K_A, l_true)
    [~, ~, Line_num, Image_num] = size(line_image);
    f_e = inf(Line_num, Image_num);
    l = inf(Line_num, Image_num);
    T_f_e = inf(Line_num, Image_num);

    flag = 1;
    threshold_fe = 10; 
    threshold_l = 0.5; 
    for n = 1 : Image_num
        for ln = 1 : Line_num
            C = line_image(:, :, ln, n);
            if norm(C) ~= 0
                C_M = K_A' * C * K_A;
                C_M = C_M / C_M(3, 3);

                % Invariant of metric projection
                a = C_M(1, 1);
                b = C_M(2, 1);
                c = C_M(2, 2);
                d = C_M(3, 1);
                e = C_M(3, 2);
                f = C_M(3, 3);
               
                f_e(ln, n) = sqrt(f * (b * f - d * e) / (d * (a * e - b * d)));
                l(ln, n) = sqrt(1 - b * f / (d * e));

                if ~isreal(f_e(ln, n))
                    f_e(ln, n) = inf;
                end
                if ~isreal(l(ln, n))
                    l(ln, n) = inf;
                end


                % Pole–polar relationship
                l_0 = C_M * [0 0 1]'; 
                l_0 = l_0 / l_0(3);

                [X, Y] = crossLandC(l_0, C_M);
                mi = [X(1) Y(1) 1]';
                
                A = [real(mi(1)^2 + mi(2)^2), 1;
                     imag(mi(1)^2 + mi(2)^2), 0];

                [~, ~, V] = svd(A); 
                v = V(:, end);
                v = v / v(end);
                T_f_e(ln, n) = sqrt(1 / v(1));

                if ~isreal(T_f_e(ln, n))
                    T_f_e(ln, n) = inf;
                end

            else 
                continue
            end
        end
    end
    
    if isempty(l_true)

        [minRowValues, rowIndices] = min(abs(T_f_e - f_e)); 
        [globalMinValue, colIndex] = min(minRowValues); 
        rowIndex = rowIndices(colIndex); 

        if globalMinValue <= threshold_fe 

            est_f_e = f_e(rowIndex, colIndex);
            est_K_M = diag([est_f_e, est_f_e, 1]);
            est_l = l(rowIndex, colIndex);

            key_line_index = [colIndex, rowIndex]; 

        else

            globalMinValue
            key_line_index = [0, 0];

            est_f_e = mean(f_e(isfinite(f_e)));

            if est_f_e == 0
                flag = 0;
                est_K_M = [];
                est_l = [];
                key_line_index = [];
                return
            end

            est_K_M = diag([est_f_e, est_f_e, 1]);

            real_l = imag(l) == 0;
            l(real_l)
            est_l = mean(l(real_l));

        end
    
    else 
        [minRowValues, rowIndices] = min(abs(l - l_true)); 
        [~, colIndex] = min(minRowValues); 
        rowIndex = rowIndices(colIndex); 
        actualMin = l(rowIndex, colIndex) - l_true;
        diffMatrix = l - l_true;
        mask = abs(diffMatrix - actualMin) <= threshold_l;                 
        selected_f_e = f_e(mask);
        finite_mask = isfinite(selected_f_e);
        filtered_f_e = selected_f_e(finite_mask);
        avg_f_e = mean(filtered_f_e);
        
        est_f_e = avg_f_e;
        est_K_M = diag([est_f_e, est_f_e, 1]);
        est_l = l_true;
        
        key_line_index = [colIndex, rowIndex]; 
        
    end

end

% Find the intersections between line L[A1;B1;D1] and conic C=[a1 b1 c1 d1 e1 1]'%
function [X,Y] = crossLandC(L,C)
    A1 = L(1);
    B1 = L(2);
    D1 = L(3);
    
    C = C/C(3,3);
    a1 = C(1,1);
    b1 = 2*C(1,2);
    c1 = C(2,2);
    d1 = 2*C(1,3);
    e1 = 2*C(2,3);
    
    X = zeros(2,1);
    Y = zeros(2,1);
    %sovle
    X(1) = -(D1 + (B1*(A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;
    
    X(2) = -(D1 - (B1*(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;
    
    Y(1) = (A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));
    
    Y(2) = -(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));

end


% Recover the extrinsic parameter from world points and their images
function [est_Rotation, est_Translation, flag] = Deng(world_point_image, world_point, K, l)
    
    [~, ~, Image_num] = size(world_point_image);
    est_Rotation = zeros(3, 3, Image_num);
    est_Translation = zeros(3, Image_num);
    flag = zeros(1, Image_num);
    
    for n = 1 : Image_num
        img_pts = world_point_image(:, :, n);

        valid_cols = ~any(isnan(img_pts), 1);

        if sum(valid_cols) < 4
            flag(n) = 0;
            est_Rotation(:, :, n) = eye(3);
            est_Translation(:, n) = zeros(3, 1);
            continue;
        end

        img_pts_valid = img_pts(:, valid_cols);
        world_pts_valid = world_point(:, valid_cols);


        [M_s, flag(n)] = back_to_sphere(img_pts_valid, K, l);

        if flag(n) == 0
            est_Rotation(:, :, n) = eye(3);
            est_Translation(:, n) = zeros(3, 1);
        end

        P = normalized_DLT(world_pts_valid, M_s); 

        [est_Rotation(:, :, n), est_Translation(:, n)] = getRT(P);
    end
    
end

% Back-projecte image point onto the unit-viewing sphere
function [M_s, flag] = back_to_sphere(world_point_image, K, l)

    flag = 1; 
    [~, r_c] = size(world_point_image); 
    M_s = ones(3, r_c); 
    
    inv_K = inv(K);
    inv_K = inv_K / inv_K(3, 3);

    image_point_inf = inv_K * world_point_image; 

    for i = 1 : r_c
        mu_i = ( l + sqrt(l^2 - (l^2 - 1) * image_point_inf(:, i)' * image_point_inf(:, i)) ) / (image_point_inf(:, i)' * image_point_inf(:, i));
        M_s(:, i) = mu_i * image_point_inf(:, i) - [0; 0; l]; 

        if ~isreal(mu_i) 
            i
            flag = 0;
            return
        end
        
        if isnan(M_s(:, i))
            world_point_image(:, i)
            image_point_inf(:, i)
            mu_i
        end

    end

end


% Compute the projection matrix from 3D/2D world points to unit sphere points using normalized DLT.
function P = normalized_DLT(world_point, M_s)
    [P_r, r_c] = size(world_point); 

    if norm(world_point(3, :)) == 0 
        world_point(3, :) = [];
        P_r = 3;
    end

    T1 = compute_normalization_transform(world_point);
    world_point_normalized = T1 * world_point;
    
    T2 = compute_normalization_transform(M_s);
    M_s_normalized = T2 * M_s;
    
    % DLT
    A = []; 
    for i = 1 : r_c
        Xi = world_point_normalized(:, i);
        A = [A;
            zeros(1, P_r)   , - M_s_normalized(3, i) * Xi', M_s_normalized(2, i) * Xi';
            M_s_normalized(3, i) * Xi' , zeros(1, P_r)    , - M_s_normalized(1, i) * Xi'];
    end
    [~, ~, V] = svd(A);
    p = V(:,end);
    p = p / p(end);

    P_normalized = zeros(3, P_r);
    for i = 1 : 3
        P_normalized(i, :) = p((P_r * (i - 1) + 1) : P_r * i);
    end

    P = T2 \ P_normalized * T1;

end

% Compute the normalization transform for DLT
function T = compute_normalization_transform(points)
    [dimension, ~] = size(points);
    if dimension == 3
        target_distance = sqrt(2);
    elseif dimension == 4
        target_distance = sqrt(3);
    else
        error('The dimension parameter must be 3 or 4.');
    end
            
    points_cartesian = points(1 : dimension - 1, :);
    
    centroid = mean(points_cartesian, 2);
    
    points_centered = points_cartesian - centroid;
    
    distances = sqrt(sum(points_centered.^2, 1));
    mean_distance = mean(distances);
    
    if mean_distance > 0
        scale = target_distance / mean_distance;
    else
        scale = 1.0;
    end
    
    T = [eye(dimension - 1)     , - scale * centroid;
         zeros(1, dimension - 1), 1                ];

end

% Recover R and t from P
function [Rotation, Translation] = getRT(P)
    [~, P_c] = size(P); 
    
    mu = 0;
    for i = 1 : P_c - 1
        mu = mu + 1 / norm(P(:, i));
    end
    mu = mu / (P_c - 1);

    R_c1 = mu * P(:, 1); 
    R_c2 = mu * P(:, 2); 

    if P_c == 3
        R_c3 = cross(R_c1, R_c2); 
    end 
    if P_c == 4
        R_c3 = mu * P(:, 3); 
    end

    R = [R_c1, R_c2, R_c3];
    
    [U, ~, V] = svd(R); 
    Rotation = U * V';

    Translation = mu * P(:, end);

end

% Get the image on the plane at infinity 
function world_point_inf_image = world_to_inf(world_point, Rotation, Translation, l)

    [~, ~, Image_num] = size(Rotation);
    [~, r_c] = size(world_point);
    world_point_inf_image = zeros(3, r_c, Image_num); %Ms 非齐次坐标

    for n = 1 : Image_num
        for i = 1 : r_c
            Xs = [Rotation(:, :, n), Translation(:, n); zeros(1,3), 1] * world_point(:, i); % 世界坐标系到单位视球坐标系
            Ms = Xs(1:3) / norm(Xs(1:3));
            world_point_inf_image(:, i, n) = [Ms(1) / (Ms(3) + l), Ms(2) / (Ms(3) + l), 1]';
        end
    end
end


% Recover the radial distortion coefficient
function [est_k, est_K, est_K_M] = fit_k(world_point_inf_image, world_point_image, K_A, K_M, order)

    if order == 0
        est_k = [0; 0];
        est_K_M = K_M;
        est_K = K_A * K_M;
        return
    end

    [~, r_c, Image_num] = size(world_point_image);
    inv_K_A = inv(K_A);
    inv_K_A = inv_K_A / inv_K_A(3, 3);

    A = [];
    b = [];

    for n = 1 :Image_num

        img_pts = world_point_image(:, :, n);

        valid_cols = ~any(isnan(img_pts), 1);

        if sum(valid_cols) < 4
            continue;
        end

        img_pts_valid = img_pts(:, valid_cols);
        world_pts_inf_valid = world_point_inf_image(:, valid_cols);

        d_world_point_metric_image = inv_K_A * img_pts_valid; 

        for i = 1 : sum(valid_cols)
            m_i_inf = world_pts_inf_valid(:, i);
            d_m_i_m = d_world_point_metric_image(:, i);
    
            r_i_m = m_i_inf(1)^2 + m_i_inf(2)^2;
            R_i_m = 1;
    
            for o = 1 : order
                R_i_m = [R_i_m, r_i_m^o];
            end
        
            A = [A;
                m_i_inf(1) * R_i_m ;
                m_i_inf(2) * R_i_m ];
        
            b = [b;
                d_m_i_m(1);
                d_m_i_m(2)];
        end
    end
        
    est_X = inv(A' * A) * A' * b;

    if sum(isnan(est_X)) == 0 && est_X(1) > 0 
        est_f_e = est_X(1);
        est_K_M = diag([est_f_e, est_f_e, 1]);
        est_K = K_A * est_K_M;
    
        est_k = est_X(2:end) / est_f_e;
    else 
        est_k = [0; 0];
        est_K_M = K_M;
        est_K = K_A * K_M;
        return
    end

end


% Recover d_0
function est_d_0 = find_d_0(mirror_contour, K, l, k)

    Q = K' * mirror_contour * K;
    Q = Q / Q(3, 3);
    e = eig(Q);
    [~, j] = max(abs(median(sign(e))-sign(e))); 
    e(j) = [];
    lambda = mean(e);
    d_r_inf = sqrt(-1 / lambda);

    if norm(k) == 0
        r_inf = d_r_inf;
    else
        order = size(k);
        r = d_r_inf;
        iteration = 10; 
        epsilon = 10^-8; 
        
        for i = 1 : iteration      
            f_i = r(end) - r(1);
            f_d_i = 1;
            for o = 1 : order
                f_i = f_i + k(o) * r(end)^(2 * o + 1);
                f_d_i = f_d_i + (2 * o + 1) * k(o) * r(end)^(2 * o);
            end
            r = [r; r(end) - f_i / f_d_i];
            if abs(r(end) - r(end - 1)) < epsilon
                break
            end
        end
        r_inf = r(end);
    end

    est_d_0 = (r_inf^2 * l - sqrt(r_inf^2 - r_inf^2 * l^2 + 1)) / (r_inf^2 + 1);

end


% % Compute mean reprojection errors from real vs estimated image points.
function mre = find_mre(real_point_image, est_point_image)
    [~, r_c, Image_num] = size(real_point_image); 
    mre = zeros(1, Image_num);
    for n = 1 : Image_num

        real_img_pts = real_point_image(:, :, n);
        est_img_pts = est_point_image(:, :, n);

        valid_cols = ~any(isnan(real_img_pts), 1);

        if sum(valid_cols) < 1
            error('没有有效点数');
        end

        img_pts_valid = real_img_pts(:, valid_cols);
        est_img_pts_valid = est_img_pts(:, valid_cols);


        RE = zeros(1, sum(valid_cols)); 
        for i = 1 : length(RE)
            RE(i) = norm(img_pts_valid(1:2, i) - est_img_pts_valid(1:2, i));
        end
        mre(n) = mean(RE);

    end
end


function validate_world_points(img_pts, world_pts, ImgNum)
    if ndims(world_pts) ~= 2
        error('world_point must be a 2D array (4xN).');
    end
    if size(world_pts, 1) ~= 4
        error('The first dimension of world_point must be 4 (homogeneous coordinates).');
    end
    N_wld = size(world_pts, 2);
    if N_wld < 1
        error('world_point requires at least one point.');
    end

    sz_img = size(img_pts);
    if ndims(img_pts) == 2
        if ImgNum ~= 1
            error('ImageNum = %d > 1, but world_point_image is a 2D array, which should be 3 × N × %d.', ImgNum, ImgNum);
        end
        if sz_img(1) ~= 3
            error('The first dimension of world_point_image must be 3 (homogeneous coordinates).');
        end
        N_img = sz_img(2);
    elseif ndims(img_pts) == 3
        if ImgNum == 1
            error('ImageNum = 1, but world_point_imageis a 3D array, which should be 3 × N.');
        end
        if sz_img(1) ~= 3
            error('The first dimension of world_point_image must be 3 (homogeneous coordinates).');
        end
        if sz_img(3) ~= ImgNum
            error('The third dimension of world_point_image is %d, but ImageNum = %d', sz_img(3), ImgNum);
        end
        N_img = sz_img(2);
    else
        error('world_point_image must be 2D or 3D array.');
    end

    if N_img ~= N_wld
        error('Number of points mismatch: There are %d world points and %d image points.', N_wld, N_img);
    end
end



function validate_line_images(line_img, ImgNum)

    sz = size(line_img);
    nd = ndims(line_img);

    if sz(1) ~= 3 || sz(2) ~= 3
        error('The first two dimensions of line_image must be 3 × 3.');
    end

    if ImgNum == 1
        if nd == 2
            return;
        elseif nd == 3
            return;
        elseif nd == 4 && sz(4) == 1
            return;
        else
            error('When ImageNum = 1, line_image must be a 3 × 3 matrix, 3 × 3 × N array, or 3 × 3 × N × 1 array.');
        end
    else
        if nd == 3
            if sz(3) ~= ImgNum
                error('When line_image is 3 × 3 × K, K must be equal to ImageNum = %d, which is actually %d.', ImgNum, sz(3));
            end
        elseif nd == 4
            if sz(4) ~= ImgNum
                error('The fourth dimension of line_image is %d, but ImageNum = %d.', sz(4), ImgNum);
            end
        else
            error('When ImageNum = %d > 1, line_image must be in 3 × 3 × K or 3 × 3 × N × K format.', ImgNum);
        end
    end
end