% Central Omnidirectional System Calibration Using Chessboard Point, Line and Mirror Contour
% Inputs:
%    row                        - 1-by-1                              - Row parameter of checkerboard
%    col                        - 1-by-1                              - Column parameter of checkerboard
%    Image_num                  - 1-by-1                              - Image number
%    chessboard_cornerpoint_image - 3-by-(row + col)-by-Image_num - Chessboard point image
%
% Outputs:
%    chessboard_cornerpoint_image  - 3-by-3-by-(row + col)-by-Image_num - Chessboard line image


function chessboard_line_image = fit_chessboard_line_image_fromcornerpoint(row, col, Image_num, chessboard_cornerpoint_image)

    arguments
       row                          (1, 1)        double       {mustBePositive, mustBeInteger}
       col                          (1, 1)        double       {mustBePositive, mustBeInteger}
       Image_num                    (1, 1)        double       {mustBePositive, mustBeInteger}
       chessboard_cornerpoint_image (3, :, :)  double          {mustBeChessboardCornerpointsImage(chessboard_cornerpoint_image, row, col, Image_num)}
    end

    chessboard_line_image = zeros(3, 3, row + col, Image_num);
    %% Line image
    for n = 1 : Image_num    

        for i = 1:col
            col_indices = (i-1)*row + 1 : i*row;
            line_image_point = chessboard_cornerpoint_image(:, col_indices, n);
            chessboard_line_image(:, :, i, n) = fit_Ellipse(line_image_point);
        end
        
        for r = 1:row
            row_indices = r : row : row*col;
            line_image_point = chessboard_cornerpoint_image(:, row_indices, n);
            chessboard_line_image(:, :, r+col, n) = fit_Ellipse(line_image_point);
        end

    end

end

%%
% Generate Chessboard Point on World
% Inputs:
%    X   - 3-by-n    - n image points on line image
%
% Outputs:
%    C - 3-by-3  - Line image
function C = fit_Ellipse(X)
    X_bar = [mean(X(1, :)); mean(X(2, :)); 1]; 
    T = [1, 0, - X_bar(1);
         0, 1, - X_bar(2);
         0, 0, 1];
    X_T = T * X; 
    mean_distance = mean(sqrt(sum(X_T(1:2, :).^2, 2)));
    if mean_distance > 0
        scale = sqrt(2) / mean_distance;
    else
        scale = 1.0;
    end
    S = diag([scale, scale, 1]);
    X_A = S * X_T;

    A = S * T;

    x = X_A(1, :)';
    y = X_A(2, :)';
    
    % x,y are lists of coordinates
    % Build design matrix
    D = [x.*x x.*y y.*y x y ones(size(x))];
    % Build scatter matrix
    S = D' * D;
    % Build 6x6 constraintmatrix
    C(6,6) = 0; C(1,3) = 2; C(2,2) = -1; C(3,1) = 2;
    % Solve eigensystem 
    [gevec, geval]= eig(S, C); 
    % Find the positive eigenvalue 
    [~, PosC] = find(geval>0 & ~isinf(geval));
    % Extract eigenvector corresponding to positive eigenvalue
    c = gevec(:, PosC);
    if isempty(c) 
        C = zeros(3, 3);
        return
    end
    c = c / c(end);
    C_A = [c(1)     c(2)/2  c(4)/2;
         c(2)/2   c(3)    c(5)/2;
         c(4)/2   c(5)/2  c(6)  ];
    
    C = A' * C_A * A;
    C = C / C(3, 3);
end


function mustBeChessboardCornerpointsImage(chessboard_cornerpoint_image, row, col, Image_num)
    expected_cols = row * col;
    actual_size = size(chessboard_cornerpoint_image);
    
    if actual_size(1) ~= 3 || ...
       (numel(actual_size) >= 2 && actual_size(2) ~= expected_cols) || ...
       (numel(actual_size) >= 3 && actual_size(3) ~= Image_num)
        
        eid = 'mustBeChessboardPoints:invalidSize';
        msg = sprintf('Chessboard points must be 3×%d×%d (row * col = %d+%d, Image_num = %d). Actual size: %s', ...
                      expected_cols, Image_num, row, col, Image_num, mat2str(actual_size));
        throwAsCaller(MException(eid, msg));
    end
    
    if numel(actual_size) < 2
        eid = 'mustBeChessboardPoints:insufficientDims';
        msg = sprintf('Expected at least 2 dimensions, but got %d', numel(actual_size));
        throwAsCaller(MException(eid, msg));
    end
end