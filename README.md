<h1 align="center">中心全景摄像机标定 Central Omnidirectional Camera Calibration (COC)</h1>

<p align="center">
  <a href="#-chinese">🇨🇳 中文</a> •
  <a href="#-english">🇬🇧 English</a> 
</p>

---

<a id="-chinese"></a>
## 🇨🇳 中文

### 简介

本仓库提供 MATLAB 代码，用于标定中央全向相机（折反射或透射式），利用棋盘格角点、直线图像和镜面轮廓。

### 功能特点

- 估计内参矩阵 \( K \)（焦距、主点、倾斜因子）
- 估计球面投影参数 \( l \)
- 估计径向畸变系数（阶数 0–3）
- 估计多视图外参（旋转 \( R \) 和平移 \( t \)）
- 从镜面轮廓恢复镜面参数 \( d_0 \)
- 计算重投影误差以评估标定精度
- 支持已知或未知 \( l \) 的情况

### 文件结构

| 文件 | 说明 |
|------|------|
| `cali_main.m` | 主脚本，生成合成数据（棋盘格点、直线、镜面轮廓）并执行标定 |
| `coc_cali.m` | 核心标定函数，整合点、线和轮廓约束 |
| `fit_chessboard_line_image_fromcornerpoint.m` | 从检测到的角点拟合各行/列的直线图像（二次曲线） |
| `plot_est_coc_image.m` | 可视化估计结果，显示镜面轮廓、主点和重投影点 |

> **注意：** 内部还使用了其他辅助函数（如 `generate_cos_chessboard_image`、`generate_cos_world_point_image` 等），它们应位于同一目录下。

### 依赖环境

- MATLAB R2019b 或更高版本（已在 R2020a 上测试）
- 无需额外工具箱（使用内置函数如 `svd`、`eig` 等）

### 使用方法

1. 在 MATLAB 中打开项目文件夹。
2. 运行主脚本：
   ```matlab
   cali_main
   ```
   该脚本将：
   - 生成合成的棋盘格世界点和直线采样点。
   - 在中央全向相机模型下模拟图像（可配置内参、位姿和畸变）。
   - 调用标定函数 `coc_cali` 恢复所有参数。
   - 绘制估计图像并显示重投影点。

3. 若要标定自己的数据，请修改 `coc_cali` 的输入参数：
   - `world_point_image` ：观测到的图像点，3×N×M 数组（齐次坐标）
   - `world_point`        ：3D 世界点，4×N 数组（齐次坐标）
   - `line_image`         ：直线图像二次曲线，3×3×L×M 数组（每个 3×3 矩阵）
   - `mirror_contour`     ：镜面轮廓，3×3 二次曲线矩阵
   - `radial_distortion_order` ：0、1、2 或 3
   - `l`                  ：可选，若已知则传入标量，否则留空

   详细说明请参阅 `coc_cali.m` 中的函数注释。

### 引用

若您在研究中使用本代码，请引用以下工作（请根据实际参考文献补充）。

### 许可证

本项目采用 MIT 许可证，详情请见 [LICENSE](LICENSE) 文件。

---

<a id="-english"></a>
## 🇬🇧 English

### Introduction

This repository provides MATLAB code for calibrating central omnidirectional cameras (catadioptric or dioptric) using world points, line images, and the mirror contour.

### Features

- Estimate intrinsic parameters \( K \) (focal lengths, principal point, skew)
- Estimate spherical projection parameter \( l \)
- Estimate radial distortion coefficients (order 0–3)
- Estimate extrinsic parameters (rotation \( R \) and translation \( t \)) for multiple views
- Recover mirror parameter \( d_0 \) from the mirror contour
- Compute reprojection errors for evaluation
- Support both known and unknown \( l \) cases

### File Structure

| File | Description |
|------|-------------|
| `cali_main.m` | Main script that generates synthetic data (chessboard points, lines, mirror contour) and performs calibration |
| `coc_cali.m` | Core calibration function that integrates point, line and contour constraints |
| `fit_chessboard_line_image_fromcornerpoint.m` | Fits line images (conics) from detected corner points on each chessboard row/column |
| `plot_est_coc_image.m` | Visualizes the estimated projection results with mirror contour, principal point, and reprojected points |

> **Note:** Additional helper functions (e.g., `generate_cos_chessboard_image`, `generate_cos_world_point_image`, etc.) are used internally; they are expected to be present in the same directory.

### Dependencies

- MATLAB R2019b or later (tested with R2020a)
- No additional toolboxes are strictly required (uses built-in functions like `svd`, `eig`, etc.)

### Usage

1. Open MATLAB and navigate to the project folder.
2. Run the main script:
   ```matlab
   cali_main
   ```
   This will:
   - Generate synthetic checkerboard world points and line samples.
   - Simulate images under a central omnidirectional camera model (with configurable intrinsic parameters, pose, and distortion).
   - Call the calibration routine `coc_cali` to recover all parameters.
   - Plot the estimated images with reprojected points.

3. To calibrate your own data, modify the input arguments to `coc_cali`:
   - `world_point_image` : 3×N×M array of observed image points (homogeneous coordinates)
   - `world_point`        : 4×N array of 3D world points (homogeneous coordinates)
   - `line_image`         : 3×3×L×M array of line image conics (each 3×3 matrix)
   - `mirror_contour`     : 3×3 conic representing the mirror contour
   - `radial_distortion_order` : 0, 1, 2, or 3
   - `l`                  : optional (if known, provide a scalar; otherwise leave empty)

   See the function header in `coc_cali.m` for full details.

### Citation

If you use this code in your research, please cite the following work (to be added based on your reference).

### License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
