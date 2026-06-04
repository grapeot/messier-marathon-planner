# 梅西耶马拉松观测规划工具

给定观测地点、日期和设备，自动规划梅西耶天体的拍摄方案：计算可见性、面亮度、综合评分和推荐曝光时间，并输出可视化报告。

## 功能

- **110 个梅西耶天体**的完整星表（SEDS 数据源）
- 基于 **astropy** 的可见性计算（高度角 + 太阳位置）
- **面亮度计算**和曝光时间估算（口径归一化 + 目标类型校正）
- **综合评分**：高度(30%)、亮度(35%)、可见时长(10%)、FOV 匹配(25%)
- 支持 **ZWO SeeStar S50** 和 **S30 Pro** 双设备
- 6 张验证可视化图表：太阳+Top10 高度图、天图、时间热力图、Pareto 前沿、FOV 对比、时间线
- 报告附录提供 **公开验证链接**（Stellarium Web / TheSkyLive / In-The-Sky）

## 快速开始

```bash
# 1. 创建虚拟环境
uv venv
uv pip install astropy numpy matplotlib pytest

# 2. 运行规划
.venv/bin/python src/planner.py

# 3. 生成额外验证图
.venv/bin/python src/verify_viz.py

# 4. 查看报告
open output/report.md

# 5. 运行测试
.venv/bin/python -m pytest tests/ -v
```

## 最重要的验证图

`output/top10_sun_altitude.png` 是最适合人工查证的图：横轴是北京时间，纵轴是高度角，红色虚线是太阳，灰色虚线是太阳 -12° 暗夜线。Top 10 推荐目标的高度曲线会叠在同一张图上，曲线圆点表示暗夜窗口内的最高位置。

这张图能直接回答三个问题：

- 推荐目标是否真的在暗夜窗口内足够高
- 目标之间的最佳拍摄时间是否冲突
- 太阳高度和目标高度的组合是否支持报告中的排序

## 设备参数

| 参数 | S50 | S30 Pro |
|------|-----|---------|
| 口径 | 50mm | 30mm |
| 焦比 | f/5 | f/5 |
| 焦距 | 250mm | 150mm |
| 传感器 | IMX462 (5.57×3.13mm) | IMX585 (11.14×6.26mm) |
| 视场 (H×V) | 1.28°×0.72° | 4.25°×2.39° |
| 像素分辨率 | 2.39"/px | 3.99"/px |

## 数据来源

- 梅西耶星表：SEDS Messier Catalog
- 设备参数：ZWO 官方网站
- 天球计算：astropy (AltAz, get_sun)

## 许可

MIT License
