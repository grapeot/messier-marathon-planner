# 测试策略

## 测试范围

- **单元测试**：星表数据完整性、面亮度公式、FOV 计算、曝光公式
- **集成测试**：可见性计算的端到端验证（依赖 astropy）
- **不测 E2E**：报告生成和可视化图表依赖 matplotlib，输出为文件，人工验证即可

## 测试命令

```bash
cd adhoc_jobs/messier_marathon
.venv/bin/python -m pytest tests/ -v
```

## 验证清单

| 检查项 | 方法 | 预期结果 |
|--------|------|----------|
| 110 个梅西耶天体全部有数据 | `test_catalog_completeness` | 110 条记录 |
| M31 坐标接近 0h42m +41° | `test_catalog_spot_check` | RA≈10.7°, Dec≈41.3° |
| 面亮度公式数学正确 | `test_surface_brightness` | M42 面亮度约 19-20 |
| FOV 计算正确 | `test_fov_calculation` | S50 FOV ≈ 1.28°×0.72° |
| 可见性返回合理数量 | `test_visibility_count` | 50-100 个目标有窗口 |
| Sun altitude 计算有日出日落 | `test_sun_altitude` | 北京 6 月有日出和日落 |

## Skill 文档人工检查

每次修改 `skills/messier_marathon_planner.md` 后，人工检查：

- 是否有明确验收标准
- 是否说明输出文件和验证方式
- 是否避免写成死板 SOP
- 是否不含私有路径、邮箱、token、服务器地址
