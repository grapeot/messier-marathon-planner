# AGENTS.md — 梅西耶马拉松观测规划工具

## 项目结构

```
src/messier_catalog.py   # 梅西耶星表（110天体）
src/planner.py           # 主程序：可见性、评分、报告、图表
src/verify_viz.py        # 验证可视化
tests/                   # pytest 测试
docs/                    # prd / rfc / working / test
output/                  # 生成产物（gitignored）
```

## 环境

- Python 3.12+，使用 `uv` 管理 venv
- 创建环境：`uv venv && uv pip install astropy numpy matplotlib pytest`
- 运行：`.venv/bin/python src/planner.py`
- 测试：`.venv/bin/python -m pytest tests/ -v`

## 规则

- 修改代码后先跑测试，再重新生成报告
- 星表数据变更后手动验证 M31、M42、M57 三颗标准星的坐标
- 每次修改后更新 `docs/working.md` 的 Changelog
- 生成产物在 `output/` 目录，不提交到 git
