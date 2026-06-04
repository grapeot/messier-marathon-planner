# Skill：梅西耶马拉松观测规划

## 元数据

- **类型**：Workflow
- **适用场景**：用户想规划某个地点、某个日期、某套设备下的梅西耶天体拍摄/观测顺序
- **项目仓库**：https://github.com/grapeot/messier-marathon-planner
- **默认输出**：Markdown 报告、JSON 排名、PNG 可视化图表

## 目标

为指定地点、日期和设备，生成一份可查证的梅西耶马拉松规划报告。报告必须同时回答三个问题：今晚能拍哪些、为什么推荐这些、读者怎样用公开工具核对计算结果。

## 验收标准

一次规划任务完成时，至少满足以下条件：

- 报告列出目标地点、日期、设备参数、太阳暗夜窗口和可见性阈值
- 对每个可见梅西耶深空目标给出：编号、名称、类型、视星等、角大小、平均面亮度、暗夜窗口最高时间、高度、可见时长、建议曝光、综合评分
- Top 5 推荐目标附带公开交叉验证入口，至少包含 Stellarium Web 的使用说明和一个天体信息网页链接
- 报告包含太阳 + Top10 推荐目标高度随时间变化图，横轴为本地时间，纵轴为高度角，太阳 -12° 暗夜线清楚可见
- 运行 `pytest` 通过，报告和图表可重新生成
- 若发布到网页，正文必须包含源代码仓库链接，方便读者复现

## 可用资源

项目仓库提供以下核心文件：

- `src/messier_catalog.py`：梅西耶星表数据
- `src/planner.py`：主规划程序，负责可见性、评分、报告和主验证图
- `src/verify_viz.py`：额外验证图，包括 RA-Dec 天图、高度热力图、Pareto 图、FOV 对比和时间线
- `tests/`：星表、FOV、面亮度、曝光估算、可见性集成测试
- `docs/rfc.md`：公式、权重和设计边界

## 方法论建议

优先把结果做成可查证，而不是只给一份漂亮排名。推荐目标时不要只看视星等；视场匹配、平均面亮度、暗夜窗口高度和目标类型都要同时考虑。

对于 Seestar 这类智能望远镜，暗夜窗口可以使用太阳高度 < -12° 的航海晨昏阈值。传统目视观测或严肃深空摄影可以改成 -18°，但报告中必须明确写出阈值。

平均面亮度用于排序，而不是精确曝光预测。星系核心、球状星团核心和发射星云亮区都可能比平均值亮得多。曝光时间应作为保守建议，并标注为累计曝光时间。

## 已知陷阱

- **赤纬负号容易错**：南天天体如 M6、M7、M16 等必须正确处理负赤纬。M6 在北京最高只有约 18°，如果被算到 80° 以上，说明 Dec 符号解析错了。
- **M40/M73/M24 不适合作为常规推荐目标**：M40 是双星，M73 是星群，M24 是人马座恒星云。它们可以保留在星表里，但推荐排序应默认过滤。
- **物理事实不是交叉验证**：距离、年龄、构成等事实可以有趣，但不能验证今晚北京的可见性。可见性验证必须给 Stellarium Web、In-The-Sky、Fourmilab Your Sky 等公开工具入口。
- **“最高点”要限定窗口**：报告里写的是暗夜可拍摄窗口内最高点，不一定是天体全天的中天高度。

## 输出规格

推荐输出目录：`output/`。

必须生成：

- `output/report.md`
- `output/results.json`
- `output/top10_sun_altitude.png`

建议生成：

- `output/altitude_curves.png`
- `output/scoring_bars.png`
- `output/verify_sky_map.png`
- `output/verify_altitude_heatmap.png`
- `output/verify_pareto.png`
- `output/verify_fov_comparison.png`
- `output/verify_timeline.png`

## 安装方式

把仓库链接交给 AI coding agent，并说明：

> 请安装并使用这个 skill 帮我规划梅西耶马拉松：https://github.com/grapeot/messier-marathon-planner

Agent 应从目标 workspace 的 `AGENTS.md` / `CLAUDE.md` / skills index 入手，把本文件作为一个可发现的 workflow skill。若目标 workspace 有 `rules/skills/INDEX.md` 或 `skills/INDEX.md`，应添加一条指向本文件的入口；如果没有 skills 目录，则在项目级 `AGENTS.md` 中加入指针。

## 最小运行命令

```bash
uv venv
uv pip install astropy numpy matplotlib pytest
.venv/bin/python src/planner.py
.venv/bin/python src/verify_viz.py
.venv/bin/python -m pytest tests/ -v
```
