# 工作日志

## Changelog

### 2026-06-04

- 创建项目骨架（src/, docs/, tests/, output/）
- 编写 messier_catalog.py：从 SEDS 数据构建完整 110 天体星表
- 编写 planner.py：基于 astropy 的可见性计算、面亮度、综合评分、曝光估算
- 编写 verify_viz.py：天图、高度热力图、Pareto 前沿、FOV 对比、观测时间线
- Sub-agent 独立 fact check 并修复以下问题：
  - Dec 解析 bug：负号处理不正确导致 M6 (Dec=-32°) 被当作北天目标
  - M73 类型更正为 AST，加入过滤列表
  - Top5 验证文案改为动态生成，消除硬编码错误
- 补充完整项目文档（prd/rfc/test）
- 初始化独立 git repo，设置 .gitignore
- 创建 GitHub public repo 并推送
- 新增太阳 + Top10 推荐目标高度图，作为报告里的主验证图，便于人工核对目标是否在暗夜窗口内足够高
- 新增 `skills/messier_marathon_planner.md`，把项目沉淀为可安装的 AI workflow skill
- 在报告正文补充 GitHub 仓库和 skill 使用说明，方便读者直接复用到其他城市、日期和设备

## Lessons Learned

### Dec 符号解析

`dec_sign` 字段需要同时考虑 `dec_d` 本身为负的情况。原始代码只检查 `dec_sign`，但 M6 的 Dec = -32°13' 在数据录入时 `dec_d=32, dec_sign=+1`（对北天的写法）。修复方案：在 `ra_dec_to_deg` 中同时检查 `dec_d < 0` 和 `dec_sign < 0`。

**教训**：星表数据录入时应统一使用 `dec_sign` 明确符号，避免依赖 `dec_d` 的符号（有些来源对南天天体使用绝对值 + 负号标志）。

### 评分系统对低 Dec 目标的处理

M6 是评分较高的目标（视星等 4.2），但在北京（N39.9°）的 Dec=-32° 意味着它最高只能到约 18°。原来的评分没有考虑这种极限情况。修复：通过正确的可见性计算（astropy AltAz transform）自动处理——M6 的 max_alt 正确计算为 18°，排名随之下降。

### 物理事实 vs 交叉验证

第一版报告用了"物理事实"（距离、年龄、构成），但用户实际要的是可公开验证的数据对照。最终版改为提供 Stellarium Web / TheSkyLive / In-The-Sky 的直接链接，让读者可以独立核对报告中的高度、位置计算。
