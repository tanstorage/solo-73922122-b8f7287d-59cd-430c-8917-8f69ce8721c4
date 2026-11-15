# 家系 SNP 传递核查服务（Pedigree Haplotype Transmission Audit）

在子代 SNP 进入相位分析前，核查一段染色体上「每个位点由父母各传递一个单倍型」的
方案是否可行，并求整段父、母单倍型编号切换次数之和的最小值。

- Python + FastAPI，JSON API，交互文档：<http://localhost:8000/docs>
- Docker Compose 一键启动

## 启动

```bash
docker compose up --build
```

然后打开 <http://localhost:8000/docs>，展开 `POST /api/v1/analyze`，点击
**Try it out** 提交请求；结果直接显示在 **Response** 中，也可用：

```bash
curl -sS -X POST http://localhost:8000/api/v1/analyze \
  -H 'Content-Type: application/json' -d @sample_request.json
```

健康检查：`GET /healthz`。

本地（不使用 Docker）运行：

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
pytest -q          # 需先 pip install -r requirements-dev.txt
```

## 请求

```json
{
  "chromosome": "chr1",
  "max_switches": 1,
  "sites": [
    {"position": 1000, "father": "0|1", "mother": "1|0", "child": "0/0"},
    {"position": 2000, "father": "0|1", "mother": "1|0", "child": "0/1"},
    {"position": 3000, "father": "0|1", "mother": "0|1", "child": "1/1"}
  ]
}
```

| 字段 | 约束 |
| --- | --- |
| `sites` | 2～40 个位点，`position` 必须严格递增 |
| `father` / `mother` | 已定相，形如 `0|1`、`1|0` |
| `child` | 未定相，形如 `0/1`、`1/1` |
| `max_switches` | 0～78 的整数（40 个位点、双亲各最多 39 次切换，故总和上限为 78） |
| `chromosome` | 可选标签，原样回显 |

违反格式、顺序或范围时返回 `422` 及具体错误信息。

## 判定规则

每个位点从父母各取一个等位；父母所取单倍型的编号（0 或 1）即「来源」。
相邻位点之间，父、母来源各自变化一次各计 1 次切换，求和得到整段代价。

- **min_switches**：所有可行传递方案中的最小切换次数。
- **确定（determined）**：存在可行方案，且**所有**最优方案在该位点对某亲本给出
  相同的来源编号。只有父、母两侧都确定，该位点才标记为 `determined`。
- **待判（ambiguous）**：该位点在最优方案集中，至少一个亲本的来源编号不唯一。
  注意：同一位点单独看可能两个来源都合法，但受全段最优性约束后仍可能是确定的。
- **冲突（conflict）**：某位点子代基因型无法由父母各一个等位组成，如 `0|0` ×
  `0|0` → `0/1`。此时 `feasible=false`，`min_switches=null`，冲突位点列在
  `conflict_positions`，其余位点标为 `unavailable`。
- **整段待复核（needs_review）**：`min_switches > max_switches` 时
  `status="needs_review"`、`over_limit=true`，并给出超出量 `limit_excess`。
  各最优方案的逐位点判定仍会返回，便于定位。

> 等位与来源是两件事：纯合亲本两个等位相同，但来源编号仍可能不唯一，此时
> `allele` 有值而 `source` 为 `null`，`status="ambiguous"`。

## 响应要点

```json
{
  "chromosome": "chr1",
  "site_count": 3,
  "max_switches": 1,
  "feasible": true,
  "status": "ok",
  "min_switches": 1,
  "over_limit": false,
  "limit_excess": null,
  "conflict_positions": [],
  "ambiguous_positions": [],
  "determined_positions": [1000, 2000, 3000],
  "sites": [
    {
      "index": 0,
      "position": 1000,
      "status": "determined",
      "father": {"status": "determined", "source": 0, "allele": 0,
                 "options": [{"source": 0, "allele": 0}]},
      "mother": {"status": "determined", "source": 1, "allele": 0,
                 "options": [{"source": 1, "allele": 0}]}
    }
  ],
  "example_assignment": {
    "father_switches": 1, "mother_switches": 0, "total_switches": 1,
    "sites": [{"index": 0, "position": 1000, "father_source": 0,
               "mother_source": 1, "father_allele": 0, "mother_allele": 0}]
  },
  "message": "Minimum switches 1 within the limit; every site has a unique parental source."
}
```

- `options` 列出所有最优方案中该亲本实际取过的 `(source, allele)` 组合。
- `example_assignment` 是一份具体的最优传递方案（并列最优时按来源编号最小选取），
  便于核对；待判位点也可能存在其他同样最优的来源。
- `message` 为可直接展示的中文/英文摘要（英文，避免编码歧义）。

## 算法

状态空间只有 4 个：`(父来源, 母来源) ∈ {0,1}²`。先由子代基因型筛出每位点的合法
状态，再用动态规划求最小切换，最后用前向/后向代价合成判断每个状态是否属于某个全局
最优解，从而逐位点、逐亲本给出「确定/待判」。复杂度 `O(位点数)`，40 个位点可忽略不计。

## 测试

```bash
pip install -r requirements-dev.txt
pytest -q
```

测试包含与「枚举所有传递方案」的暴力解随机对拍，覆盖并列最优导致的待判、纯合亲本
的来源不可判、冲突位点、超限、输入校验与 `/docs` 可达性。
