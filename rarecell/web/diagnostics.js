"use strict";

const diagDialog = $("diagnostics");
let diagTask = null, diagPage = 0, diagVersion = 0;
let diagnosticText = "";

function diagFacts(entries) {
  const list = element("dl");
  for (const [label, value] of entries) {
    if (value !== undefined && value !== null && value !== "") {
      list.append(element("dt", label), element("dd", String(value)));
    }
  }
  return list;
}

function diagTable(headers, rows) {
  const wrapper = element("div", undefined, "diag-table");
  const table = element("table"), head = table.createTHead().insertRow(), body = table.createTBody();
  for (const label of headers) head.append(element("th", label));
  for (const values of rows) {
    const row = body.insertRow();
    for (const value of values) row.append(element("td", value ?? "—"));
  }
  wrapper.append(table);
  return wrapper;
}

function disclosure(label, content) {
  const box = element("details", undefined, "diag-disclosure");
  box.append(element("summary", label), content);
  return box;
}

function showOverview(data) {
  const box = $("diag-overview"), s = data.summary, a = data.annotation, u = data.usage, scores = data.evaluation;
  box.replaceChildren();
  if (s.detail) {
    const title = s.detail === "stalled" ? "发现已停滞（stalled）" : "发现结果 · " + s.detail;
    const note = element("div", undefined, s.submitted ? "diag-outcome" : "diag-outcome attention");
    note.append(element("strong", title));
    if (s.detail === "stalled") note.append(element("p", "连续多轮动作都没有被控制器计为新进展，运行已停止。"));
    if (s.issue) note.append(element("p", s.issue));
    box.append(note);
  }
  box.append(diagFacts([
    ["准备 / 注释", [data.prepare.status, a.status || s.annotation_status].map(v => v ?? "暂无记录").join(" / ")],
    ["注释详情", a.detail], ["已提交分群", s.submitted == null ? null : s.submitted ? "是" : "否"],
    ["未决细胞 / 总数", s.n_cells == null ? null : (s.unresolved_cells ?? "—") + " / " + s.n_cells],
    ["待调查种子", s.n_pending_seeds], ["回退规则", s.fallback_view_resolution?.join(" · 分辨率 ")],
    ["连续无进展 / 停机阈值", data.trace.stall != null ?
      data.trace.stall + " / " + (data.stall_limit ?? "—") : null]
  ]));
  const usage = diagFacts([
    ["注释模型", data.annotation_model], ["已注释分组", a.annotated],
    ["新完成 / 缓存注释批次", a.newly_completed_batches != null ? a.newly_completed_batches + " / " + (a.cached_batches ?? "—") : null],
    ["注释 API 次数 / 费用（USD）", a.annotated != null ?
      (a.provider_calls_observed_at_wire ?? "未采集") + " / " + (a.annotation_cost_usd ?? "未采集") : null],
    ["发现轮数", u.turns], ["发现 API 请求", u.provider_dispatches], ["新科学查询", u.new_queries],
    ["协议错误", u.protocol_errors], ["动作错误", u.action_errors]
  ]);
  if (usage.childElementCount) box.append(disclosure("模型调用与错误计数", usage));
  if (data.pending_seeds.length) box.append(disclosure("待调查种子 · 显示 " + data.pending_seeds.length + " 个",
    diagTable(["候选 ID", "细胞 ID"], data.pending_seeds.map(row => [row.candidate_id, row.cell_id]))));
  if (data.prepare.stages) box.append(disclosure("准备阶段明细", diagTable(["阶段", "状态", "结果 / 原因"],
    Object.entries(data.prepare.stages).map(([name, row]) => [name, row.status, row.reason ?? row.shape?.join(" × ") ??
      (row.informative != null ? row.informative + " 个有效细胞" : "")]))));
  if (a.errors?.length || a.skipped?.length) box.append(disclosure("注释错误与跳过的分组",
    element("pre", JSON.stringify({errors:a.errors, skipped:a.skipped}, null, 2))));
  if (scores) {
    const source = {fallback:"全部来自回退分群", mixed:"已选候选与回退分群的组合", claims:"全部来自已选候选", unknown:"分群来源未记录"};
    const base = scores.frozen.fallback_view_resolution;
    const lines = [
      ["分群来源", source[scores.source] + (base && ["fallback","mixed"].includes(scores.source) ? "（" + base.join(" / ") + "）" : "")],
      ["评估时已提交", scores.frozen.submitted == null ? null : scores.frozen.submitted ? "是" : "否"],
      ["ARI / NMI / Macro F1", [scores.ari,scores.nmi,scores.macro_f1].map(v=>v?.toFixed(3) ?? "—").join(" / ")]
    ];
    const targets = diagTable(["目标","细胞数","TP","FP","FN","F1"],
      scores.targets.map(r => [r.target,r.n_target,r.tp,r.fp,r.fn,r.f1?.toFixed(3)]));
    targets.classList.add("target-scores");
    box.append(element("h3", "评估使用的冻结分群"), diagFacts(lines), disclosure("各目标 F1", targets));
  } else if (s.detail) {
    box.append(element("p", "暂无与本任务关联的评估结果。", "sub"));
  }
}

function turnDisclosure(id, turn) {
  const body = element("div"), number = "T" + String(turn.number).padStart(4, "0");
  const box = disclosure(number + " · 原始回复、请求参数与完整回执", body);
  let started = false;
  box.ontoggle = async () => {
    if (!box.open || started) return;
    started = true;
    body.textContent = "读取中…";
    try {
      const data = await api("/api/tasks/" + id + "/diagnostics/turn/" + turn.number);
      if (!box.isConnected || diagTask !== id) return;
      body.replaceChildren(element("h4", "模型原始回复"), element("pre", data.message.raw || "尚无可用的回复正文"));
      if (data.message.complete === false) body.append(element("p","传输记录标记为未完成。"));
      for (const [label, value] of [["已解析的请求", data.commands], ["动作回执", data.receipts], ["服务端调用错误", data.failures]]) {
        if (value.length) body.append(element("h4", label),element("pre", JSON.stringify(value, null, 2)));
      }
    } catch (err) {
      body.textContent = err.message;
      started = false;
    }
  };
  return box;
}

function showTrace(id, trace) {
  const box = $("diag-turns");
  box.replaceChildren();
  for (const turn of trace.turns || []) {
    const section = element("section", undefined, "diag-turn");
    const label = turn.progress == null ? "尚无轮次结束记录" : turn.progress ? "有进展" : "无新进展";
    section.append(element("h4", "T" + String(turn.number).padStart(4, "0") + " · " + label,
      turn.progress === false ? "no-progress" : ""));
    if (turn.actions.length) section.append(diagTable(["动作","回执状态","进展","原因"],
      turn.actions.map(a => [a.action ?? "未解析", a.status, a.changed ? "是" : "否", a.error ||
        (a.status === "no_expansion" ? "没有新增细胞" : a.changed === false ? "控制器未计入新进展" : "")])));
    section.append(turnDisclosure(id, turn));
    box.append(section);
  }
  if (!trace.turns?.length) box.append(element("p", trace.notice || "暂无可读取的发现轮次。"));
  diagPage = trace.page || 0;
  $("diag-pager").hidden = !(trace.pages > 1);
  $("diag-page").textContent = diagPage + " / " + (trace.pages || 1);
  $("diag-prev").disabled = diagPage <= 1;
  $("diag-next").disabled = diagPage >= trace.pages;
}

async function loadDiagnostics(page) {
  const id = diagTask, version = ++diagVersion;
  $("diag-status").textContent = "正在读取运行记录…";
  try {
    const data = await api("/api/tasks/" + id + "/diagnostics?page=" + page);
    if (version !== diagVersion || !diagDialog.open) return;
    showOverview(data);
    showTrace(id, data.trace);
    const lines = ["任务 " + id, "发现 " + (data.summary.detail ?? "暂无结束记录") +
      " · submitted=" + (data.summary.submitted ?? "未知")];
    for (const t of data.trace.turns || []) {
      const label = "T" + String(t.number).padStart(4, "0");
      if (!t.actions.length) lines.push(label + " 暂无动作回执");
      for (const a of t.actions)
        lines.push([label, a.action ?? "未解析", a.status, "progress=" + a.changed, a.error].filter(Boolean).join(" "));
    }
    diagnosticText = lines.join("\n");
    $("diag-copy").disabled = false;
    $("diag-status").textContent = "已读取 · " + new Date().toLocaleTimeString("zh-CN");
  } catch (err) {
    if (version === diagVersion) $("diag-status").textContent = err.message;
  }
}
function openDiagnostics(id, name) {
  diagTask = id;
  diagPage = 0;
  diagnosticText = "";
  $("diag-copy").disabled = true;
  $("diagnostic-task").textContent = name;
  $("diag-overview").replaceChildren();
  $("diag-turns").replaceChildren();
  $("diag-pager").hidden = true;
  diagDialog.showModal();
  loadDiagnostics(0);
}
$("open-diagnostics").onclick = () => openDiagnostics(current.id, current.name);
$("diag-refresh").onclick = () => loadDiagnostics(0);
$("diag-copy").onclick = async () => {
  try {
    await navigator.clipboard.writeText(diagnosticText);
    $("diag-status").textContent = "已复制本页诊断";
  } catch {
    const text = element("textarea");
    text.value = diagnosticText;
    text.readOnly = true;
    text.setAttribute("aria-label", "本页诊断文本");
    $("diag-status").replaceChildren("可手动复制所选文本。", text);
    text.select();
  }
};
$("diag-prev").onclick = () => loadDiagnostics(diagPage - 1);
$("diag-next").onclick = () => loadDiagnostics(diagPage + 1);
$("close-diagnostics").onclick = () => diagDialog.close();
diagDialog.onclose = () => { diagTask = null; diagVersion++; };
