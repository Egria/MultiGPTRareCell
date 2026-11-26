"use strict";
const $ = (id) => document.getElementById(id);
const form = $("form");
const kinds = {
  pipeline: "完整发现流程",
  prepare: "准备快照",
  annotate: "注释快照",
  discover: "从快照发现",
  resume: "恢复发现",
  evaluate: "评估结果",
  demo: "离线 Core demo",
  inspect: "查看已有运行"
};
const states = {
  queued: "排队中",
  running: "运行中",
  stopping: "正在停止",
  succeeded: "已完成",
  partial: "部分完成",
  blocked: "受阻",
  failed: "失败",
  stopped: "已停止",
  interrupted: "已中断",
  paused: "已暂停",
  linked: "已关联"
};
const stages = {
  prepare: "准备",
  annotate: "注释",
  discover: "发现",
  evaluate: "评估",
  resume: "恢复",
  demo: "Core demo"
};
const metricNames = {
  n_cells: "细胞数",
  n_candidates: "候选数",
  n_groups: "群体数",
  unresolved_cells: "未决细胞",
  provider_dispatches: "发现 API 请求",
  reported_total_tokens: "已报告 token"
};
const timeText = (ts) =>
  ts
    ? new Date(ts * 1000).toLocaleString("zh-CN", {
        month: "2-digit",
        day: "2-digit",
        hour: "2-digit",
        minute: "2-digit",
        hour12: false
      })
    : "—";
let tasks = [],
  selected = null,
  current = null,
  busy = false,
  detailVersion = 0,
  loading = false;
const key = "rarecell-session";
if (location.hash.length > 1) {
  sessionStorage.setItem(key, location.hash.slice(1));
  history.replaceState(null, "", location.pathname);
}
const token = sessionStorage.getItem(key) || "";

async function api(path, body, download = false) {
  const options = { headers: { Authorization: "Bearer " + token } };
  if (body !== undefined) {
    options.method = "POST";
    options.headers["Content-Type"] = "application/json";
    options.body = JSON.stringify(body);
  }
  const res = await fetch(path, options);
  if (!res.ok) {
    const failure = await res.json().catch(() => ({}));
    throw new Error(failure.error || "请求失败（HTTP " + res.status + "）");
  }
  return download ? res.blob() : res.json();
}
function element(tag, text, cls) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  if (cls) node.className = cls;
  return node;
}
function message(id, text = "") {
  $(id).textContent = text;
  $(id).hidden = !text;
}
function badge(status) {
  return element("span", states[status] || status, "badge " + status);
}
function render() {
  const groups = {
    active: ["running", "stopping"],
    done: ["succeeded"],
    attention: ["partial", "blocked", "failed", "interrupted", "paused"]
  };
  $("total").textContent = tasks.length;
  for (const k of ["queued", "active", "done"])
    $(k).textContent = tasks.filter((t) => (groups[k] || [k]).includes(t.status)).length;
  const filter = $("filter").value,
    query = $("search").value.trim().toLowerCase();
  const rows = tasks.filter(
    (t) =>
      (filter === "all" || (groups[filter] || [filter]).includes(t.status)) &&
      t.name.toLowerCase().includes(query)
  );
  const fragment = document.createDocumentFragment();
  for (const task of rows) {
    const row = element("tr"),
      cell = element("td"),
      link = element("button", task.name, "task-link");
    link.type = "button";
    link.onclick = () => openDetail(task.id);
    cell.append(link, element("div", kinds[task.kind], "sub"));
    row.append(cell);
    const mode = element("td");
    mode.append(task.mode ? element("span", task.mode, "tag") : "—");
    row.append(mode, element("td", stages[task.stage] || "—"));
    const state = element("td");
    state.append(badge(task.status));
    row.append(state, element("td", timeText(task.created)));
    fragment.append(row);
  }
  const focused = document.activeElement;
  const focusName =
    focused && focused.classList.contains("task-link")
      ? [...$("rows").querySelectorAll("button")].indexOf(focused)
      : -1;
  $("rows").replaceChildren(fragment);
  if (focusName >= 0) $("rows").querySelectorAll("button")[focusName]?.focus({ preventScroll: true });
  $("empty").hidden = rows.length > 0;
  $("empty-title").textContent = tasks.length ? "没有匹配的任务" : "队列中还没有任务";
  $("empty-detail").textContent = tasks.length
    ? "可以调整搜索词或状态筛选。"
    : "从数据集开始一次发现，或选择已有快照继续分析。";
  $("empty-action").textContent = tasks.length ? "清除筛选" : "运行离线示例";
  $("count").textContent =
    rows.length === tasks.length
      ? "共 " + tasks.length + " 个任务"
      : "显示 " + rows.length + " / " + tasks.length + " 个任务";
}
let savedModel = null;
function displayModel(value) {
  const input = $("project-model");
  if (savedModel === null || input.value === savedModel) input.value = value;
  savedModel = value;
  input.disabled = false;
  $("save-model").disabled = false;
  const hint = value ? "全局模型：" + value : "沿用配置文件中的模型";
  for (const name of ["model", "annotation_model"]) form.elements[name].placeholder = hint;
}
$("project-settings").onsubmit = async (event) => {
  event.preventDefault();
  const input = $("project-model");
  const model = input.value.trim();
  try {
    await api("/api/model", { model });
    input.value = model;
    displayModel(model);
    $("model-status").textContent = model ? "已保存" : "已清除";
  } catch (err) {
    $("model-status").textContent = err.message;
  }
};
$("project-model").oninput = () => ($("model-status").textContent = "");
async function refresh() {
  if (loading) return;
  loading = true;
  try {
    const data = await api("/api/tasks");
    if (data.ui_revision && document.documentElement.dataset.uiRevision !== data.ui_revision &&
        !$("create").open && !busy && (savedModel === null || $("project-model").value === savedModel)) {
      sessionStorage.setItem("rarecell-reopen", JSON.stringify({id:selected, diagnostics:$("diagnostics").open,
        filter:$("filter").value, search:$("search").value}));
      location.reload();
      return;
    }
    tasks = data.tasks;
    displayModel(data.model);
    render();
    $("workspace").textContent = data.workspace;
    $("connection").lastElementChild.textContent = data.host;
    $("connection").classList.add("online");
    $("new").disabled = false;
    message("notice");
    if ($("detail").open && selected) await renderDetail(selected);
  } catch (err) {
    message("notice", err.message);
    $("connection").classList.remove("online");
    $("connection").lastElementChild.textContent = "连接中断";
  } finally {
    loading = false;
  }
}
function field(id, show, required = false) {
  const box = $(id);
  box.hidden = !show;
  for (const input of box.querySelectorAll("input,select")) {
    input.disabled = !show;
    input.required = show && required;
  }
}
function updateForm() {
  const kind = $("kind").value,
    mode = $("mode").value,
    full = kind === "pipeline";
  for (const el of document.querySelectorAll("[data-kinds]")) {
    const show = el.dataset.kinds.split(" ").includes(kind);
    el.hidden = !show;
    el.querySelectorAll("input,select").forEach((n) => (n.disabled = !show));
  }
  const needsDataset = kind === "prepare" || (full && $("start").value === "dataset");
  field("dataset-field", needsDataset, true);
  field("imported-field", needsDataset && $("imported-view").options.length > 1);
  field(
    "snapshot-field",
    ["annotate", "discover"].includes(kind) || (full && $("start").value === "snapshot"),
    true
  );
  const annotated = ["A1", "A3", "A3_adaptive", "A3_routed"].includes(mode);
  field("routed-workers-field", ((full || kind === "discover") && mode === "A3_routed") || kind === "resume");
  field("seed-review-field", (full || kind === "discover") && ["A2", "A3", "A3_adaptive"].includes(mode));
  const showAnnotations = (full || kind === "discover") && annotated;
  field("annotations-field", showAnnotations, kind === "discover");
  $("annotations-hint").textContent = full ? "可选，留空自动注释" : "必填";
  const annModel = kind === "annotate" || (full && annotated && !form.elements.annotations.value.trim());
  const model = full || kind === "discover";
  $("models").hidden = !model && !annModel;
  field("model-field", model);
  field("ann-model-field", annModel);
  field("singleton-annotation-field", annModel);
  field("panel-mode-field", annModel);
  field("reuse-annotation-field", annModel && (kind === "annotate" || $("start").value === "snapshot"));
  form.elements.gt.required = kind === "evaluate";
  $("gt-hint").textContent = kind === "evaluate" ? "必填" : "可选，结束后评估";
  const hints = {
    inspect: "填写含 checkpoint.json 的发现运行目录，只读查看已有记录。",
    demo: "使用内置合成数据运行并评估，不调用模型 API。",
    pipeline: "模型从服务端环境读取密钥。相对路径以启动控制台的目录为准。",
    prepare: "使用配置中指定的视角生成快照，结果保存在本任务目录。",
    annotate: "将调用配置中的 OntoAnno 模型，可在后续发现任务中复用输出目录。",
    discover: "使用真实主 Agent 调查候选，API 调用使用服务端环境中的密钥。",
    resume: "恢复原发现进度。A3_routed 因预算停止时，请提高已耗尽的总上限。",
    evaluate: "从冻结结果和完整的 cell_id、label 真值表生成评估。"
  };
  $("form-hint").textContent = hints[kind];
  $("enqueue").textContent = kind === "inspect" ? "查看诊断" : "加入队列";
}
const imported = $("imported-view");
let autoGt = "";
function clearAutoGt() {
  if (autoGt && form.elements.gt.value === autoGt) form.elements.gt.value = "";
  autoGt = "";
}
imported.onchange = () => {
  clearAutoGt();
  if (imported.value) {
    form.elements.dataset.value = imported.selectedOptions[0].dataset.manifest;
    autoGt = imported.selectedOptions[0].dataset.gt;
    form.elements.gt.value = autoGt;
  }
};
form.elements.dataset.addEventListener("input", () => {
  imported.value = "";
  clearAutoGt();
});
$("start").addEventListener("change", () => {
  if ($("start").value === "snapshot") {
    imported.value = "";
    clearAutoGt();
  }
});
async function loadImported() {
  try {
    const views = await api("/api/datasets");
    const options = [new Option("选择已导入的 view，或手动填写路径", "")];
    for (const view of views) {
      const title = view.profile || (view.gt ? "含评估 GT" : "仅表达");
      const option = new Option(view.label + " · " + title, view.profile ? view.gt : view.dataset);
      option.dataset.manifest = view.dataset;
      option.dataset.gt = view.gt;
      options.push(option);
    }
    imported.replaceChildren(...options);
    updateForm();
  } catch (err) {
    message("form-error", err.message);
  }
}
function openForm(kind = "pipeline", preset = {}) {
  form.reset();
  imported.replaceChildren();
  autoGt = "";
  loadImported();
  $("kind").value = kind;
  for (const [k, v] of Object.entries(preset)) if (v && form.elements[k]) form.elements[k].value = v;
  message("form-error");
  updateForm();
  $("create").showModal();
}
async function renderDetail(id) {
  const seq = ++detailVersion;
  try {
    const t = await api("/api/tasks/" + id);
    if (selected !== id || seq !== detailVersion || !$("detail").open) return;
    current = t;
    $("open-diagnostics").disabled = false;
    $("detail-title").textContent = t.name;
    $("detail-type").textContent = kinds[t.kind] + (t.mode ? " / " + t.mode : "");
    $("detail-status").replaceChildren(badge(t.status));
    const facts = [
      ["当前阶段", stages[t.stage]],
      ["创建", timeText(t.created)],
      ["开始", t.started && timeText(t.started)],
      ["结束", t.ended && timeText(t.ended)],
      ["发现模型", t.models.provider],
      ["注释模型", t.models.annotation],
      ["任务目录", t.output],
      ["运行目录", t.run_path],
      ["状态详情", t.detail],
      ...Object.entries(t.metrics).map(([k, v]) => [metricNames[k], String(v)])
    ];
    const fragment = document.createDocumentFragment();
    for (const [k, v] of facts) if (v) fragment.append(element("dt", k), element("dd", v));
    $("facts").replaceChildren(fragment);
    const outputs = $("outputs");
    outputs.replaceChildren();
    if (t.files.length) outputs.append(element("h3", "结果文件"));
    for (const name of t.files) {
      const button = element("button", name, "file secondary");
      button.onclick = async () => {
        try {
          const data = await api("/api/tasks/" + id + "/files/" + name, undefined, true);
          const url = URL.createObjectURL(data),
            a = element("a");
          a.href = url;
          a.download = name;
          a.click();
          setTimeout(() => URL.revokeObjectURL(url), 30000);
        } catch (err) {
          message("detail-error", err.message);
        }
      };
      outputs.append(button);
    }
    $("events").hidden = !t.events.length;
    $("event-log").textContent = t.events.join("\n");
    const log = $("log"),
      atEnd = log.scrollHeight - log.scrollTop - log.clientHeight < 12;
    const text = t.log || (t.status === "queued" ? "任务启动后显示执行日志。" : "暂无日志。");
    if (log.textContent !== text) {
      log.textContent = text;
      if (atEnd) log.scrollTop = log.scrollHeight;
    }
    $("stop").hidden = !["queued", "running"].includes(t.status);
    $("stop").textContent = t.status === "queued" ? "取消排队" : "停止任务";
    $("resume").hidden = !t.resumable;
  } catch (err) {
    if (selected === id) message("detail-error", err.message);
  }
}
function openDetail(id) {
  selected = id;
  current = null;
  $("open-diagnostics").disabled = true;
  message("detail-error");
  $("detail-title").textContent = "任务详情";
  $("detail-type").textContent = "";
  for (const id of ["facts", "detail-status", "outputs", "log"]) $(id).replaceChildren();
  $("events").hidden = $("stop").hidden = $("resume").hidden = true;
  $("detail").showModal();
  renderDetail(id);
}
$("form").onsubmit = async (e) => {
  e.preventDefault();
  if (busy) return;
  busy = true;
  $("enqueue").disabled = true;
  try {
    const result = await api("/api/tasks", Object.fromEntries(new FormData(form)));
    $("create").close();
    $("filter").value = "all";
    $("search").value = "";
    await refresh();
    openDetail(result.id);
    if ($("kind").value === "inspect") openDiagnostics(result.id, form.elements.name.value || "已有运行");
  } catch (err) {
    message("form-error", err.message);
  } finally {
    busy = false;
    $("enqueue").disabled = false;
  }
};
$("stop").onclick = async () => {
  if (!current) return;
  if (!window.confirm("确认停止此任务？如果已经发送 API 请求，该请求可能仍会计费。")) return;
  try {
    await api("/api/tasks/" + selected + "/stop", {});
    await refresh();
  } catch (err) {
    message("detail-error", err.message);
  }
};
$("resume").onclick = () => {
  if (!current) return;
  const preset = { name: current.name + " · 恢复", run: current.run_path, gt: current.spec.gt };
  $("detail").close();
  openForm("resume", preset);
};
$("detail").addEventListener("close", () => {
  selected = null;
  current = null;
});
$("new").onclick = () => openForm();
$("empty-action").onclick = () => {
  if (!tasks.length) openForm("demo");
  else {
    $("search").value = "";
    $("filter").value = "all";
    render();
  }
};
for (const id of ["mode", "kind", "start"]) $(id).onchange = updateForm;
form.elements.annotations.oninput = updateForm;
$("search").oninput = render;
$("filter").onchange = render;
document.querySelectorAll("[data-close]").forEach((b) => (b.onclick = () => $(b.dataset.close).close()));
window.addEventListener("DOMContentLoaded", () => {
  const saved = sessionStorage.getItem("rarecell-reopen");
  sessionStorage.removeItem("rarecell-reopen");
  if (saved) {
    const view = JSON.parse(saved);
    $("filter").value = view.filter;
    $("search").value = view.search;
    if (view.id) {
      openDetail(view.id);
      if (view.diagnostics) openDiagnostics(view.id, view.id);
    }
  }
  refresh();
  setInterval(refresh, 2500);
});
