# Bridge for OntoAnno 2.0.0, ceda74c7d74af67a02093b07dd87623d846a59fe.
# One fresh R process handles at most one upstream batch. No clustering here.
args <- commandArgs(trailingOnly = TRUE)
pin <- "ceda74c7d74af67a02093b07dd87623d846a59fe"
if (!requireNamespace("OntoAnno", quietly = TRUE)) {
  stop("Install OntoAnno and its imports in this R library")
}
desc <- utils::packageDescription("OntoAnno")
if (!identical(as.character(desc$Version), "2.0.0") ||
    !pin %in% c(desc$RemoteSha, desc$GithubSHA1)) {
  stop("Install the pinned commit via remotes::install_github; RemoteSha mismatch")
}
required <- c("input", "tissue_name", "model", "topgenenumber", "add_cl_prompt",
              "restrict_to", "parent_celltype", "llm_config", "min_log2fc")
if (!all(required %in% names(formals(OntoAnno::gptcelltype)))) {
  stop("Incompatible OntoAnno annotation formals")
}
if (identical(args, "--probe")) {
  cat("OntoAnno", desc$Version, pin, "ellmer", as.character(packageVersion("ellmer")))
  quit(status = 0L)
}
if (length(args) != 2L) stop("Provide request.json and result.json")
job <- jsonlite::read_json(args[[1]], simplifyVector = FALSE)
settings <- job$settings
llm <- list(provider = settings$provider, model = settings$model,
            system_prompt = job$system_prompt)
if (isTRUE(settings$single_line_stop)) {
  # Only when the operator declares stop support; exactly one group is in
  # the batch. The provider ends generation at the label's line boundary.
  llm$params <- ellmer::params(stop_sequences = "\n")
}
if (!is.null(settings$api_key_env)) llm$api_key <- Sys.getenv(settings$api_key_env)
if (identical(settings$provider, "vllm")) {
  # Author call_llm uses ellmer::chat_vllm without api_key; ellmer instead
  # reads this variable. Keep it scoped to this private R process.
  Sys.setenv(VLLM_API_KEY = llm$api_key)
}
if (!is.null(settings$base_url)) llm$api_url <- settings$base_url
markers <- lapply(job$groups, function(group) unlist(group$genes, use.names = FALSE))
names(markers) <- vapply(job$groups, function(group) group$candidate_id, character(1))
stopifnot(length(markers) >= 1L, length(markers) <= 50L, !anyDuplicated(names(markers)))

# The public function returns parsed labels. Bind its transport in a private
# lexical environment so the adapter can save the unmodified provider reply.
# The upstream function body, prompt, parser, model transport and package
# namespace stay intact.
annotator <- OntoAnno::gptcelltype
audit <- new.env(parent = environment(annotator))
audit$prompt <- NULL
audit$reply <- NULL
audit$error <- NULL
audit$call_llm <- function(prompt, ...) {
  audit$prompt <- prompt
  value <- tryCatch(OntoAnno::call_llm(prompt = prompt, ...), error = function(e) {
    msg <- conditionMessage(e)
    http <- regmatches(msg, regexpr("HTTP [1-5][0-9]{2}", msg, perl = TRUE))
    audit$error <- if (length(http)) http else if (grepl("env var|API key", msg, ignore.case = TRUE)) {
      "credential configuration"
    } else if (grepl("timeout|timed out", msg, ignore.case = TRUE)) {
      "transport timeout"
    } else "upstream transport error"
    stop(e)
  })
  audit$reply <- as.character(value)
  value
}
environment(annotator) <- audit
ans <- annotator(input = markers, tissue_name = job$tissue, model = settings$model,
                 topgenenumber = as.integer(settings$top_genes),
                 add_cl_prompt = TRUE, restrict_to = NULL, parent_celltype = NULL,
                 llm_config = llm, min_log2fc = NULL)
if (is.null(audit$reply)) {
  reason <- if (is.null(audit$error)) "no upstream response" else audit$error
  jsonlite::write_json(list(reason = reason), paste0(args[[2]], ".failure.json"), auto_unbox = TRUE)
  stop("OntoAnno transport did not return a complete response")
}
parsed <- getFromNamespace(".parse_llm_annotation_lines", "OntoAnno")(audit$reply)
ok <- length(parsed) == length(markers) && setequal(names(ans), names(markers))
rows <- list()
if (ok) {
  for (id in names(markers)) {
    label <- trimws(unname(ans[[id]]))
    mapped <- OntoAnno::map_celltypes_to_cl(label, verbose = FALSE)
    rows[[length(rows) + 1L]] <- list(candidate_id = id, label = label,
                                      ontology_id = mapped$clid[[1]])
  }
}
result <- list(valid = ok, raw_response = audit$reply, prompt = audit$prompt,
               rows = rows, R = as.character(getRversion()),
               ellmer = as.character(packageVersion("ellmer")), commit = pin)
tmp <- paste0(args[[2]], ".tmp")
jsonlite::write_json(result, tmp, auto_unbox = TRUE, na = "null", null = "null")
if (!file.rename(tmp, args[[2]])) stop("Could not commit response")
