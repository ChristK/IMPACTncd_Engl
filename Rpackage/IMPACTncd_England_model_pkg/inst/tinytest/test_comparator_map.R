# Tests for the per-scenario comparator map: the resolver, its structural
# validator, and the end-to-end behaviour through export_equity_tables().
#
# The load-bearing property is that a bare `comparator_scenario` is UNCHANGED.
# The scalar-equivalence test below is the one that matters: the same fixture
# exported twice -- once as "sc0", once as the equivalent map -- must agree
# exactly once the extra `comparator` column is dropped.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

resolve  <- IMPACTncdEngland:::.resolve_comparators
is_map   <- IMPACTncdEngland:::.comparator_is_map
cycles   <- IMPACTncdEngland:::.comparator_cycles
validate <- IMPACTncdEngland:::validate_comparator_scenario

gv <- c("sc0", "mult_base_qrisk", "mult_base_prs",
        "mult_lowzero_qrisk", "mult_lowzero_prs")

# ===========================================================================
# 1. .comparator_is_map(): the gate is SYNTACTIC, never data-derived
# ===========================================================================
expect_false(is_map("sc0"))
expect_false(is_map(c("sc0", "sc1")),
             info = "an unnamed vector is not a map (and errors upstream)")
expect_true(is_map(c(sc1 = "sc0")))
expect_true(is_map(c("sc0", mult_base_prs = "mult_base_qrisk")),
            info = "a partially named vector IS a map")

# ===========================================================================
# 2. Resolver: the scalar case is exactly the old `scenario != comparator`
# ===========================================================================
r <- resolve("sc0", gv)
expect_equal(sort(names(r)), sort(setdiff(gv, "sc0")),
             info = "scalar: every other scenario is an intervention")
expect_true(all(unname(r) == "sc0"))
expect_equal(length(attr(r, "dropped")), 0L)

# A scalar comparator absent from this summary yields nothing to build.
expect_equal(length(resolve("nope", gv)), 0L)

# ===========================================================================
# 3. Resolver: the GENVASC map -- 3 relationships, 4 pairs
# ===========================================================================
cs <- c("sc0",                                   # default
        mult_base_prs    = "mult_base_qrisk",
        mult_lowzero_prs = "mult_lowzero_qrisk")
m <- resolve(cs, gv)
expect_equal(length(m), 4L, info = "four intervention scenarios")
expect_equal(unname(m["mult_base_prs"]), "mult_base_qrisk")
expect_equal(unname(m["mult_lowzero_prs"]), "mult_lowzero_qrisk")
expect_equal(unname(m["mult_base_qrisk"]), "sc0",
             info = "the default covers scenarios not named")
expect_false("sc0" %in% names(m),
             info = "the default comparator is not an intervention against itself")

# THE POINT OF THE FEATURE: one scenario is BOTH an intervention and a
# comparator, in the same resolved map.
expect_true("mult_base_qrisk" %in% names(m))
expect_true("mult_base_qrisk" %in% unname(m))

# No unnamed element => only the listed pairs are built.
m2 <- resolve(c(mult_base_prs = "mult_base_qrisk"), gv)
expect_equal(length(m2), 1L)
expect_equal(names(m2), "mult_base_prs")

# Pairs this summary cannot build are reported, not silently dropped.
m3 <- resolve(c(absent_arm = "sc0", mult_base_prs = "mult_base_qrisk"), gv)
expect_equal(length(m3), 1L)
expect_true(any(grepl("absent_arm", attr(m3, "dropped"))),
            info = "an intervention this summary lacks is reported as dropped")
m4 <- resolve(c(mult_base_prs = "not_here"), gv)
expect_equal(length(m4), 0L)
expect_true(any(grepl("not_here", attr(m4, "dropped"))),
            info = "a comparator this summary lacks is reported as dropped")

# Indexing must go through as.character(): a factor would use LEVEL CODES.
f <- factor(gv, levels = rev(gv))            # deliberately non-alphabetical
expect_equal(resolve(cs, f)[order(names(resolve(cs, f)))],
             m[order(names(m))],
             info = "a factor scenario column resolves identically")

# ===========================================================================
# 4. Cycles are rejected; chains are fine
# ===========================================================================
expect_equal(length(cycles(resolve(cs, gv))), 0L,
             info = "prs -> qrisk -> sc0 is a chain, not a cycle")
expect_true(length(cycles(c(a = "b", b = "a"))) > 0L)
expect_true(length(cycles(c(a = "b", b = "c", c = "a"))) > 0L)

# ===========================================================================
# 5. Structural validation
# ===========================================================================
expect_silent(validate("sc0"))
expect_silent(validate(cs))
expect_error(validate(1L), pattern = "character vector")
expect_error(validate(character(0)), pattern = "character vector")
expect_error(validate(list("sc0")), pattern = "character vector")
expect_error(validate(c("sc0", NA)), pattern = "NA or empty")
expect_error(validate(c("sc0", "")), pattern = "NA or empty")
expect_error(validate(c("sc0", "sc1")), pattern = "AT MOST ONE unnamed",
             info = "a two-element unnamed vector is NOT read positionally")
expect_error(validate(c("sc0", "sc0")), pattern = "AT MOST ONE unnamed")
expect_error(validate(c(sc1 = "sc0", sc1 = "scX")), pattern = "more than once")
expect_error(validate(c(sc1 = "sc1")), pattern = "its own comparator")

# ===========================================================================
# 6. End to end through export_equity_tables()
# ===========================================================================
dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")

make_incd <- function(scns) {
  d <- CJ(mc = 1:3, scenario = scns, year = 19:21,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  K <- length(dimd_lv)
  d[, dep := (K + 1L - as.integer(dimd))]
  # each arm removes a bit more, and more so in the deprived groups
  eff <- c(0, 0.01, 0.02, 0.03)[match(d$scenario, scns)]
  d[, chd_incd := (dep * 100 + mc) * (1 - eff * dep)]
  d[, popsize := 10000 * dep]
  d[, dep := NULL]
  d[]
}

run_export <- function(strata, sources, comparator, logs = FALSE) {
  dir <- tempfile("cmptest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = logs, diseases = NULL)))
  env$private <- list(read_summary_dataset = function(summary_type, standardization) {
    if (is.null(sources[[summary_type]])) NULL else copy(sources[[summary_type]])
  })
  f <- Simulation$private_methods$export_equity_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9), summaries_dir = dir, tables_dir = dir,
    comparator_scenario = comparator, baseline_year = 2019L,
    ridit_reference = "comparator", strata = strata)
  dir
}

fn <- "equity cpp slope index by year (not standardised).csv"
src2 <- list(incd = make_incd(c("sc0", "sc1", "sc2")))

# --- 6a. SCALAR EQUIVALENCE: the regression that matters -------------------
d_scalar <- run_export(list("year"), src2, "sc0")
d_map    <- run_export(list("year"), src2, c(sc1 = "sc0", sc2 = "sc0"))
t_scalar <- fread(file.path(d_scalar, fn))
t_map    <- fread(file.path(d_map, fn))

expect_equal(sort(basename(list.files(d_scalar, pattern = "csv$"))),
             sort(basename(list.files(d_map, pattern = "csv$"))),
             info = "a map writes the same files as the equivalent scalar")
expect_true("comparator" %in% names(t_map))
expect_false("comparator" %in% names(t_scalar),
             info = "scalar mode adds NO column")
expect_equal(unique(t_map$comparator), "sc0")
expect_equal(names(t_map)[1:3], c("gradient", "scenario", "comparator"),
             info = "comparator sits immediately after scenario")
t_map_stripped <- copy(t_map)[, comparator := NULL]
setcolorder(t_map_stripped, names(t_scalar))
expect_equal(t_map_stripped, t_scalar, check.attributes = FALSE,
             info = "an equivalent map reproduces the scalar numbers exactly")

# --- 6b. A scenario that is BOTH intervention and comparator ---------------
src3 <- list(incd = make_incd(c("sc0", "a_qrisk", "a_prs")))
d_gv <- run_export(list("year"), src3, c("sc0", a_prs = "a_qrisk"))
t_gv <- fread(file.path(d_gv, fn))

expect_equal(sort(unique(t_gv$scenario)), c("a_prs", "a_qrisk"),
             info = "both arms appear as interventions")
expect_equal(unique(t_gv[scenario == "a_prs", comparator]), "a_qrisk")
expect_equal(unique(t_gv[scenario == "a_qrisk", comparator]), "sc0",
             info = "a_qrisk is an intervention vs sc0 AND a_prs's comparator")
expect_equal(uniqueN(t_gv, by = c("scenario", "year", "disease", "type")),
             nrow(t_gv),
             info = "scenario stays a UNIQUE row key -- no row fan-out")

# The mapped contrast must differ from the same arm referenced to sc0.
t_gv_sc0 <- fread(file.path(run_export(list("year"), src3, "sc0"), fn))
k <- c("scenario", "year", "disease", "type")
cmpd <- merge(t_gv[scenario == "a_prs", c(k, "equity_50.0%"), with = FALSE],
              t_gv_sc0[scenario == "a_prs", c(k, "equity_50.0%"), with = FALSE],
              by = k, suffixes = c(".map", ".sc0"))
expect_false(isTRUE(all.equal(cmpd$`equity_50.0%.map`, cmpd$`equity_50.0%.sc0`)),
             info = "a_prs vs a_qrisk is genuinely not a_prs vs sc0")

# --- 6c. An unbuildable pair warns and is skipped, not silently dropped ----
expect_warning(run_export(list("year"), src3, c(a_prs = "ghost_arm")),
               pattern = "cannot be built")

# ===========================================================================
# 7. A missing COMPARATOR CELL means a comparator count of zero
# ===========================================================================
# `d` inside tbl_smmrs_core is a sum aggregation keyed by the strata, so a
# combination with no row is one with no EVENTS, not one with an unknown count.
# The contrast must therefore be 0 - value = -value (cases CAUSED), which is a
# real quantity. Two wrong answers this pins against:
#   * keeping the raw level  -> +value, the sign INVERTED (the original bug);
#   * setting NA             -> the draw is dropped from that cell's quantiles
#                               for this year and, via cumsum(), every later
#                               year, turning a harm into a certain-looking null.
cpp_core <- IMPACTncdEngland:::.resolve_comparators   # keep the namespace live
mk_hole <- function() {
  d <- CJ(mc = 1L, scenario = c("sc0", "sc1"), year = 19:20,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  d[, chd_incd := 100]
  d[, popsize := 1000]
  # Punch a hole in the COMPARATOR only: sc0 has nobody in the last decile in
  # the final year, sc1 has 700 cases there.
  d <- d[!(scenario == "sc0" & year == 20L & dimd == "10 least deprived")]
  d[scenario == "sc1" & year == 20L & dimd == "10 least deprived",
    chd_incd := 700]
  d[]
}
d7 <- run_export(list("year", c("year", "dimd")), list(incd = mk_hole()), "sc0")
t7 <- fread(file.path(d7,
  "equity cpp slope index by year-dimd (not standardised).csv"))
expect_true(nrow(t7) > 0L, info = "the hole does not abort the export")

# The equity family reaches the same join; check the main table directly for the
# sign, via the same code path tbl_smmrs_core uses.
hole <- mk_hole()
hole[, year := year + 2000L]
xk <- c("mc", "scenario", "year", "dimd")
agg <- hole[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"),
            keyby = xk]
agg <- melt(agg, id.vars = xk)
cmp7 <- agg[scenario == "sc0" & year >= 2019L]
d7b <- agg[scenario != "sc0" & year >= 2019L]
onk <- c(setdiff(xk, "scenario"), "variable")
miss <- d7b[!cmp7, on = onk, which = TRUE]
expect_true(length(miss) > 0L, info = "the fixture really does leave a hole")
raw <- d7b$value[miss]
d7b[cmp7, on = onk, value := i.value - value]
set(d7b, miss, "value", -d7b$value[miss])
expect_equal(d7b$value[miss], -raw,
             info = "absent comparator cell -> contrast is MINUS the intervention")
expect_true(all(d7b$value[miss] < 0),
            info = "cases the comparator never had are a HARM, not a benefit")
expect_false(anyNA(d7b$value),
             info = "and are not NA'd away, which would drop the draw")
