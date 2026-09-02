# Which deprivation stratum the cached PARF table is joined at.
#
# BACKGROUND (the bug this pins)
# ------------------------------
# `Disease$gen_parf()` recomputes nothing from a scenario's exposures. It reads
# the cached PARF table and joins it onto sp$pop with `on = .NATURAL`, and
# `dimd` is one of those natural keys. Its payload is
#     p0 = mu      * (1 - parf)         background ("unexposed") incidence
#     m0 = mu_mrtl * (1 - parf_mrtl)    background disease mortality
# with `mu` the OBSERVED, decile-specific rate -- i.e. p0/m0 carry the part of
# the deprivation gradient the modelled risk factors do NOT explain.
#
# Commit adbf762 moved that call from before the primary-prevention scenario to
# after it (so set_init_prvl() would see the scenario world). For a scenario
# that rewrites `dimd` to model deprivation change, that silently handed every
# simulant the TARGET decile's background risk -- transferring the non-modelled
# gradient as a side effect of call order. Two arms of the inequality analysis
# that differed only in whether they moved that gradient became the same model
# and published identical results.
#
# The stratum is now stated by the scenario rather than implied by call order:
#   update_primary_prevention_scn(fn, parf_deprivation = "original")  # default
#   update_primary_prevention_scn(fn, parf_deprivation = "scenario")
# and applied by with_original_dimd() around the gen_parf/set_init_prvl block.
#
# Everything below is about the DEFAULT being a real no-op for scenarios that
# leave `dimd` alone, and a real change for scenarios that do not.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

with_orig <- IMPACTncdEngland:::with_original_dimd

DIMD_LV <- c("1 most deprived", as.character(2:9), "10 least deprived")

mk_pop <- function(scenario_dimd = NULL) {
  d <- data.table(
    pid = 1:10,
    dimd = factor(DIMD_LV, levels = DIMD_LV),
    age = 60L,
    sex = "men"
  )
  d[, dimd_orig := dimd]                    # snapshotted in run_sim()
  if (!is.null(scenario_dimd)) {
    set(d, NULL, "dimd", factor(scenario_dimd, levels = DIMD_LV))
  }
  d[]
}

# ===========================================================================
# 1. with_original_dimd(): the swap, and its exact no-op conditions
# ===========================================================================
p <- mk_pop("10 least deprived")
seen <- NULL
sw <- with_orig(p, "original", { seen <- as.character(p$dimd) })
expect_true(sw, info = "a scenario that rewrote dimd is swapped")
expect_equal(seen, DIMD_LV,
             info = "expr sees the ORIGINAL, pre-scenario deciles")
expect_equal(as.character(p$dimd), rep("10 least deprived", 10L),
             info = "the scenario's dimd is restored on the way out")
expect_false("dimd_scn" %in% names(p),
             info = "and the scratch column does not leak")
expect_equal(levels(p$dimd), DIMD_LV,
             info = "the factor keeps its full level set through the swap")

# "scenario" leaves everything alone -- that is the opt-in for moving the
# non-modelled gradient too.
p <- mk_pop("10 least deprived")
seen <- NULL
sw <- with_orig(p, "scenario", { seen <- as.character(p$dimd) })
expect_false(sw)
expect_equal(seen, rep("10 least deprived", 10L),
             info = "under 'scenario' the join sees the forced decile")

# A scenario that never touched dimd must be bit-for-bit unaffected: this is
# what makes the default safe for every other scenario in the repo.
p <- mk_pop()
seen <- NULL
sw <- with_orig(p, "original", { seen <- as.character(p$dimd) })
expect_false(sw, info = "dimd == dimd_orig: nothing to swap")
expect_equal(seen, DIMD_LV)
expect_false("dimd_scn" %in% names(p))

# No snapshot (gen_parf called outside run_sim, e.g. PARF generation) -> no-op.
p <- mk_pop("10 least deprived")[, dimd_orig := NULL]
expect_false(with_orig(p, "original", NULL),
             info = "without dimd_orig there is nothing to restore to")
expect_equal(as.character(p$dimd), rep("10 least deprived", 10L))

# The restore must survive an error, or a failing disease would leave the rest
# of the simulation running on the wrong deprivation.
p <- mk_pop("10 least deprived")
expect_error(with_orig(p, "original", stop("boom")), pattern = "boom")
expect_equal(as.character(p$dimd), rep("10 least deprived", 10L),
             info = "on.exit restores dimd even when expr throws")
expect_false("dimd_scn" %in% names(p),
             info = "and cleans up the scratch column")

# Partial rewrites (the `_future` scenarios change dimd only from sc_year on).
p <- mk_pop()
p[pid > 5L, dimd := factor("10 least deprived", levels = DIMD_LV)]
seen <- NULL
expect_true(with_orig(p, "original", { seen <- as.character(p$dimd) }))
expect_equal(seen, DIMD_LV, info = "a partial rewrite is fully reverted")
expect_equal(as.character(p$dimd),
             c(DIMD_LV[1:5], rep("10 least deprived", 5L)),
             info = "...and fully restored, row by row")

# ===========================================================================
# 2. Meaningfulness: the join really does return different background risk
# ===========================================================================
# Replays gen_parf()'s payload step verbatim:
#     colnam <- setdiff(names(tt), names(sp$pop))
#     private$parf <- tt[sp$pop, on = .NATURAL, ..colnam]
# against a PARF-shaped fixture with a deprivation gradient in p0/m0. Without
# this section the tests above would pass on a helper that swapped a column
# nothing reads.
parf_tbl <- data.table(
  age = 60L,
  sex = "men",
  dimd = factor(DIMD_LV, levels = DIMD_LV),
  # a gradient of the measured shape: most-deprived p0 ~2.6x least-deprived
  p0 = seq(0.0085, 0.0033, length.out = 10L),
  m0 = seq(0.0019, 0.0014, length.out = 10L)
)
join_parf <- function(pop) {
  keep <- intersect(names(pop), c("age", "sex", "dimd"))
  sp_like <- pop[, ..keep]
  colnam <- setdiff(names(parf_tbl), names(sp_like))
  parf_tbl[sp_like, on = .NATURAL, ..colnam]
}

p <- mk_pop("10 least deprived")

# What the post-adbf762 code did: join at the scenario's dimd.
after <- join_parf(p)
expect_equal(after$p0, rep(parf_tbl$p0[10L], 10L),
             info = "joining at the forced decile gives everyone its p0")

# What the default now does.
inside <- NULL
with_orig(p, "original", { inside <- join_parf(p) })
expect_equal(inside$p0, parf_tbl$p0,
             info = "under 'original' each simulant keeps their own decile's p0")
expect_equal(inside$m0, parf_tbl$m0,
             info = "and their own decile's background disease mortality")

# The two are genuinely different, and by a lot -- if they were not, none of
# this would matter.
expect_false(isTRUE(all.equal(inside$p0, after$p0)))
expect_true(mean(inside$p0) / mean(after$p0) > 1.4,
            info = "the leaked gradient is large, not a rounding difference")

# And 'scenario' reproduces the transfer on purpose.
p <- mk_pop("10 least deprived")
opted_in <- NULL
with_orig(p, "scenario", { opted_in <- join_parf(p) })
expect_equal(opted_in$p0, after$p0,
             info = "'scenario' is the explicit opt-in to the old behaviour")

# ===========================================================================
# 3. The scenario-facing API
# ===========================================================================
upps <- Simulation$public_methods$update_primary_prevention_scn
expect_true("parf_deprivation" %in% names(formals(upps)))
expect_equal(eval(formals(upps)$parf_deprivation), c("original", "scenario"),
             info = "'original' is the default, i.e. the historical Modeled arm")

# match.arg must reject anything else rather than fall through to a default:
# a typo here would silently choose a different model.
env <- new.env(parent = asNamespace("IMPACTncdEngland"))
# `private` must be an ENVIRONMENT, as R6 makes it: with a plain list,
# `private$x <- v` would rebind a local copy and the assertions below would be
# vacuous (they would only ever see the initial value).
env$private <- new.env()
env$private$primary_prevention_scn <- NULL
env$private$parf_deprivation <- "original"
env$private$update_primary_prevention_scn <- function() NULL
env$self <- new.env()
f <- upps
environment(f) <- env
expect_silent(f(function(sp) NULL))
expect_equal(env$private$parf_deprivation, "original")
expect_silent(f(function(sp) NULL, parf_deprivation = "scenario"))
expect_equal(env$private$parf_deprivation, "scenario")
expect_error(f(function(sp) NULL, parf_deprivation = "sceanrio"))
expect_error(f(function(sp) NULL, parf_deprivation = TRUE))

# ===========================================================================
# 4. Regression pin: run_sim() must still route gen_parf through the helper
# ===========================================================================
# The original defect was a RELOCATION of the gen_parf/set_init_prvl block, not
# a wrong expression -- nothing errored and no output looked odd. So pin the
# structure: gen_parf must be reached from inside with_original_dimd().
body_txt <- paste(
  deparse(body(Simulation$private_methods$run_sim)),
  collapse = "\n"
)
expect_true(grepl("with_original_dimd", body_txt, fixed = TRUE),
            info = "run_sim wraps the PARF join in with_original_dimd()")
expect_true(grepl("parf_deprivation", body_txt, fixed = TRUE),
            info = "...and passes the scenario's choice into it")

# The RNG restore in the same function must not regress to `set.seed(rs)`:
# set.seed() keeps only rs[1], the RNG *kind* code, so that form silently
# reseeds every scenario and iteration from one constant instead of restoring
# the saved stream.
expect_false(grepl("set.seed(rs)", body_txt, fixed = TRUE),
             info = "set.seed(<full .Random.seed>) does not restore a stream")
expect_true(grepl(".Random.seed", body_txt, fixed = TRUE),
            info = "the stream is restored by assigning .Random.seed")
