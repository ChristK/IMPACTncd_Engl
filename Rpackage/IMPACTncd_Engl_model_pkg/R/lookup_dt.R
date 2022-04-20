## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.



fct_to_int <- function(x, byref = FALSE) {
  # converts factor to integer
  if (!byref) x <- copy(x)
  setattr(x, name = "levels", value = NULL)
  setattr(x, name = "class", value = NULL)
  x
}

# ensure integers start from 1 in lookup_tbl
starts_from_1 <- function(tbl, on, i, min_lookup, cardinality) {
  minx <- min_lookup[[i]]
  if (is.integer(tbl[[on[[i]]]])) {
    if (minx == 1L) {
      out <- tbl[[on[[i]]]]
    } else {
      minx <- minx - 1L
      out <- tbl[[on[[i]]]] - minx
    }
    # Account for rows in tbl with no match on dt_table because integer key cols are out of range
    # TODO consider similar treat for factors. All we need is the levels in
    # lookp_tbl to be included in the levels in tbl, without gaps.
    # TODO use C++ to replace by reference the line below
    out[out < 1L | out > cardinality[[i]]] <- NA_integer_
    return(out)
  } else {
    if (minx == 1L) {
      return(fct_to_int(tbl[[on[[i]]]]))
    } else { # is this necessary? minx = 1 always for factors
      minx <- minx - 1L
      return(fct_to_int(tbl[[on[[i]]]] - minx))
    }
  }
}

# lookup_tbl = CJ(b=1:4, a = factor(letters[1:4]))[, c:=rep(1:4, 4)]
# tbl = data.table(b=0:5, a = factor(letters[1:4]))

#' @export
lookup_dt <- function(tbl,
  lookup_tbl,
  merge = TRUE,
  exclude_col = NULL,
  check_lookup_tbl_validity = FALSE) {
  # algo assumes all keys are factors or integers that start from 1 and increase
  # by 1
  # lookup_tbl needs to have unique combination of keys
  # TODO check that works for special case when only one key
  # find common cols
  if (!is.data.table(tbl))
    stop("tbl need to be a data.table")
  if (!is.data.table(lookup_tbl))
    stop("lookup_tbl need to be a data.table")

  nam_x <- names(tbl)
  nam_i <- names(lookup_tbl)
  on <- sort(setdiff(intersect(nam_x, nam_i), exclude_col))
  # Ensure year is always the first key if present
  on <- on[order(match(on, "year"))]

  return_cols_nam <- setdiff(nam_i, on)
  return_cols <- which(nam_i %in% return_cols_nam)

  if (length(on) == 0L) stop("No common keys in the two tables")
  if (length(on) == length(nam_i))
    stop("No value cols identified in lookup_tbl. Most likely all column names in lookup_tbl are present in tbl. Consider using arg exclude_col")

  if (check_lookup_tbl_validity) is_valid_lookup_tbl(lookup_tbl, on)

  # prepare lookup_tbl
  setkeyv(lookup_tbl, cols = on) # alphabetic order (breaks early if already keyed on on)
  cardinality <- vector("integer", length(on))
  names(cardinality) <- on
  min_lookup <- cardinality

  for (j in on) {
    if (is.factor(lookup_tbl[[j]])) {
      lv <- levels(lookup_tbl[[j]])
      if (check_lookup_tbl_validity &&
          !identical(lv, levels(tbl[[j]])))
        stop(j, " has different levels in tbl and lookup_tbl!")
      cardinality[[j]] <- length(lv)
      min_lookup[[j]] <- 1L
    } else {
      # only works for sorted int that increase by 1 (uniqueN is much slower)
      xmax <- last(lookup_tbl[[j]])
      xmin <- first(lookup_tbl[[j]])
      if (check_lookup_tbl_validity &&
         (min(tbl[[j]]) < xmin || max(tbl[[j]]) > xmin))
        warning(j, " has rows in tbl without a match in lookup_tbl!")
      cardinality[[j]] <- xmax - xmin + 1L
      min_lookup[[j]] <- xmin
    }
  }
  cardinality_prod <-
    shift(rev(cumprod(rev(cardinality))), -1, fill = 1L)



  # core algo
  rownum <-
    as.integer(starts_from_1(tbl, on, 1L, min_lookup, cardinality) * cardinality_prod[[1L]])
  if (length(on) > 1L) {
    for (i in 2:length(on)) {
      rownum <- as.integer(rownum -
          (cardinality[[i]] - starts_from_1(tbl, on, i, min_lookup, cardinality)) * cardinality_prod[[i]])
    }
  }


  # absorb value cols into tbl
  if (merge) {
    # lookup_tbl[, (on) := NULL]
    # # tbl[, (return_cols) := lookup_tbl[rownum]]
    # tbl[, (return_cols) := ss(lookup_tbl, rownum)]
    # for (j in return_cols) {
    #   set(tbl, NULL, j, subset_vec(lookup_tbl[[j]], as.integer(rownum)))
    #   } # Slower!
    tbl[, (return_cols_nam) := dtsubset(lookup_tbl, rownum, return_cols)]
    return(invisible(tbl))
  } else {
    return(invisible(dtsubset(lookup_tbl, rownum, return_cols)))
    # tbl[, (return_cols) := ss(lookup_tbl, rownum, return_cols)]
    # settransform(tbl, ss(lookup_tbl, rownum, return_cols)) # Not by reference
  }
}


#' @export
is_valid_lookup_tbl <- function(lookup_tbl, keycols) {
  if (!is.data.table(lookup_tbl)) stop("lookup_tbl should be a data.table.")
  if (missing(keycols) || length(keycols) == 0L) stop("keycols argument is missing.")
  keycols <- sort(keycols)
  # Ensure year is always the first key if present
  keycols <- keycols[order(match(keycols, "year"))]

  if (any(duplicated(lookup_tbl, by = keycols)))
    stop("Lookup table need to have unique combination of key columns!") # test unique keycols

  expected_rows <- 1L
  # l <- list()
  for (j in keycols) {
    # l[[j]] <- unique(lookup_tbl[[j]])
    # check type of keys
    if (typeof(lookup_tbl[[j]]) != "integer")
      stop(
        paste0(
          "Lookup table needs to have key cols of type integer (class integer or factor). ",
          j,
          " is not an integer!"
        )
      )
    expected_rows <- fifelse(is.integer(lookup_tbl[[j]]), uniqueN(lookup_tbl[[j]]),
                             length(levels(lookup_tbl[[j]]))) * expected_rows

    # check integer keys have all integer values between min and max
    if (is.integer(lookup_tbl[[j]])) {
      x <- sort.int(lookup_tbl[[j]])
      x <- x - shift(x, 1, fill = first(x))
      if (max(x) > 1L)
        stop(
          paste0(
            "Lookup table key cols of class integer need to have all integer values between min and max. ",
            j,
            " does not have all values between min(", j, ") and max(", j, ")!"
          )
        )
    }

    if (!identical(key(lookup_tbl), keycols)) warning("for best performance, set key to lookup_tbl to ", keycols)
  }

  if (nrow(lookup_tbl) != expected_rows) stop(paste0("If all possible combinations of key columns would be present in the lookup table it should have ", expected_rows, " rows. This table has ", nrow(lookup_tbl), " rows."))

  return(TRUE)
}

