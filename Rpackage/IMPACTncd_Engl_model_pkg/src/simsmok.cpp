/* IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
 Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
 Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
 funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
 Research Simulation Environment.  The views expressed are those of the
 authors and not necessarily those of the NHS, the NIHR or the Department of
 Health.

 Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos

 IMPACTncdEngl is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later
 version. This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details. You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/> or write
 to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA 02110-1301 USA. */

#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
void simsmok(
    DataFrame& df,
  const NumericMatrix& pr_relapse,
  const int& relapse_cutoff) {

  if (pr_relapse.ncol() < relapse_cutoff) stop("relapse_cutoff should be smaller than the number of columns of pr_relapse matrix.");

  //access the df columns
  IntegerVector smok_status     = df["smok_status"];
  NumericVector prb_smok_incid  = df["prb_smok_incid"];
  NumericVector prb_smok_cess   = df["prb_smok_cess"];
  NumericVector rn_smok         = df["rankstat_smok"];
  LogicalVector new_pid         = df["pid_mrk"];
  IntegerVector sex             = df["sex"];
  IntegerVector qimd            = df["qimd"];
  IntegerVector smok_quit_yrs   = df["smok_quit_yrs"];
  IntegerVector smok_dur        = df["smok_dur"];

  // id should be sorted by year
  const int n = df.nrows();
  int nrow = 0;

  for (int i = 0; i < n; i++)
  {
    if (!new_pid[i]) // if not a new simulant
    {
      if (smok_status[i-1] == 1)
      { // never smoker the previous year
        if (rn_smok[i] < prb_smok_incid[i])
        {
          smok_status[i] = 4;
          smok_dur[i] = 1;
        }
        else
        {
          smok_status[i] = 1;
        }
      }

      if (smok_status[i-1] == 4)
      { //current smoker the previous year
        if (rn_smok[i] < prb_smok_cess[i])
        {
          smok_status[i] = 3;
          smok_quit_yrs[i] = 1;
          smok_dur[i] = smok_dur[i-1];
        }
        else
        {
          smok_status[i] = 4;
          smok_dur[i] = smok_dur[i-1] + 1;
        }
      }

      if ((smok_status[i-1] == 2 || smok_status[i-1] == 3) && (smok_quit_yrs[i-1] <= relapse_cutoff))
        {
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (rn_smok[i] < (pr_relapse(nrow, smok_quit_yrs[i-1] - 1)))
        {
          smok_status[i] = 4;
          smok_quit_yrs[i] = 0;
          smok_dur[i] = smok_dur[i-1] + 1;
        }
        else
        {
          smok_status[i] = smok_status[i-1];
          smok_quit_yrs[i] = smok_quit_yrs[i-1] + 1;
          smok_dur[i] = smok_dur[i-1];
        }
      }
      if ((smok_status[i-1] == 2 || smok_status[i-1] == 3) && (smok_quit_yrs[i-1] > relapse_cutoff))
      {
        smok_status[i] = smok_status[i-1];
        smok_quit_yrs[i] = smok_quit_yrs[i-1] + 1;
        smok_dur[i] = smok_dur[i-1];
      }
    }
  }
}

// Same as simsmok but cols are for social policies scenarios
//' @export
// [[Rcpp::export]]
void simsmok_sc(
    DataFrame& df,
    const NumericMatrix& pr_relapse,
    const int& relapse_cutoff,
    const IntegerVector& row_sel) {

  if (pr_relapse.ncol() < relapse_cutoff) stop("relapse_cutoff should be smaller than the number of columns of pr_relapse matrix.");

  //access the df columns
  IntegerVector smok_status     = df["smok_status_sc"];
  NumericVector prb_smok_incid  = df["prb_smok_incid_sc"];
  NumericVector prb_smok_cess   = df["prb_smok_cess_sc"];
  NumericVector rn_smok         = df["rankstat_smok"];
  LogicalVector new_pid         = df["pid_mrk_sc"];
  IntegerVector sex             = df["sex"];
  IntegerVector qimd            = df["qimd_sc"];
  IntegerVector smok_quit_yrs   = df["smok_quit_yrs_sc"];
  IntegerVector smok_dur        = df["smok_dur_sc"];

  // id should be sorted by year
  const int n = row_sel.size();
  int nrow = 0;
  int i = 0;
  for (int j = 0; j < n; j++)
  {
    i = row_sel[j] - 1;
    if (!new_pid[i]) // if not a new simulant
    {
      if (smok_status[i-1] == 1)
      { // never smoker the previous year
        if (rn_smok[i] < prb_smok_incid[i])
        {
          smok_status[i] = 4;
          smok_dur[i] = 1;
        }
        else
        {
          smok_status[i] = 1;
        }
      }

      if (smok_status[i-1] == 4)
      { //current smoker the previous year
        if (rn_smok[i] < prb_smok_cess[i])
        {
          smok_status[i] = 3;
          smok_quit_yrs[i] = 1;
          smok_dur[i] = smok_dur[i-1];
        }
        else
        {
          smok_status[i] = 4;
          smok_dur[i] = smok_dur[i-1] + 1;
        }
      }

      if ((smok_status[i-1] == 2 || smok_status[i-1] == 3) && (smok_quit_yrs[i-1] <= relapse_cutoff))
      {
        switch (sex[i]) {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (rn_smok[i] < (pr_relapse(nrow, smok_quit_yrs[i-1] - 1)))
        {
          smok_status[i] = 4;
          smok_quit_yrs[i] = 0;
          smok_dur[i] = smok_dur[i-1] + 1;
        }
        else
        {
          smok_status[i] = smok_status[i-1];
          smok_quit_yrs[i] = smok_quit_yrs[i-1] + 1;
          smok_dur[i] = smok_dur[i-1];
        }
      }
      if ((smok_status[i-1] == 2 || smok_status[i-1] == 3) && (smok_quit_yrs[i-1] > relapse_cutoff))
      {
        smok_status[i] = smok_status[i-1];
        smok_quit_yrs[i] = smok_quit_yrs[i-1] + 1;
        smok_dur[i] = smok_dur[i-1];
      }
    }
  }
}


//' @export
// [[Rcpp::export]]
void simsmok_postcalibration(
    DataFrame& df) {

  //access the df columns
  IntegerVector smok_status     = df["smok_status"];
  LogicalVector new_pid         = df["pid_mrk"];
  IntegerVector smok_quit_yrs   = df["smok_quit_yrs"];
  IntegerVector smok_dur        = df["smok_dur"];
  IntegerVector smok_cig        = df["smok_cig"];


  // id should be sorted by year
  const int n = df.nrows();
  int nrow = 0;

  for (int i = 0; i < n; i++)
  {
    if (!new_pid[i]) // if not a new simulant
    {
      if (smok_status[i-1] == smok_status[i] && (smok_status[i] == 2 || smok_status[i] == 3))
      {
        smok_quit_yrs[i] = smok_quit_yrs[i - 1] + 1;
        smok_dur[i] = smok_dur[i - 1];
        smok_cig[i] = smok_cig[i - 1];
      }

      if (smok_status[i-1] == smok_status[i] && smok_status[i] == 4)
      {
        smok_quit_yrs[i] = 0; // happens anyway in R side
        smok_dur[i] = smok_dur[i - 1] + 1;
      }


    }
  }
}



//' @export
// [[Rcpp::export]]
void simsmok_cig(DataFrame& df) {

  //access the df columns
  IntegerVector smok_status     = df["smok_status"];
  IntegerVector smok_cig        = df["smok_cig"];
  LogicalVector new_pid         = df["pid_mrk"];

  // id should be sorted by year
  const int n = df.nrows();

  for (int i = 0; i < n; i++)
  {
    if (!new_pid[i]) // if not a new simulant
    { // if smok_status[i] == 3 then previous year was either 3 or 4. In both cases smog_cig should carry forward
      if (smok_status[i] == 3) smok_cig[i] = smok_cig[i-1];
    }
  }
}

//' @export
// [[Rcpp::export]]
void simsmok_cig_sc(DataFrame& df, const IntegerVector& row_sel) {

  //access the df columns
  IntegerVector smok_status     = df["smok_status_sc"];
  IntegerVector smok_cig        = df["smok_cig_sc"];
  LogicalVector new_pid         = df["pid_mrk_sc"];

  // id should be sorted by year
  const int n = row_sel.size();
  int i = 0;
  for (int j = 0; j < n; j++)
  {
    i = row_sel[j] - 1;
    if (!new_pid[i]) // if not a new simulant
    { // if smok_status[i] == 3 then previous year was either 3 or 4. In both cases smog_cig should carry forward
      if (smok_status[i] == 3) smok_cig[i] = smok_cig[i-1];
    }
  }
}

//' @export
// [[Rcpp::export]]
List simsmok_cessation(const IntegerVector& smok_status,
                       const IntegerVector& smok_quit_yrs,
                       const IntegerVector& smok_dur,
                       const IntegerVector& sex,
                       const IntegerVector& qimd,
                       const LogicalVector& new_pid,
                       const IntegerVector& hc_eff,
                       const NumericVector& relapse_rn,
                       const NumericMatrix& pr_relapse,
                       const int&           relapse_cutoff)
  {
  if (pr_relapse.ncol() < relapse_cutoff) stop("relapse_cutoff should be smaller than the number of columns of pr_relapse matrix.");

  // id should be sorted by year
  const int n = smok_status.size();
  IntegerVector out_status = clone(smok_status);
  IntegerVector out_quit_yrs = clone(smok_quit_yrs);
  IntegerVector out_dur = clone(smok_dur);
  bool marker = false;
  int nrow = 0;

  for (int i = 0; i < n; i++)
  {
    if (!new_pid[i]) // if not a new simulant
    {
      if (out_status[i] == 4 && hc_eff[i] == 1) // out_status[i] == 4 not ensured in R side because hc_eff carried forward
      {
        marker = true; // marker == true only for pids that quit after a  HC
        out_status[i] = 3;
        out_quit_yrs[i] = 1;
        out_dur[i] = out_dur[i-1];
      }
      if (marker && out_status[i-1] == 3 && out_quit_yrs[i-1] <= relapse_cutoff)
      {
        switch (sex[i])
        {
        case 1: nrow = qimd[i] - 1; break;
        case 2: nrow = qimd[i] + 4; break;
        }
        if (relapse_rn[i] < (pr_relapse(nrow, smok_quit_yrs[i-1] - 1)))
        {
          out_status[i] = 4;
          out_quit_yrs[i] = 0;
          out_dur[i] = out_dur[i-1] + 1;
        }
        else
        {
          out_status[i] = 3;
          out_quit_yrs[i] = out_quit_yrs[i-1] + 1;
          out_dur[i] = out_dur[i-1];
        }
        if (out_status[i-1] == 3 && marker && smok_quit_yrs[i-1] > relapse_cutoff)
          {
          out_status[i] = 3;
          out_quit_yrs[i] = out_quit_yrs[i-1] + 1;
          out_dur[i] = out_dur[i-1];
          }
      }
    }
    else // if new pid
    {
      marker = false; // reset market for each pid
      if (out_status[i] == 4 && hc_eff[i] == 1) // out_status[i] == 4 not ensured in R side because hc_eff carried forward
      {
        marker = true; // marker == true only for pids that quit after a  HC
        out_status[i] = 3;
        out_quit_yrs[i] = 1;
        out_dur[i] = smok_dur[i] - 1;
      }
    }
  }
  return List::create(_["smok_status"]= out_status,
                      _["smok_quit_yrs"]= out_quit_yrs,
                      _["smok_dur"]= out_dur);
}

//' @export
 // [[Rcpp::export]]
 List simsmok_complete_cessation(const IntegerVector& smok_status,
                        const IntegerVector& smok_quit_yrs,
                        const IntegerVector& smok_dur,
                        const LogicalVector& new_pid,
                        const IntegerVector& year,
                        const int policy_first_year)
 {
   // id should be sorted by year
   const int n = smok_status.size();
   IntegerVector out_status = clone(smok_status);
   IntegerVector out_quit_yrs = clone(smok_quit_yrs);
   IntegerVector out_dur = clone(smok_dur);
   bool marker = false;
   int nrow = 0;

   for (int i = 0; i < n; i++)
   {
     if (year[i] < policy_first_year) continue; // skip years before policy

     if (!new_pid[i]) // if not a new simulant
     {
       if (out_status[i] == 4)
       {
         out_status[i] = 3;
         out_quit_yrs[i] = 1;
         out_dur[i] = out_dur[i-1];
       }
       if (out_status[i] = 3)
       {
         {
           out_quit_yrs[i] = out_quit_yrs[i-1] + 1;
           out_dur[i] = out_dur[i-1];
         }
       }
     }
     else // if new pid
     {
       if (out_status[i] == 4)
       {
         out_status[i] = 3;
         out_quit_yrs[i] = 1;
         out_dur[i] = smok_dur[i] - 1;
       }
     }
   }
   return List::create(_["smok_status"]= out_status,
                       _["smok_quit_yrs"]= out_quit_yrs,
                       _["smok_dur"]= out_dur);
 }


// ex smokers are turning into smokers. Used in structural smoking when smoking
// is increasing.
//' @export
// [[Rcpp::export]]
List simsmok_policy_impact_incr(const IntegerVector& smok_status,
                                const IntegerVector& smok_quit_yrs,
                                const IntegerVector& smok_dur,
                                const LogicalVector& new_pid,
                                const IntegerVector& hc_eff)
{
  // id should be sorted by year
  const int n = smok_status.size();
  IntegerVector out_status = clone(smok_status);
  IntegerVector out_quit_yrs = clone(smok_quit_yrs);
  IntegerVector out_dur = clone(smok_dur);
  bool marker = false;

  for (int i = 0; i < n; i++)
  {
    if (hc_eff[i] == 1) // out_status[i] == 3 ensured in R side
    {
      marker = true; // marker == true only for pids that quit after a  HC
      out_status[i] = 4;
      out_quit_yrs[i] = 0;
      out_dur[i] = smok_dur[i] + 1;
    }

    if (marker && hc_eff[i] == 0 && !new_pid[i])
    {
      if (out_status[i] == 4)
      {
        out_quit_yrs[i] = 0;
        out_dur[i] = out_dur[i-1] + 1;
      }
      if (out_status[i] == 3)
      {
        out_quit_yrs[i] = out_quit_yrs[i - 1] + 1;
        out_dur[i] = out_dur[i-1];
      }
    }

    if (new_pid[i]) marker = false; // reset market for each pid
  }
  return List::create(_["smok_status"]= out_status,
                      _["smok_quit_yrs"]= out_quit_yrs,
                      _["smok_dur"]= out_dur);
}


// smokers are turning into ex smokers. Used in structural smoking when smoking
// is increasing.
//' @export
// [[Rcpp::export]]
List simsmok_policy_impact_decr(const IntegerVector& smok_status,
                                const IntegerVector& smok_quit_yrs,
                                const IntegerVector& smok_dur,
                                const IntegerVector& smok_cig,
                                const LogicalVector& new_pid,
                                const IntegerVector& hc_eff)
{
  // id should be sorted by year
  const int n = smok_status.size();
  IntegerVector out_status = clone(smok_status);
  IntegerVector out_quit_yrs = clone(smok_quit_yrs);
  IntegerVector out_dur = clone(smok_dur);
  IntegerVector out_cig = clone(smok_cig);
  bool marker = false;

  for (int i = 0; i < n; i++)
  {
    if (hc_eff[i] == 1 && new_pid[i]) // out_status[i] == 4 ensured in R side
    {
      marker = true; // marker == true only for pids that quit after a  HC
      out_status[i] = 3; // NOTE to avoid smokers that go from status 1 to 3 without 4
      out_quit_yrs[i] = 1;
      out_dur[i] = smok_dur[i] - 1;
    }

    if (hc_eff[i] == 1 && !new_pid[i]) // out_status[i] == 4 ensured in R side
    {
      marker = true; // marker == true only for pids that quit after a  HC
      if (out_status[i - 1] != 1)
      {
        out_status[i] = 3; // NOTE to avoid smokers that go from status 1 to 3 without 4
        out_quit_yrs[i] = 1;
        out_dur[i] = smok_dur[i] - 1;
      }
      else // if previous status was 1
      {
        out_status[i] = 1; // NOTE to avoid smokers that go from status 1 to 3 without 4
        out_quit_yrs[i] = 0;
        out_dur[i] = 0;
        out_cig[i] = 0;
      }

    }

    if (marker && hc_eff[i] == 0 && !new_pid[i])
    {
      if (out_status[i-1] == 1)
      {
        out_status[i] = 1; // NOTE to avoid smokers that go from status 1 to 3 without 4
        out_quit_yrs[i] = 0;
        out_dur[i] = 0;
        out_cig[i] = 0;
      }
      if (out_status[i] == 4) // NOTE out_status[i] not out_status[i-1]
      {
        out_quit_yrs[i] = 0;
        out_dur[i] = out_dur[i-1] + 1;
      }
      if (out_status[i] == 3) // NOTE out_status[i] not out_status[i-1]
      {
        out_quit_yrs[i] = out_quit_yrs[i - 1] + 1;
        out_dur[i] = out_dur[i-1];
      }
    }

    if (new_pid[i]) marker = false; // reset market for each pid
  }
  return List::create(_["smok_status"]= out_status,
                      _["smok_quit_yrs"]= out_quit_yrs,
                      _["smok_dur"]= out_dur,
                      _["smok_cig"]= out_cig);
}
