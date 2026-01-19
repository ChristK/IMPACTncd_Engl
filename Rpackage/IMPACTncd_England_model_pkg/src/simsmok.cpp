/* IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by Chris
 Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
 Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
 funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngland: Health Outcomes
 Research Simulation Environment.  The views expressed are those of the
 authors and not necessarily those of the NHS, the NIHR or the Department of
 Health.

 Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos

 IMPACTncdEngland is free software; you can redistribute it and/or modify it under
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

//' Simulate smoking status transitions over time (baseline scenario)
//' 
//' Simulates smoking behavior transitions through a lifecourse, including:
//' - Smoking initiation (never smoker -> current smoker)
//' - Smoking cessation (current smoker -> ex-smoker)
//' - Relapse (ex-smoker -> current smoker)
//' 
//' SMOKING STATUS CODES:
//' 1 = Never smoker (has never smoked regularly)
//' 2 = Occasional smoker (rare/social smoker)
//' 3 = Ex-smoker (former regular smoker who has quit)
//' 4 = Current smoker (currently smoking regularly)
//' 
//' LOGIC:
//' For each person-year observation (must be sorted by pid and year):
//' 1. If status at time t-1 equals 1 (never smoker):
//'    - May initiate smoking based on prb_smok_incid probability
//'    - If initiated: status at time t = 4, smok_dur at time t = 1
//'    - Otherwise: status at time t = 1 (remains never smoker)
//' 
//' 2. If status at time t-1 equals 4 (current smoker):
//'    - May quit based on prb_smok_cess probability
//'    - If quit: status at time t = 3, smok_quit_yrs at time t = 1, smok_dur carries forward
//'    - Otherwise: status at time t = 4, smok_dur at time t = smok_dur at time t-1 plus 1
//' 
//' 3. If status at time t-1 equals 2 or 3 (occasional/ex-smoker):
//'    a) If quit_yrs at time t-1 <= relapse_cutoff:
//'       - May relapse based on pr_relapse matrix (stratified by sex, qimd, quit years)
//'       - If relapse: status at time t = 4, smok_quit_yrs at time t = 0, smok_dur increments
//'       - Otherwise: status remains same, smok_quit_yrs increments, smok_dur carries forward
//'    b) If quit_yrs at time t-1 > relapse_cutoff:
//'       - No relapse risk, status carries forward, smok_quit_yrs increments
//' 
//' @param df DataFrame with person-year observations (sorted by pid, year)
//' @param pr_relapse Matrix of relapse probabilities (10 rows, relapse_cutoff columns)
//'                   Rows 1-5: males by qimd quintile (1=most deprived, 5=least)
//'                   Rows 6-10: females by qimd quintile
//'                   Columns: years since quit (1, 2, ..., relapse_cutoff)
//' @param relapse_cutoff Maximum years since quit with relapse risk (typically 5-10 years)
//' 
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
  NumericVector rn_smok         = df["rankstat_simsmok"];
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

//' Simulate smoking status transitions (scenario variant)
//' 
//' Identical logic to simsmok() but operates on scenario-specific columns (_sc suffix)
//' and processes only selected rows via row_sel parameter. Used for policy scenarios
//' where only a subset of the population is affected.
//' 
//' See simsmok() documentation for detailed logic explanation.
//' 
//' @param df DataFrame with scenario columns (smok_status_sc, prb_smok_incid_sc, etc.)
//' @param pr_relapse Relapse probability matrix (same as simsmok)
//' @param relapse_cutoff Maximum years with relapse risk
//' @param row_sel 1-based row indices to process (typically policy-affected individuals)
//' 
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
  NumericVector rn_smok         = df["rankstat_simsmok"];
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


//' Post-calibration smoking variable updates
//' 
//' Called after initial smoking status assignment (from calibration/input data) to ensure
//' consistency of derived variables (smok_quit_yrs, smok_dur, smok_cig) across person-years.
//' 
//' LOGIC:
//' For each continuing person (not new_pid):
//' 1. If status unchanged and = 2 or 3 (occasional/ex-smoker):
//'    - Increment smok_quit_yrs by 1
//'    - Carry forward smok_dur (doesn't increase while quit)
//'    - Carry forward smok_cig (last number smoked before quitting)
//' 
//' 2. If status unchanged and = 4 (current smoker):
//'    - Reset smok_quit_yrs to 0
//'    - Increment smok_dur by 1 (adds year of smoking)
//' 
//' This ensures temporal consistency when smoking status is externally assigned rather
//' than simulated via incidence/cessation probabilities.
//' 
//' @param df DataFrame with smoking variables (assumes status already set)
//' 
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



//' Carry forward cigarette consumption for ex-smokers (baseline)
//' 
//' When a person becomes an ex-smoker (status = 3), their cigarette consumption
//' (smok_cig) should remain at the level from when they quit. This function ensures
//' that smok_cig carries forward from the previous year for ex-smokers.
//' 
//' LOGIC:
//' - If statusat time t = 3 (ex-smoker): smok_cigat time t = smok_cigat time t-1
//' - Previous year could have been status 3 (continuing ex-smoker) or 4 (just quit)
//' - In both cases, preserve the cigarette count from when actively smoking
//' 
//' @param df DataFrame with smok_status and smok_cig columns
//' 
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
    { // if smok_status[i] == 3 then previous year was either 3 or 4. In both cases smok_cig should carry forward
      if (smok_status[i] == 3) smok_cig[i] = smok_cig[i-1];
    }
  }
}

//' Carry forward cigarette consumption for ex-smokers (scenario variant)
//' 
//' Scenario-specific version of simsmok_cig() operating on _sc columns and selected rows.
//' See simsmok_cig() for logic details.
//' 
//' @param df DataFrame with scenario columns
//' @param row_sel 1-based row indices to process
//' 
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
    { // if smok_status[i] == 3 then previous year was either 3 or 4. In both cases smok_cig should carry forward
      if (smok_status[i] == 3) smok_cig[i] = smok_cig[i-1];
    }
  }
}

//' Carry forward cigarette consumption for complete cessation scenario
//' 
//' Variant used with complete cessation policies where all smokers quit.
//' Operates on _curr_xps exposure columns.
//' 
//' @param df DataFrame with current exposure columns
//' 
//' @export
 // [[Rcpp::export]]
 void simsmok_complete_cessation_cig(DataFrame& df) {

   //access the df columns
   IntegerVector smok_status     = df["smok_status_curr_xps"];
   IntegerVector smok_cig        = df["smok_cig_curr_xps"];
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

//' Simulate healthcare-induced smoking cessation with relapse
//' 
//' Models smoking cessation interventions (e.g., stop-smoking services, medication)
//' where certain individuals receive treatment (hc_eff = 1) that causes them to quit.
//' Tracks these individuals and simulates relapse probability over subsequent years.
//' 
//' LOGIC:
//' A boolean 'marker' tracks individuals who quit due to healthcare intervention:
//' 
//' For continuing persons (not new_pid):
//' 1. If status = 4 (current smoker) and hc_eff = 1:
//'    - Set marker = TRUE (flag this person as intervention-quitter)
//'    - Change status to 3 (ex-smoker)
//'    - Set quit_yrs = 1 (first year quit)
//'    - Carry forward smoking duration
//' 
//' 2. If marker = TRUE and previous status = 3 (ex-smoker):
//'    a) If quit_yrs at time t-1 <= relapse_cutoff:
//'       - Calculate relapse probability from pr_relapse matrix (by sex, qimd, quit years)
//'       - If relapse: status = 4, quit_yrs = 0, increment smoking duration
//'       - Otherwise: status = 3, increment quit_yrs, carry forward duration
//'    
//'    b) If quit_yrs at time t-1 > relapse_cutoff:
//'       - Beyond relapse risk period
//'       - Status = 3, increment quit_yrs, carry forward duration
//' 
//' For new persons (new_pid = TRUE):
//' - Reset marker = FALSE
//' - If status = 4 and hc_eff = 1: apply cessation as above
//' 
//' @param smok_status Smoking status vector (1=never, 2=occasional, 3=ex, 4=current)
//' @param smok_quit_yrs Years since quitting
//' @param smok_dur Years of smoking (cumulative)
//' @param sex Sex (1=male, 2=female)
//' @param qimd Index of Multiple Deprivation quintile (1=most deprived, 5=least)
//' @param new_pid Logical vector marking first year of each person
//' @param hc_eff Healthcare effectiveness indicator (1=receives intervention, 0=not)
//' @param relapse_rn Random number for relapse (compared to pr_relapse)
//' @param pr_relapse Relapse probability matrix (10 rows, relapse_cutoff columns)
//' @param relapse_cutoff Maximum years with relapse risk
//' 
//' @return List with modified smok_status, smok_quit_yrs, smok_dur
//' 
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
        if (relapse_rn[i] < (pr_relapse(nrow, out_quit_yrs[i-1] - 1)))
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
      }
      if (marker && out_status[i-1] == 3 && out_quit_yrs[i-1] > relapse_cutoff)
      {
        out_status[i] = 3;
        out_quit_yrs[i] = out_quit_yrs[i-1] + 1;
        out_dur[i] = out_dur[i-1];
      }
    }
    else // if new pid
    {
      marker = false; // reset marker for each pid
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

//' Simulate complete smoking cessation policy (tobacco endgame scenario)
//' 
//' Models a policy where ALL current smokers quit from a specified year onwards
//' (e.g., tobacco ban, smoking age increase to 100). No relapse is simulated - 
//' all smokers become permanent ex-smokers.
//' 
//' LOGIC:
//' Only processes years greater than or equal to policy_first_year.
//' 
//' For continuing persons (not new_pid):
//' 1. If status equals 4 (current smoker):
//'    - Change to status 3 (ex-smoker)
//'    - Set quit_yrs to 1
//'    - Carry forward smoking duration from previous year
//' 
//' 2. If status equals 3 (ex-smoker):
//'    - Remain status 3
//'    - Increment quit_yrs
//'    - Carry forward smoking duration
//' 
//' For new persons (new_pid = TRUE):
//' - If status = 4: convert to status 3, set quit_yrs = 1, duration = input - 1
//' 
//' @param smok_status Smoking status vector
//' @param smok_quit_yrs Years since quitting
//' @param smok_dur Years of smoking
//' @param new_pid Logical vector marking first year of each person
//' @param year Simulation year for each observation
//' @param policy_first_year First year policy takes effect (observations before this are skipped)
//' 
//' @return List with modified smok_status, smok_quit_yrs, smok_dur
//' 
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
       if (out_status[i] == 3)
       {
         out_quit_yrs[i] = out_quit_yrs[i-1] + 1;
         out_dur[i] = out_dur[i-1];
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


//' Simulate policy impact when smoking prevalence INCREASES (structural changes)
//' 
//' Used for scenarios where smoking is increasing (e.g., reversal of tobacco control,
//' increased availability, marketing). Ex-smokers (status = 3) are converted back to
//' current smokers (status = 4) based on hc_eff indicator.
//' 
//' LOGIC:
//' A boolean 'marker' tracks individuals affected by the policy:
//' 
//' 1. If hc_eff = 1 (policy affects this person):
//'    - Assumes input status = 3 (ex-smoker, ensured in R side)
//'    - Set marker = TRUE
//'    - Change status to 4 (current smoker - relapse due to policy)
//'    - Set quit_yrs = 0 (no longer quit)
//'    - Increment smoking duration by 1
//' 
//' 2. If marker = TRUE and hc_eff = 0 (subsequent years after policy effect):
//'    - If status = 4: maintain quit_yrs = 0, increment duration
//'    - If status = 3: increment quit_yrs, carry forward duration
//' 
//' 3. Reset marker = FALSE when encountering new person (new_pid = TRUE)
//' 
//' @param smok_status Smoking status vector (input should have status = 3 where hc_eff = 1)
//' @param smok_quit_yrs Years since quitting
//' @param smok_dur Years of smoking
//' @param new_pid Logical vector marking first year of each person
//' @param hc_eff Policy effect indicator (1=convert to smoker, 0=no effect)
//' 
//' @return List with modified smok_status, smok_quit_yrs, smok_dur
//' 
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

    if (new_pid[i]) marker = false; // reset marker for each pid
  }
  return List::create(_["smok_status"]= out_status,
                      _["smok_quit_yrs"]= out_quit_yrs,
                      _["smok_dur"]= out_dur);
}


//' Simulate policy impact when smoking prevalence DECREASES (structural changes)
//' 
//' Used for scenarios where smoking is decreasing (e.g., tobacco control policies,
//' price increases, smoking bans). Current smokers (status = 4) are converted to
//' ex-smokers (status = 3) based on hc_eff indicator.
//' 
//' LOGIC:
//' A boolean 'marker' tracks individuals affected by the policy:
//' 
//' For new persons (new_pid = TRUE):
//' 1. If hc_eff = 1:
//'    - Assumes input status = 4 (current smoker, ensured in R side)
//'    - Set marker = TRUE
//'    - Change status to 3 (ex-smoker)
//'    - Set quit_yrs = 1
//'    - Set duration = input duration - 1
//'    - NOTE: Prevents illegal transition 1 -> 3 (never -> ex without being current)
//' 
//' For continuing persons (not new_pid):
//' 1. If hc_eff = 1:
//'    - Set marker = TRUE
//'    - If previous status != 1 (not never smoker):
//'      * Change status to 3, quit_yrs = 1, duration = input - 1
//'    - If previous status = 1 (never smoker):
//'      * Keep status = 1, reset all smoking variables to 0
//'      * Prevents illegal 1 -> 3 transition
//' 
//' 2. If marker = TRUE and hc_eff = 0 (years following policy effect):
//'    - If previous status = 1: maintain status 1, all variables = 0
//'    - If status = 4: maintain quit_yrs = 0, increment duration
//'    - If status = 3: increment quit_yrs, carry forward duration
//' 
//' 3. Reset marker = FALSE when encountering new person
//' 
//' @param smok_status Smoking status vector (input should have status = 4 where hc_eff = 1)
//' @param smok_quit_yrs Years since quitting
//' @param smok_dur Years of smoking
//' @param smok_cig Cigarettes per day
//' @param new_pid Logical vector marking first year of each person
//' @param hc_eff Policy effect indicator (1=convert to ex-smoker, 0=no effect)
//' 
//' @return List with modified smok_status, smok_quit_yrs, smok_dur, smok_cig
//' 
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

    if (new_pid[i]) marker = false; // reset marker for each pid
  }
  return List::create(_["smok_status"]= out_status,
                      _["smok_quit_yrs"]= out_quit_yrs,
                      _["smok_dur"]= out_dur,
                      _["smok_cig"]= out_cig);
}
