// [[Rcpp::depends(dqrng, BH, sitmo)]]
// [[Rcpp::plugins(cpp11)]]

#include <mystdint.h>
#include <dqrng_generator.h>
#include <dqrng_distribution.h>
#include <pcg_random.hpp>
#include <convert_seed.h>
#include <R_randgen.h>
#include <minimal_int_set.h>
#include <dqrng.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// combine 2 uint to a uint64
long long u32tou64(const unsigned int low, const unsigned int high)
{
  return (((uint64_t) high) << 32) | ((uint64_t) low);
}


// Prepare the rng
namespace
{
  dqrng::rng64_t init()
  {
    Rcpp::RNGScope rngScope;
    Rcpp::IntegerVector seed(2, dqrng::R_random_int);
    return dqrng::generator(dqrng::convert_seed<uint64_t>(seed));
  }
  dqrng::rng64_t rng = init();

  using generator = double(*)();
  dqrng::uniform_distribution uniform{};
  generator runif_impl = [] () {return uniform(*rng);};
}

struct infl
{
  vector<IntegerVector> disease_prvl;
  vector<NumericVector> mltp;
  vector<int> lag;
};

struct disease_epi
{
  string type;
  IntegerVector prvl;
  NumericVector prbl1; // for 1st year case fatality
  NumericVector prbl2; // for prevalent cases case fatality
  infl influenced_by;
  CharacterVector aggregate;
  double mm_wt;
  bool can_recur;
  bool flag; // set true the first time incidence occurs & set true to cure in the next year
  int cure;
  int death_code;
};

struct disease_meta
{
  disease_epi incd;
  disease_epi dgns;
  disease_epi mrtl;
  bool mrtl1flag;
  int seed;
};

struct simul_meta
{
  int init_year;
  int age_low;
  IntegerVector pid;
  IntegerVector year;
  IntegerVector age;
  IntegerVector dead;
  IntegerVector mm_count;
  NumericVector mm_score;
};

simul_meta get_simul_meta(const List l, DataFrame dt)
{
  simul_meta out = {};
  out.init_year  = as<int>(l["init_year"]);
  out.age_low  = as<int>(l["ageL"]);
  out.pid = dt[as<string>(l["pids"])];
  out.year = dt[as<string>(l["years"])];
  out.age = dt[as<string>(l["ages"])];
  out.dead = dt[as<string>(l["all_cause_mrtl"])];
  out.mm_count = dt[as<string>(l["cms_count"])];
  out.mm_score = dt[as<string>(l["cms_score"])];
  return out;
}

disease_meta get_disease_meta(const List l, DataFrame dt)
{
  // l is the disease list. I.e.l$disease$chd
  disease_meta out = {};

  List incd, dgns, mrtl;


  // incidence
  if (l.containsElementNamed("incidence"))
  {
    incd = l["incidence"];
    out.incd.type =  as<string>(incd["type"]);

    if (incd.containsElementNamed("prevalence")) out.incd.prvl = dt[as<string>(incd["prevalence"])];
    if (incd.containsElementNamed("probability"))out.incd.prbl1 = dt[as<string>(incd["probability"])];

    if (incd.containsElementNamed("aggregate")) out.incd.aggregate = as<CharacterVector>(incd["aggregate"]);

    if (incd.containsElementNamed("influenced_by"))
    {
      List ib = incd["influenced_by"];
      CharacterVector tmps= ib.names();
      int n = ib.length();
      List ibb;
      if (out.incd.type == "Type0")
      {
        for (int i = 0; i < n; ++i)
        {
          ibb = ib[i];
          out.incd.influenced_by.disease_prvl.push_back(dt[as<string>(tmps[i])]);
          out.incd.influenced_by.lag.push_back(as<int>(ibb["lag"])); // set to 0 from R side
        }
      }
      else
      {
        for (int i = 0; i < n; ++i)
        {
          ibb = ib[i];
          out.incd.influenced_by.disease_prvl.push_back(dt[as<string>(tmps[i])]);
          out.incd.influenced_by.mltp.push_back(dt[as<string>(ibb["multiplier"])]);
          out.incd.influenced_by.lag.push_back(as<int>(ibb["lag"]));
        }
      }


    }
    out.incd.flag = false;
    if (incd.containsElementNamed("can_recur")) out.incd.can_recur = as<bool>(incd["can_recur"]);
    else out.incd.can_recur = false;
  }

  // diagnosis
  if (l.containsElementNamed("diagnosis"))
  {
    dgns = l["diagnosis"];

    out.dgns.type  = as<string>(dgns["type"]);
    out.dgns.mm_wt = as<double>(dgns["mm_wt"]);
    if (dgns.containsElementNamed("diagnosed")) out.dgns.prvl  = dt[as<string>(dgns["diagnosed"])];
    if (dgns.containsElementNamed("probability")) out.dgns.prbl1 = dt[as<string>(dgns["probability"])];

    if (out.dgns.type == "Type0")
    {
      List ib = dgns["influenced_by"];
      CharacterVector tmps= ib.names();
      int n = ib.length();
      List ibb;
      for (int i = 0; i < n; ++i)
      {
        ibb = ib[i];
        out.dgns.influenced_by.disease_prvl.push_back(dt[as<string>(tmps[i])]);
      }
    }

    out.dgns.flag = false;
  }

  // mortality
  if (l.containsElementNamed("mortality"))
  {
    mrtl = l["mortality"];
    out.mrtl.type = as<string>(mrtl["type"]);

    if (mrtl.containsElementNamed("probability"))
    {
      out.mrtl.prbl2 = dt[as<string>(mrtl["probability"])];
      out.mrtl1flag = false;
    }
    if (mrtl.containsElementNamed("probability1styear"))
    {
      out.mrtl.prbl1 =  dt[as<string>(mrtl["probability1styear"])];
      out.mrtl1flag = true;
    }
    if (mrtl.containsElementNamed("cure")) out.mrtl.cure = as<int>(mrtl["cure"]);
    if (mrtl.containsElementNamed("code")) out.mrtl.death_code = as<int>(mrtl["code"]);

    if (mrtl.containsElementNamed("influenced_by"))
    {
      List ib = mrtl["influenced_by"];
      CharacterVector tmps= ib.names();
      int n = ib.length();
      List ibb;
      for (int i = 0; i < n; ++i)
      {
        ibb = ib[i];
        out.mrtl.influenced_by.disease_prvl.push_back(dt[as<string>(tmps[i])]);
        out.mrtl.influenced_by.mltp.push_back(dt[as<string>(ibb["multiplier"])]);
        out.mrtl.influenced_by.lag.push_back(as<int>(ibb["lag"]));
      }
    }

    out.mrtl.flag = false;

  }


  out.seed = as<int>(l["seed"]);
  return out;
}


//' @export
// [[Rcpp::export]]
void simcpp(DataFrame dt, const List l, const int mc) {
  uint64_t seed = rng->operator()();
  rng =  dqrng::generator<pcg64>(seed);
  using parm_t = decltype(uniform)::param_type;
  uniform.param(parm_t(0.0, 1.0));

  uint64_t _stream = 0;
  uint64_t _seed = 0;

  const int n = dt.nrow();

  simul_meta meta = get_simul_meta(l, dt);

  List diseases = l["diseases"];
  const int dn = diseases.length();
  vector<disease_meta> dsmeta(dn);
  for (int j = 0; j < dn; ++j) dsmeta[j] = get_disease_meta(diseases[j], dt);


  double mltp = 1.0; // to be used for influence_by multiplier
  vector<int> tempdead; // temporary dead possibly from multiple causes
  double rn1, rn2;
  int pid_buffer = meta.pid[0]; // flag for when new pid to reset other flags. Holds the last pid
  bool pid_mrk = true;

  for (int i = 0; i < n; ++i) // loop over dt rows (and resolve all diseases for each row before you move on)
  {

    // set true when new participant (need to be here and remain true until the
    // next if

    if (meta.year[i] >= meta.init_year && meta.age[i] >= meta.age_low && meta.dead[i - 1] == 0) // if year >= init_year & alive
    {
      // NOTE values of i - x are certainly inbound and belong to the same
      // individual as long as x <= max_lag
      // NA_INTEGER == 0 is false

      if (meta.pid[i] == pid_buffer)
      {
        pid_mrk = false;
      }
      else {
        pid_mrk = true;
        pid_buffer = meta.pid[i];
      }

      _seed = u32tou64(meta.pid[i], meta.year[i]);

      for (int j = 0; j < dn; ++j) // loop over diseases
      {
        _stream = u32tou64(mc, dsmeta[j].seed);
        rng->seed(_seed, _stream);

        // Incidence ------------------------------------------------
        // Generate RN irrespective of whether will be used. Crucial for
        // reproducibility and to remove stochastic noise between scenarios
        rn1 = runif_impl();

        // reset flags for new simulants
        if (pid_mrk)
        {
          dsmeta[j].incd.flag = dsmeta[j].incd.prvl[i] > 0 ? true : false; // denotes that incd occurred
          dsmeta[j].mrtl.flag = false; // denotes cure
        }


        // reset prvl if cure and this is not a new patient
        if (dsmeta[j].mrtl.flag)
        {
          dsmeta[j].incd.prvl[i] = 0;
          // dsmeta[j].dgns.prvl[i] = 0;
          if ((dsmeta[j].dgns.type == "Type0" || dsmeta[j].dgns.type == "Type1") && dsmeta[j].dgns.mm_wt > 0.0 && dsmeta[j].dgns.prvl[i - 1] > 0) meta.mm_score[i] -= dsmeta[j].dgns.mm_wt;
          if ((dsmeta[j].dgns.type == "Type0" || dsmeta[j].dgns.type == "Type1") && dsmeta[j].dgns.mm_wt > 0.0 && dsmeta[j].dgns.prvl[i - 1] > 0) meta.mm_count[i]--;

          dsmeta[j].mrtl.flag = false;
        }

        if (dsmeta[j].incd.type == "Type0")
        {
          for (int k = 0; k < dsmeta[j].incd.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
          {
            if (dsmeta[j].incd.influenced_by.disease_prvl[k][i] > dsmeta[j].incd.prvl[i])
            {
              dsmeta[j].incd.prvl[i] = dsmeta[j].incd.influenced_by.disease_prvl[k][i];
            }
          }
        }

        if (dsmeta[j].incd.type == "Type1")
        { // NOTE Type 1 doesn't need to use flags for recurrence
          if (dsmeta[j].incd.can_recur)
          {
            if (dsmeta[j].incd.prbl1[i] == 1.0) // logic overwrites prvl for init year
              dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
          }
          else // if no recurrence assuming no cure
          {
            if (dsmeta[j].incd.prvl[i - 1] == 0 && dsmeta[j].incd.prbl1[i] == 1.0)
              dsmeta[j].incd.prvl[i] = 1;
            if (dsmeta[j].incd.prvl[i - 1] > 0)
              dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
          }
        }

        else if (dsmeta[j].incd.type == "Type2")
        {
          if (dsmeta[j].incd.can_recur) // no need to check incd.flag
          {
            if (dsmeta[j].incd.prvl[i - 1] == 0 && rn1 <= dsmeta[j].incd.prbl1[i])
            {
              dsmeta[j].incd.prvl[i] = 1;
            }

            if (dsmeta[j].mrtl.type == "Type2" || dsmeta[j].mrtl.type == "Type4")
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0 &&
                  dsmeta[j].incd.prvl[i - 1] < dsmeta[j].mrtl.cure)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
            }
            else
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
              // Above holds irrespective of if disease can be cured as
              // dsmeta[j].incd.prvl[i - 1] cannot go above duration for cure
            }

          }
          else // no recurrence
          {
            if (!dsmeta[j].incd.flag && dsmeta[j].incd.prvl[i - 1] == 0 &&
                rn1 <= dsmeta[j].incd.prbl1[i])
            {
              dsmeta[j].incd.prvl[i] = 1;
              dsmeta[j].incd.flag = true; // Set flag
            }

            if (dsmeta[j].mrtl.type == "Type2" || dsmeta[j].mrtl.type == "Type4")
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0 &&
                  dsmeta[j].incd.prvl[i - 1] < dsmeta[j].mrtl.cure)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
            }
            else
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
              // Above holds irrespective of if disease can be cured as
              // dsmeta[j].incd.prvl[i - 1] cannot go above duration for cure
            }
          }
        }

        else if (dsmeta[j].incd.type == "Type3") // Don't need this type. Can be replaced by a flag to notify disease dependence
        {
          for (int k = 0; k < dsmeta[j].incd.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
          {
            // if lag > 0 look back. For diseases depending on self, lag = 0 and that triggers the use of the .incd.flag
            if ((dsmeta[j].incd.influenced_by.lag[k] > 0 && dsmeta[j].incd.influenced_by.disease_prvl[k][i - dsmeta[j].incd.influenced_by.lag[k]] > 0) ||
                (dsmeta[j].incd.influenced_by.lag[k] == 0 && dsmeta[j].incd.flag))
            {
              mltp *= dsmeta[j].incd.influenced_by.mltp[k][i]; // no lag here
            }
          }

          if (dsmeta[j].incd.can_recur) // no need to check incd.flag
          {
            if (dsmeta[j].incd.prvl[i - 1] == 0 && rn1 <= dsmeta[j].incd.prbl1[i] * mltp)
            {
              dsmeta[j].incd.prvl[i] = 1;
            }

            if (dsmeta[j].mrtl.type == "Type2" || dsmeta[j].mrtl.type == "Type4")
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0 &&
                  dsmeta[j].incd.prvl[i - 1] < dsmeta[j].mrtl.cure)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
            }
            else
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
              // Above holds irrespective of if disease can be cured as
              // dsmeta[j].incd.prvl[i - 1] cannot go above duration for cure
            }

          }
          else // no recurrence
          {
            if (!dsmeta[j].incd.flag && dsmeta[j].incd.prvl[i - 1] == 0 &&
                rn1 <= dsmeta[j].incd.prbl1[i] * mltp)
            {
              dsmeta[j].incd.prvl[i] = 1;
              dsmeta[j].incd.flag = true; // Set flag
            }

            if (dsmeta[j].mrtl.type == "Type2" || dsmeta[j].mrtl.type == "Type4")
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0 &&
                  dsmeta[j].incd.prvl[i - 1] < dsmeta[j].mrtl.cure)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
            }
            else
            {
              if (dsmeta[j].incd.prvl[i - 1] > 0)
                dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
              // Above holds irrespective of if disease can be cured as
              // dsmeta[j].incd.prvl[i - 1] cannot go above duration for cure
            }
          }

          mltp = 1.0;
        }

        else if (dsmeta[j].incd.type == "Type4")
        {
          // TODO
        }

        else if (dsmeta[j].incd.type == "Type5")
        {
          // TODO
        }

        // diagnosis & multimorbidity --------------------------------
        rn1 = runif_impl();

        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].dgns.type == "Type0")
        {
          for (int k = 0; k < dsmeta[j].dgns.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
          {
            if (dsmeta[j].dgns.influenced_by.disease_prvl[k][i] > dsmeta[j].dgns.prvl[i])
            {
              dsmeta[j].dgns.prvl[i] = dsmeta[j].dgns.influenced_by.disease_prvl[k][i];
            }
          }
        }

        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].dgns.type == "Type1") // enter branch only for prevalent cases
        {
          if (dsmeta[j].dgns.prvl[i - 1] == 0 && rn1 <= dsmeta[j].dgns.prbl1[i])
          {
            dsmeta[j].dgns.prvl[i] = 1;
          }
          if (dsmeta[j].dgns.prvl[i - 1] > 0)
          {
            dsmeta[j].dgns.prvl[i] = dsmeta[j].dgns.prvl[i - 1] + 1;
          }
        }

        if ((dsmeta[j].dgns.type == "Type0" || dsmeta[j].dgns.type == "Type1") && dsmeta[j].dgns.prvl[i] > 0 && dsmeta[j].dgns.mm_wt > 0.0) meta.mm_score[i] += dsmeta[j].dgns.mm_wt;
        if ((dsmeta[j].dgns.type == "Type0" || dsmeta[j].dgns.type == "Type1") && dsmeta[j].dgns.prvl[i] > 0 && dsmeta[j].dgns.mm_wt > 0.0) meta.mm_count[i]++;

        // mortality ----------------------------------------------
        rn1 = runif_impl();

        if (dsmeta[j].incd.type == "Universal" || dsmeta[j].incd.prvl[i] > 0) // enter branch only for prevalent cases or Universal incidence
        {
          // Type 1 mortality (no cure, no disease dependency)
          if (dsmeta[j].mrtl.type == "Type1")
          {
            if (dsmeta[j].mrtl1flag)  // Type universal never enters this branch
            {
              if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl1[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
              if (dsmeta[j].incd.prvl[i] > 1 && rn1 < dsmeta[j].mrtl.prbl2[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
            }
            else if (rn1 < dsmeta[j].mrtl.prbl2[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
          }// End Type 1


          // Type 2 mortality (cure, no disease dependency). Not compatible
          // with universal incidence (would be meaningless as Type 2 mortality
          // allows cure).
          // NOTE cure with disease dependence is type 4
          else if (dsmeta[j].mrtl.type == "Type2")
          {
            if (dsmeta[j].incd.prvl[i] < dsmeta[j].mrtl.cure) // Valid since not universal incidence
            {
              if (dsmeta[j].mrtl1flag)  // if fatality for incident cases estimated separately
              {
                if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl1[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
                if (dsmeta[j].incd.prvl[i] > 1 && rn1 < dsmeta[j].mrtl.prbl2[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
              }
              else if (rn1 < dsmeta[j].mrtl.prbl2[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
            }
            else dsmeta[j].mrtl.flag = true; // if alive cure after defined period
          } // End Type 2 mortality

          // Type 3 mortality (no cure, disease dependency)
          else if (dsmeta[j].mrtl.type == "Type3")
          {
            for (int k = 0; k < dsmeta[j].mrtl.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
            {
              if (dsmeta[j].mrtl.influenced_by.disease_prvl[k][i - dsmeta[j].mrtl.influenced_by.lag[k]] > 0)
              {
                mltp *= dsmeta[j].mrtl.influenced_by.mltp[k][i]; // no lag here
              }
            }

            if (dsmeta[j].mrtl1flag)  // if fatality for incident cases estimated separately
            {
              if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl1[i])        tempdead.push_back(dsmeta[j].mrtl.death_code); // mltp not needed here
              if (dsmeta[j].incd.prvl[i]  > 1 && rn1 < dsmeta[j].mrtl.prbl2[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
            }
            else // if (dsmeta[j].incd.type == "Universal" || !dsmeta[j].mrtl1flag). NOTE Universal enters this branch as well.
            {
              if (rn1 < dsmeta[j].mrtl.prbl2[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
            }

            mltp = 1.0;
          } // End Type 3 mortality

          // Type 4 mortality (cure, disease dependency). Note that universal incidence is impossible to be of type 4 because there is no cure.
          else if (dsmeta[j].mrtl.type == "Type4")
          {
            for (int k = 0; k < dsmeta[j].mrtl.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
            {
              if (dsmeta[j].mrtl.influenced_by.disease_prvl[k][i - dsmeta[j].mrtl.influenced_by.lag[k]] > 0)
              {
                mltp *= dsmeta[j].mrtl.influenced_by.mltp[k][i]; // no lag here
              }
            }

            if (dsmeta[j].incd.prvl[i] < dsmeta[j].mrtl.cure) // Valid since not universal incidence
            {
              if (dsmeta[j].mrtl1flag)  // if fatality for incident cases estimated separately
              {
                if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl1[i])        tempdead.push_back(dsmeta[j].mrtl.death_code); // mltp not needed here
                if (dsmeta[j].incd.prvl[i]  > 1 && rn1 < dsmeta[j].mrtl.prbl2[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
              }
              else if (rn1 < dsmeta[j].mrtl.prbl2[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
            }
            else dsmeta[j].mrtl.flag = true; // if alive cure after defined period

            mltp = 1.0;
          } // End Type 4 mortality
        } // end if prevalent case

      } // end loop over diseases


      // resolve mortality from multiple causes in a year



      if (tempdead.size() == 1) meta.dead[i] = tempdead[0];
      if (tempdead.size() > 1)
      {
        _stream = u32tou64(mc, 1L);
        rng->seed(_seed, _stream);
        rn2 = runif_impl(); // the seed for this row defined for each mc. So is is safe/reproducible to only be calculated when needed.

        int ind = (int)(rn2 * 100000000) % tempdead.size();
        meta.dead[i] = tempdead[ind];
      }
      tempdead.clear();

    } // end year >= init_year & alive

    if ((meta.dead[i - 1] > 0 || IntegerVector::is_na(meta.dead[i - 1])) &&
        meta.year[i] >= meta.init_year && meta.age[i] >= meta.age_low)
    {
      // If dead carry forward.
      // 0 means alive, > 0 is the cod the year of death and NA is for longdeads
      // NA_INTEGER > 0 is false
      meta.dead[i] = NA_INTEGER;
    }
  } // and loop over dt rows

}
