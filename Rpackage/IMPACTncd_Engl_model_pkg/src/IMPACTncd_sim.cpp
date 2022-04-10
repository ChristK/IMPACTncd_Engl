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




// https://codescracker.com/cpp/cpp-nested-structures.htm
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
  NumericVector prbl1;
  NumericVector prbl2; // for 1st year case fatality
  infl influenced_by;
  CharacterVector aggregate;
  bool can_recur;
  int death_code;
};

struct disease_meta
{
  disease_epi incd;
  disease_epi dgns;
  disease_epi mrtl;
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
};

simul_meta get_simul_meta(const List l, DataFrame dt)
{
  simul_meta out;
  out.init_year  = as<int>(l["init_year"]);
  out.age_low  = as<int>(l["ageL"]);
  out.pid = dt[as<string>(l["pids"])];
  out.year = dt[as<string>(l["years"])];
  out.age = dt[as<string>(l["ages"])];
  out.dead = dt[as<string>(l["all_cause_mrtl"])];
  return out;
}

disease_meta get_disease_meta(const List l, DataFrame dt)
{
  // l is the disease list. I.e.l$disease$chd
  disease_meta out;

  List incd, dgns, mrtl;
  if (l.containsElementNamed("incidence")) incd = l["incidence"];
  if (l.containsElementNamed("mortality")) mrtl = l["mortality"];

  // incidence
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
    for (int i = 0; i < n; ++i)
    {
      ibb = ib[i];
      out.incd.influenced_by.disease_prvl.push_back(dt[as<string>(tmps[i])]);
      out.incd.influenced_by.mltp.push_back(dt[as<string>(ibb["multiplier"])]);
      out.incd.influenced_by.lag.push_back(as<int>(ibb["lag"]));
    }
  }

  if (incd.containsElementNamed("can_recur")) out.incd.can_recur = as<bool>(incd["can_recur"]);

  // diagnosis
  if (l.containsElementNamed("diagnosis"))
  {
    dgns = l["diagnosis"];

    out.dgns.type  = as<string>(dgns["type"]);
    out.dgns.prvl  = dt[as<string>(dgns["diagnosed"])];
    out.dgns.prbl1 = dt[as<string>(dgns["probability"])];
  }

  // mortality
  out.mrtl.type = as<string>(mrtl["type"]);

  if (mrtl.containsElementNamed("probability")) out.mrtl.prbl1 = dt[as<string>(mrtl["probability"])];
  if (mrtl.containsElementNamed("probability1styear")) out.mrtl.prbl2 =  dt[as<string>(mrtl["probability1styear"])];
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

  out.seed = as<int>(l["seed"]);
  return out;
}


//' @export
// [[Rcpp::export]]
void simcpp(DataFrame dt, const List l, const int mc) {
  uint64_t seed = rng->operator()();
  rng =  dqrng::generator<pcg64>(seed);
  using parm_t = decltype(uniform)::param_type; // can be moved into the loop
  uniform.param(parm_t(0.0, 1.0)); // can be moved into the loop

  uint64_t _stream = 0;
  uint64_t _seed = 0;

  const int n = dt.nrow();

  simul_meta meta = get_simul_meta(l, dt);

  List diseases = l["diseases"];
  const int dn = diseases.length();
  vector<disease_meta> dsmeta(dn);
  for (int j = 0; j < dn; ++j) dsmeta[j] = get_disease_meta(diseases[j], dt);


  double mltp = 1.0; // to be used for influence_by multiplier
  vector<int> tempdead; // temporary dead perhaps from multiple causes
  double rn1, rn2;
  for (int i = 0; i < n; ++i) // loop over dt rows (and resolve all diseases for each row before you move on)
  {

    if (meta.year[i] >= meta.init_year && meta.age[i] >= meta.age_low && meta.dead[i - 1] == 0) // if year >= init_year & alive
    {
      // NOTE values of i - x are certainly inbound and belong to the same
      // individual as long as x <= max_lag
      // NA_INTEGER == 0 is false

      _seed = u32tou64(meta.pid[i], meta.year[i]);

      for (int j = 0; j < dn; ++j) // loop over diseases
      {
        _stream = u32tou64(mc, dsmeta[j].seed);
        rng->seed(_seed, _stream);
        // incidence ------------------------------------------------
        // Generate RN irrespective of whether will be used. Crucial for
        // reproducibility and to remove stochastic noise between scenario
        rn1 = runif_impl();
        if (dsmeta[j].incd.type == "Type1")
        {
          // TODO
        }

        if (dsmeta[j].incd.type == "Type2")
        {
          // TODO can_recur
          if (dsmeta[j].incd.prvl[i - 1] == 0 && rn1 <= dsmeta[j].incd.prbl1[i]) dsmeta[j].incd.prvl[i] = 1;
          if (dsmeta[j].incd.prvl[i - 1] > 0) dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
        }

        if (dsmeta[j].incd.type == "Type3")
        {
          // TODO can_recur
          for (int k = 0; k < dsmeta[j].incd.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
          {
            if (dsmeta[j].incd.influenced_by.disease_prvl[k][i - dsmeta[j].incd.influenced_by.lag[k]] > 0)
            {
              mltp *= dsmeta[j].incd.influenced_by.mltp[k][i]; // no lag here
            }
          }

          if (dsmeta[j].incd.prvl[i - 1] == 0 && rn1 <= (dsmeta[j].incd.prbl1[i] * mltp)) dsmeta[j].incd.prvl[i] = 1;
          if (dsmeta[j].incd.prvl[i - 1] > 0) dsmeta[j].incd.prvl[i] = dsmeta[j].incd.prvl[i - 1] + 1;
          mltp = 1.0;
        }

        if (dsmeta[j].incd.type == "Type4")
        {
          // TODO
        }

        if (dsmeta[j].incd.type == "Type5")
        {
          // TODO
        }

        // diagnosis ----------------------------------------------
        rn1 = runif_impl();
        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].dgns.type == "Type1") // enter branch only for prevalent cases
        {
          if (dsmeta[j].dgns.prvl[i - 1] == 0 && rn1 <= dsmeta[j].dgns.prbl1[i]) dsmeta[j].dgns.prvl[i] = 1;
          if (dsmeta[j].dgns.prvl[i - 1] > 0) dsmeta[j].dgns.prvl[i] = dsmeta[j].dgns.prvl[i - 1] + 1;
        }


        // mortality ----------------------------------------------
        rn1 = runif_impl();

        if (dsmeta[j].incd.type == "Universal") // enter branch only for prevalent cases or Universal incidence
        {
          if (rn1 < dsmeta[j].mrtl.prbl1[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
        }

        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].mrtl.type == "Type1")
        {
          if (dsmeta[j].mrtl.prbl2.size() > 0)  // if fatality for incident cases estimated separately
          {
            if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl2[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
            if (dsmeta[j].incd.prvl[i] > 1 && rn1 < dsmeta[j].mrtl.prbl1[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
          }
          else
          {
            if (rn1 < dsmeta[j].mrtl.prbl1[i]) tempdead.push_back(dsmeta[j].mrtl.death_code);
          }
        }

        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].mrtl.type == "Type2")
        {
          // TODO for cancers that can be cured
        }

        if (dsmeta[j].incd.prvl[i] > 0 && dsmeta[j].mrtl.type == "Type3")
        {
          for (int k = 0; k < dsmeta[j].mrtl.influenced_by.disease_prvl.size(); ++k) // Loop over influenced by diseases
          {
            if (dsmeta[j].mrtl.influenced_by.disease_prvl[k][i - dsmeta[j].mrtl.influenced_by.lag[k]] > 0)
            {
              mltp *= dsmeta[j].mrtl.influenced_by.mltp[k][i]; // no lag here
            }
          }

          if (dsmeta[j].mrtl.prbl2.size() > 0)  // if fatality for incident cases estimated separately
          {
            if (dsmeta[j].incd.prvl[i] == 1 && rn1 < dsmeta[j].mrtl.prbl2[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
            if (dsmeta[j].incd.prvl[i] > 1 && rn1 < dsmeta[j].mrtl.prbl1[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
          }
          else
          {
            if (rn1 < dsmeta[j].mrtl.prbl1[i] * mltp) tempdead.push_back(dsmeta[j].mrtl.death_code);
          }
          mltp = 1.0;
        }

      } // end loop over diseases


      // resolve mortality from multiple causes in a year

      if (tempdead.size() == 1) meta.dead[i] = tempdead[0];
      if (tempdead.size() > 1)
      {
        _stream = u32tou64(mc, 1L); //
        rng->seed(_seed, _stream);
        rn2 = runif_impl();

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
