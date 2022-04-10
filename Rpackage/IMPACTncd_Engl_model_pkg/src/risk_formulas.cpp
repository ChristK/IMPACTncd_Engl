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
#include <math.h>
using namespace Rcpp;

// The following code was adapted from QDiabetes-2013 (http://qdiabetes.org, http://svn.clinrisk.co.uk/opensource/qdiabetes)

//' @export
// [[Rcpp::export]]
DataFrame QDiabetes(const DataFrame& df, int surv = 1)
  {
  //access the df columns
  IntegerVector age = df["age"];
  IntegerVector sex = df["sex"];
  IntegerVector b_corticosteroids = df["cst_prvl"];
  //IntegerVector b_cvd = df["b_cvd"];
  IntegerVector b_treatedhyp = df["bpmed"];
  NumericVector bmi = df["bmi"];
  IntegerVector ethrisk = df["ethnicity"];
  IntegerVector fh_diab = df["fam_t2dm"];
  IntegerVector smoke_cat = df["smoke_cat"];
  NumericVector town = df["tds"];

  const int n = df.nrows();
  NumericVector out_nocvd(n);
  NumericVector out_cvd(n);
  double a = 0.0;
  double b = 0.0;
  double dage = 0.0;
  double age_1 = 0.0;
  double age_2 = 0.0;
  double dbmi = 0.0;
  double bmi_1 = 0.0;
  double bmi_2 = 0.0;
  double dtown = 0.0;

  /* The conditional arrays */
  const double survivor_M[16] = {
    0.0,
    0.998213708400726,
    0.996353209018707,
    0.994382798671722,
    0.992213606834412,
    0.989733397960663,
    0.987064540386200,
    0.984254062175751,
    0.981255292892456,
    0.977990627288818,
    0.974455237388611,
    0.970843732357025,
    0.967315018177032,
    0.963437378406525,
    0.959633111953735,
    0.955690681934357
  };
  const double Iethrisk_M[10] = {
    0.0,
    0.0,
    1.2366090720913343000000000,
    1.4716746107789032000000000,
    1.8073235649498174000000000,
    1.2056055595936399000000000,
    0.6032369975938766100000000,
    0.9095436207452737300000000,
    0.9137604632927512900000000,
    0.7123719045990779500000000
  };
  const double Ismoke_M[5] = {
    0.0,
    0.1618238582395977700000000,
    0.1902020385619117000000000,
    0.3210636179312467100000000,
    0.4140001301797494600000000
  };
  const double survivor_F[16] = {
    0.0,
    0.998714804649353,
    0.997435748577118,
    0.996052920818329,
    0.994562506675720,
    0.992949724197388,
    0.991141080856323,
    0.989293158054352,
    0.987293541431427,
    0.985133886337280,
    0.982810735702515,
    0.980465650558472,
    0.978020071983337,
    0.975493073463440,
    0.972945988178253,
    0.970350146293640
  };

  /* The conditional arrays */

  const double Iethrisk_F[10] = {
    0.0,
    0.0,
    1.2672136244963337000000000,
    1.4277605208830098000000000,
    1.8624060798103199000000000,
    1.2379988338989651000000000,
    0.4709034172907677900000000,
    0.3476400901703160500000000,
    1.1587283467731935000000000,
    0.7335499325010315100000000
  };
  const double Ismoke_F[5] = {
    0.0,
    0.1012537024947505100000000,
    0.1915520564380613400000000,
    0.3091894136143333900000000,
    0.4646730392693820800000000
  };
  IntegerVector age_ = clamp(25, age, 84);
  NumericVector bmi_ = clamp(20, bmi, 40);
  for (int i = 0; i < n; i++)
  {
    dage = (double) age_[i];
    dbmi = bmi_[i];

    if (sex[i] == 1)
    {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage /= 10.0;
      age_1 = log(dage);
      age_2 = pow(dage, 3);
      dbmi /=10.0;
      bmi_1 = pow(dbmi, 2);
      bmi_2 = pow(dbmi, 3);

      /* Centring the continuous variables */

      age_1 += -1.496771812438965;
      age_2 += -89.149559020996094;
      bmi_1 += -6.832604885101318;
      bmi_2 += -17.859918594360352;
      dtown  = town[i] + 0.132148191332817;

      /* Start of Sum */
      a = 0.0;

      /* The conditional sums */

      a += Iethrisk_M[ethrisk[i]];
      a += Ismoke_M[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * 4.4205598323371680000000000;
      a += age_2 * -0.0041132238299394193000000;
      a += bmi_1 * 1.1169895991721528000000000;
      a += bmi_2 * -0.1793529530251269100000000;
      a += dtown * 0.0291530815903822650000000;

      /* Sum from boolean values */

      a += b_corticosteroids[i] * 0.2059811979905692400000000;
      //a += b_cvd[i] * 0.3914728454990503100000000; // a for nocvd
      a += b_treatedhyp[i] * 0.5010787979849035100000000;
      a += fh_diab[i] * 0.8385800403428993500000000;

      /* Sum from interaction terms */

      a += age_1 * bmi_1 * 0.5051031253768063500000000;
      a += age_1 * bmi_2 * -0.1375233635462656000000000;
      a += age_1 * fh_diab[i] * -1.1463560542602569000000000;
      a += age_2 * bmi_1 * -0.0015800686452772700000000;
      a += age_2 * bmi_2 * 0.0003394090057824062300000;
      a += age_2 * fh_diab[i] * 0.0018524160353981260000000;

      b = a + 0.3914728454990503100000000; // b for cvd
      /* Calculate the score itself (0 to 1)*/
      out_nocvd[i] = 1.0 - pow(survivor_M[surv], exp(a));
      out_cvd[i]   = 1.0 - pow(survivor_M[surv], exp(b));
    }
    else
    { // if (sex[i] == 2)
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage /= 10.0;
      age_1 = pow(dage, 0.5);
      age_2 = pow(dage, 3);
      dbmi /= 10.0;
      bmi_1 = dbmi;
      bmi_2 = pow(dbmi, 3);

      /* Centring the continuous variables */

      age_1 += -2.135220289230347;
      age_2 += -94.766799926757813;
      bmi_1 += -2.549620866775513;
      bmi_2 += -16.573980331420898;
      dtown  = town[i] + 0.224075347185135;

      /* Start of Sum */
      a = 0.0;

      /* The conditional sums */

      a += Iethrisk_F[ethrisk[i]];
      a += Ismoke_F[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * 4.3848331212989669000000000;
      a += age_2 * -0.0049763964406541149000000;
      a += bmi_1 * 3.3753336326064329000000000;
      a += bmi_2 * -0.0631628488667318330000000;
      a += dtown  * 0.0432726992998635970000000;

      /* Sum from boolean values */

      a += b_corticosteroids[i] * 0.2681990966241487000000000;
      //a += b_cvd * 0.3596176830984252900000000;
      a += b_treatedhyp[i] * 0.5314598436974725700000000;
      a += fh_diab[i] * 0.7315358845837640600000000;

      /* Sum from interaction terms */

      a += age_1 * bmi_1 * 1.3037832873997990000000000;
      a += age_1 * bmi_2 * -0.0708293717769046120000000;
      a += age_1 * fh_diab[i] * -0.7968266815834251800000000;
      a += age_2 * bmi_1 * -0.0067725323761278549000000;
      a += age_2 * bmi_2 * 0.0002374980728666116700000;
      a += age_2 * fh_diab[i] * 0.0017048228889394394000000;

      b = a + 0.3596176830984252900000000;
      /* Calculate the score itself */
      out_nocvd[i] = 1.0 - pow(survivor_F[surv], exp(a));
      out_cvd[i] = 1.0 - pow(survivor_F[surv], exp(b));
    }
  }
  return DataFrame::create(_["QDiabetes_nocvd"]= out_nocvd,
                           _["QDiabetes_cvd"]= out_cvd);
}

//' @export
// [[Rcpp::export]]
DataFrame QDiabetes_vec_inputs(const IntegerVector& age,
                               const IntegerVector& sex,
                               const IntegerVector& b_corticosteroids,
                               const IntegerVector& b_treatedhyp,
                               const NumericVector& bmi,
                               const IntegerVector& ethrisk,
                               const IntegerVector& fh_diab,
                               const IntegerVector& smoke_cat,
                               const NumericVector& town,
                               int surv = 1)
{


  const int n = age.size();
  NumericVector out_nocvd(n);
  NumericVector out_cvd(n);
  double a = 0.0;
  double b = 0.0;
  double dage = 0.0;
  double age_1 = 0.0;
  double age_2 = 0.0;
  double dbmi = 0.0;
  double bmi_1 = 0.0;
  double bmi_2 = 0.0;
  double dtown = 0.0;

  /* The conditional arrays */
  const double survivor_M[16] = {
    0.0,
    0.998213708400726,
    0.996353209018707,
    0.994382798671722,
    0.992213606834412,
    0.989733397960663,
    0.987064540386200,
    0.984254062175751,
    0.981255292892456,
    0.977990627288818,
    0.974455237388611,
    0.970843732357025,
    0.967315018177032,
    0.963437378406525,
    0.959633111953735,
    0.955690681934357
  };
  const double Iethrisk_M[10] = {
    0.0,
    0.0,
    1.2366090720913343000000000,
    1.4716746107789032000000000,
    1.8073235649498174000000000,
    1.2056055595936399000000000,
    0.6032369975938766100000000,
    0.9095436207452737300000000,
    0.9137604632927512900000000,
    0.7123719045990779500000000
  };
  const double Ismoke_M[5] = {
    0.0,
    0.1618238582395977700000000,
    0.1902020385619117000000000,
    0.3210636179312467100000000,
    0.4140001301797494600000000
  };
  const double survivor_F[16] = {
    0.0,
    0.998714804649353,
    0.997435748577118,
    0.996052920818329,
    0.994562506675720,
    0.992949724197388,
    0.991141080856323,
    0.989293158054352,
    0.987293541431427,
    0.985133886337280,
    0.982810735702515,
    0.980465650558472,
    0.978020071983337,
    0.975493073463440,
    0.972945988178253,
    0.970350146293640
  };

  /* The conditional arrays */

  const double Iethrisk_F[10] = {
    0.0,
    0.0,
    1.2672136244963337000000000,
    1.4277605208830098000000000,
    1.8624060798103199000000000,
    1.2379988338989651000000000,
    0.4709034172907677900000000,
    0.3476400901703160500000000,
    1.1587283467731935000000000,
    0.7335499325010315100000000
  };
  const double Ismoke_F[5] = {
    0.0,
    0.1012537024947505100000000,
    0.1915520564380613400000000,
    0.3091894136143333900000000,
    0.4646730392693820800000000
  };
  IntegerVector age_ = clamp(25, age, 84);
  NumericVector bmi_ = clamp(20, bmi, 40);
  for (int i = 0; i < n; i++)
  {
    dage = (double) age_[i];
    dbmi = bmi_[i];

    if (sex[i] == 1)
    {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage /= 10.0;
      age_1 = log(dage);
      age_2 = pow(dage, 3);
      dbmi /=10.0;
      bmi_1 = pow(dbmi, 2);
      bmi_2 = pow(dbmi, 3);

      /* Centring the continuous variables */

      age_1 += -1.496771812438965;
      age_2 += -89.149559020996094;
      bmi_1 += -6.832604885101318;
      bmi_2 += -17.859918594360352;
      dtown  = town[i] + 0.132148191332817;

      /* Start of Sum */
      a = 0.0;

      /* The conditional sums */

      a += Iethrisk_M[ethrisk[i]];
      a += Ismoke_M[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * 4.4205598323371680000000000;
      a += age_2 * -0.0041132238299394193000000;
      a += bmi_1 * 1.1169895991721528000000000;
      a += bmi_2 * -0.1793529530251269100000000;
      a += dtown * 0.0291530815903822650000000;

      /* Sum from boolean values */

      a += b_corticosteroids[i] * 0.2059811979905692400000000;
      //a += b_cvd[i] * 0.3914728454990503100000000; // a for nocvd
      a += b_treatedhyp[i] * 0.5010787979849035100000000;
      a += fh_diab[i] * 0.8385800403428993500000000;

      /* Sum from interaction terms */

      a += age_1 * bmi_1 * 0.5051031253768063500000000;
      a += age_1 * bmi_2 * -0.1375233635462656000000000;
      a += age_1 * fh_diab[i] * -1.1463560542602569000000000;
      a += age_2 * bmi_1 * -0.0015800686452772700000000;
      a += age_2 * bmi_2 * 0.0003394090057824062300000;
      a += age_2 * fh_diab[i] * 0.0018524160353981260000000;

      b = a + 0.3914728454990503100000000; // b for cvd
      /* Calculate the score itself (0 to 1)*/
      out_nocvd[i] = 1.0 - pow(survivor_M[surv], exp(a));
      out_cvd[i]   = 1.0 - pow(survivor_M[surv], exp(b));
    }
    else
    { // if (sex[i] == 2)
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage /= 10.0;
      age_1 = pow(dage, 0.5);
      age_2 = pow(dage, 3);
      dbmi /= 10.0;
      bmi_1 = dbmi;
      bmi_2 = pow(dbmi, 3);

      /* Centring the continuous variables */

      age_1 += -2.135220289230347;
      age_2 += -94.766799926757813;
      bmi_1 += -2.549620866775513;
      bmi_2 += -16.573980331420898;
      dtown  = town[i] + 0.224075347185135;

      /* Start of Sum */
      a = 0.0;

      /* The conditional sums */

      a += Iethrisk_F[ethrisk[i]];
      a += Ismoke_F[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * 4.3848331212989669000000000;
      a += age_2 * -0.0049763964406541149000000;
      a += bmi_1 * 3.3753336326064329000000000;
      a += bmi_2 * -0.0631628488667318330000000;
      a += dtown  * 0.0432726992998635970000000;

      /* Sum from boolean values */

      a += b_corticosteroids[i] * 0.2681990966241487000000000;
      //a += b_cvd * 0.3596176830984252900000000;
      a += b_treatedhyp[i] * 0.5314598436974725700000000;
      a += fh_diab[i] * 0.7315358845837640600000000;

      /* Sum from interaction terms */

      a += age_1 * bmi_1 * 1.3037832873997990000000000;
      a += age_1 * bmi_2 * -0.0708293717769046120000000;
      a += age_1 * fh_diab[i] * -0.7968266815834251800000000;
      a += age_2 * bmi_1 * -0.0067725323761278549000000;
      a += age_2 * bmi_2 * 0.0002374980728666116700000;
      a += age_2 * fh_diab[i] * 0.0017048228889394394000000;

      b = a + 0.3596176830984252900000000;
      /* Calculate the score itself */
      out_nocvd[i] = 1.0 - pow(survivor_F[surv], exp(a));
      out_cvd[i] = 1.0 - pow(survivor_F[surv], exp(b));
    }
  }
  return DataFrame::create(_["QDiabetes_nocvd"]= out_nocvd,
                           _["QDiabetes_cvd"]= out_cvd);
}



// The following code was adapted from Qrisk2-2017 (https://qrisk.org/2017/, https://qrisk.org/2017/QRISK2-2015-lgpl-source.tgz)

//' @export
// [[Rcpp::export]]
List Qrisk2 (const DataFrame& df, const bool& ignore_bmi = false,
                      const bool& ignore_sbp = false, const bool& ignore_chol = false)
{
  // TODO enable t1dm. Currently not used because diabetics are not invited in NHSHCs
  //access the df columns
  IntegerVector age = df["age"]; // 25:84
  IntegerVector sex = df["sex"]; // 1, 2
  IntegerVector b_AF = df["af_dgn"]; // 0, 1
  IntegerVector b_ra = df["ra_prvl"]; // 0, 1
  IntegerVector b_renal = df["ckd5_prvl"]; // 0, 1
  IntegerVector b_treatedhyp = df["bpmed_curr_xps"]; // 0, 1
  // b_type1 and
  IntegerVector b_type2 = df["t2dm_dgn"];
  NumericVector bmi = df["bmi_curr_xps"]; // 20:40
  IntegerVector ethrisk = df["ethnicity"]; // 1:9
  IntegerVector fh_cvd = df["famcvd"]; // 0, 1
  NumericVector rati = df["tchol_hdl_ratio"]; // 1 - 12
  NumericVector sbp = df["sbp_curr_xps"]; // 70 - 210
  IntegerVector smoke_cat = df["smoke_cat"]; // 0:4
  NumericVector town = df["tds"]; // -7 to 11

  const int n = df.nrows();
  NumericVector out(n);
  IntegerVector out_cat(n);
  const int surv = 10;

  double survivor_m[16] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.978794217109680,
    0,
    0,
    0,
    0,
    0
  };

  /* The conditional arrays */

  double Iethrisk_m[10] = {
    0,
    0,
    0.3173321430481919100000000,
    0.4738590786081115500000000,
    0.5171314655968145500000000,
    0.1370301157366419200000000,
    -0.3885522304972663900000000,
    -0.3812495485312194500000000,
    -0.4064461381650994500000000,
    -0.2285715521377336100000000
  };
  double Ismoke_m[5] = {
    0,
    0.2684479158158020200000000,
    0.6307674973877591700000000,
    0.7178078883378695700000000,
    0.8704172533465485100000000
  };


  double survivor_f[16] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0.989747583866119,
    0,
    0,
    0,
    0,
    0
  };

  /* The conditional arrays */

  double Iethrisk_f[10] = {
    0,
    0,
    0.2574099349831925900000000,
    0.6129795430571779400000000,
    0.3362159841669621300000000,
    0.1512517303224336400000000,
    -0.1794156259657768100000000,
    -0.3503423610057745400000000,
    -0.2778372483233216800000000,
    -0.1592734122665366000000000
  };
  double Ismoke_f[5] = {
    0,
    0.2119377108760385200000000,
    0.6618634379685941500000000,
    0.7570714587132305600000000,
    0.9496298251457036000000000
  };

  IntegerVector age_ = clamp(25, age, 84);
  NumericVector bmi_ = clamp(20, bmi, 40);
  NumericVector sbp_ = clamp(70, sbp, 210);
  NumericVector rati_ = clamp(1, rati, 12);
  NumericVector town_ = clamp(-7, town, 11);

  double dage = 0.0;
  double age_1 = 0.0;
  double age_2 = 0.0;
  double dbmi = 0.0;
  double bmi_2 = 0.0;
  double bmi_1 = 0.0;

  double dsbp = 0.0;
  double drati = 0.0;
  double dtown = 0.0;
  int af = 0;
  int t2dm = 0;

  for (int i = 0; i < n; i++)
  {
    if (b_AF[i] >= 1) af = 1; else af = 0;
    if (b_type2[i] >= 1) t2dm = 1; else t2dm = 0;


    if (sex[i] == 1)
    {
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */
      dage = (double)age_[i]/10.0;
      age_1 = pow(dage,-1);
      age_2 = pow(dage,2);
      dbmi=bmi_[i]/10.0;
      double bmi_2 = pow(dbmi,-2)*log(dbmi);
      double bmi_1 = pow(dbmi,-2);

      /* Centring the continuous variables */

      age_1 = age_1 - 0.233734160661697;
      age_2 = age_2 - 18.304403305053711;
      bmi_1 = bmi_1 - 0.146269768476486;
      bmi_2 = bmi_2 - 0.140587374567986;
      drati = rati_[i] - 4.321151256561279;
      dsbp  = sbp_[i] - 130.589752197265620;
      dtown = town_[i] - 0.551009356975555;

      /* Start of Sum */
      double a=0.0;

      /* The conditional sums */

      a += Iethrisk_m[ethrisk[i]];
      a += Ismoke_m[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * -18.0437312550377270000000000;
      a += age_2 * 0.0236486454254306940000000;
      if (!ignore_bmi)
      {
        a += bmi_1 * 2.5388084343581578000000000;
        a += bmi_2 * -9.1034725871528597000000000;
      }
      if (!ignore_chol) a += drati * 0.1684397636136909500000000;
      if (!ignore_sbp) a += dsbp * 0.0105003089380754820000000;
      a += dtown * 0.0323801637634487590000000;

      /* Sum from boolean values */

      a += af * 1.0363048000259454000000000;
      a += b_ra[i] * 0.2519953134791012600000000;
      a += b_renal[i] * 0.8359352886995286000000000;
      a += b_treatedhyp[i] * 0.6603459695917862600000000;
      // a += b_type1[i] * 1.3309170433446138000000000;
      a += t2dm * 0.9454348892774417900000000;
      a += fh_cvd[i] * 0.5986037897136281500000000;

      /* Sum from interaction terms */

      a += age_1 * (smoke_cat[i]==1) * 0.6186864699379683900000000;
      a += age_1 * (smoke_cat[i]==2) * 1.5522017055600055000000000;
      a += age_1 * (smoke_cat[i]==3) * 2.4407210657517648000000000;
      a += age_1 * (smoke_cat[i]==4) * 3.5140494491884624000000000;
      a += age_1 * af * 8.0382925558108482000000000;
      a += age_1 * b_renal[i] * -1.6389521229064483000000000;
      a += age_1 * b_treatedhyp[i] * 8.4621771382346651000000000;
      // a += age_1 * b_type1 * 5.4977016563835504000000000;
      a += age_1 * t2dm * 3.3974747488766690000000000;
      if (!ignore_bmi)
      {
        a += age_1 * bmi_1 * 33.8489881012767600000000000;
        a += age_1 * bmi_2 * -140.6707025404897100000000000;
      }
      a += age_1 * fh_cvd[i] * 2.0858333154353321000000000;
      if (!ignore_sbp) a += age_1 * dsbp * 0.0501283668830720540000000;
      a += age_1 * dtown * -0.1988268217186850700000000;
      a += age_2 * (smoke_cat[i]==1) * -0.0040893975066796338000000;
      a += age_2 * (smoke_cat[i]==2) * -0.0056065852346001768000000;
      a += age_2 * (smoke_cat[i]==3) * -0.0018261006189440492000000;
      a += age_2 * (smoke_cat[i]==4) * -0.0014997157296173290000000;
      a += age_2 * af * 0.0052471594895864343000000;
      a += age_2 * b_renal[i] * -0.0179663586193546390000000;
      a += age_2 * b_treatedhyp[i] * 0.0092088445323379176000000;
      // a += age_2 * b_type1 * 0.0047493510223424558000000;
      a += age_2 * t2dm * -0.0048113775783491563000000;
      if (!ignore_bmi)
      {
        a += age_2 * bmi_1 * 0.0627410757513945650000000;
        a += age_2 * bmi_2 * -0.2382914909385732100000000;
      }
      a += age_2 * fh_cvd[i] * -0.0049971149213281010000000;
      if (!ignore_sbp) a += age_2 * dsbp * -0.0000523700987951435090000;
      a += age_2 * dtown * -0.0012518116569283104000000;

      /* Calculate the score itself (0 to 1)*/
      out[i] = 1.0 - pow(survivor_m[surv], exp(a));
    }
    else
    { // if (sex[i] == 2)
      /* Applying the fractional polynomial transforms */
      /* (which includes scaling)                      */

      dage = (double)age_[i]/10.0;
      double age_2 = dage;
      double age_1 = pow(dage,.5);
      dbmi=bmi_[i]/10.0;
      double bmi_2 = pow(dbmi,-2)*log(dbmi);
      double bmi_1 = pow(dbmi,-2);

      /* Centring the continuous variables */

      age_1 = age_1 - 2.086397409439087;
      age_2 = age_2 - 4.353054523468018;
      bmi_1 = bmi_1 - 0.152244374155998;
      bmi_2 = bmi_2 - 0.143282383680344;
      drati = rati_[i] - 3.506655454635620;
      dsbp  = sbp_[i] - 125.040039062500000;
      dtown = town_[i] - 0.416743695735931;

      /* Start of Sum */
      double a=0;

      /* The conditional sums */

      a += Iethrisk_f[ethrisk[i]];
      a += Ismoke_f[smoke_cat[i]];

      /* Sum from continuous values */

      a += age_1 * 4.4417863976316578000000000;
      a += age_2 * 0.0281637210672999180000000;
      if (!ignore_bmi)
      {
        a += bmi_1 * 0.8942365304710663300000000;
        a += bmi_2 * -6.5748047596104335000000000;
      }
      if (!ignore_chol) a += drati * 0.1433900561621420900000000;
      if (!ignore_sbp) a += dsbp * 0.0128971795843613720000000;
      a += dtown * 0.0664772630011438850000000;

      /* Sum from boolean values */

      a += af * 1.6284780236484424000000000;
      a += b_ra[i] * 0.2901233104088770700000000;
      a += b_renal[i] * 1.0043796680368302000000000;
      a += b_treatedhyp[i] * 0.6180430562788129500000000;
      // a += b_type1 * 1.8400348250874599000000000;
      a += t2dm * 1.1711626412196512000000000;
      a += fh_cvd[i] * 0.5147261203665195500000000;

      /* Sum from interaction terms */

      a += age_1 * (smoke_cat[i]==1) * 0.7464406144391666500000000;
      a += age_1 * (smoke_cat[i]==2) * 0.2568541711879666600000000;
      a += age_1 * (smoke_cat[i]==3) * -1.5452226707866523000000000;
      a += age_1 * (smoke_cat[i]==4) * -1.7113013709043405000000000;
      a += age_1 * af * -7.0177986441269269000000000;
      a += age_1 * b_renal[i] * -2.9684019256454390000000000;
      a += age_1 * b_treatedhyp[i] * -4.2219906452967848000000000;
      // a += age_1 * b_type1 * 1.6835769546040080000000000;
      a += age_1 * t2dm * -2.9371798540034648000000000;
      if (!ignore_bmi)
      {
        a += age_1 * bmi_1 * 0.1797196207044682300000000;
        a += age_1 * bmi_2 * 40.2428166760658140000000000;
      }
      a += age_1 * fh_cvd[i] * 0.1439979240753906700000000;
      if (!ignore_sbp) a += age_1 * dsbp * -0.0362575233899774460000000;
      a += age_1 * dtown * 0.3735138031433442600000000;
      a += age_2 * (smoke_cat[i]==1) * -0.1927057741748231000000000;
      a += age_2 * (smoke_cat[i]==2) * -0.1526965063458932700000000;
      a += age_2 * (smoke_cat[i]==3) * 0.2313563976521429400000000;
      a += age_2 * (smoke_cat[i]==4) * 0.2307165013868296700000000;
      a += age_2 * af * 1.1395776028337732000000000;
      a += age_2 * b_renal[i] * 0.4356963208330940600000000;
      a += age_2 * b_treatedhyp[i] * 0.7265947108887239600000000;
      // a += age_2 * b_type1 * -0.6320977766275653900000000;
      a += age_2 * t2dm * 0.4023270434871086800000000;
      if (!ignore_bmi)
      {
        a += age_2 * bmi_1 * 0.1319276622711877700000000;
        a += age_2 * bmi_2 * -7.3211322435546409000000000;
      }
      a += age_2 * fh_cvd[i] * -0.1330260018273720400000000;
      if (!ignore_sbp) a += age_2 * dsbp * 0.0045842850495397955000000;
      a += age_2 * dtown * -0.0952370300845990780000000;

      /* Calculate the score itself */
      out[i] = 1.0 - pow(survivor_f[surv], exp(a));
    }
    if (out[i] <= 0.1) out_cat[i] = 1;
    else if (out[i] > 0.1 && out[i] <= 0.2) out_cat[i] = 2;
    else out_cat[i] = 3;

    CharacterVector ch = {"low","mid","high"};
    out_cat.attr("class") = "factor";
    out_cat.attr("levels") = ch;
  }
  return List::create(_["QRisk2"]= out,
                      _["Qrisk2_cat"]= out_cat);
}

