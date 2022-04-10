## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

observeEvent(input$cancel_program_checkbox_sc1, {
  if (input$cancel_program_checkbox_sc1) {
    lapply(
      c(
        "frequency_eligibility_slider_sc1",
        "age_eligibility_slider_sc1",
        "invite_known_hypertensives_checkbox_sc1",
        "invite_known_diabetics_checkbox_sc1",
        "coverage_detailed_checkbox_sc1",
        "coverage_qimd0_slider_sc1",
        "coverage_qimd1_slider_sc1",
        "coverage_qimd2_slider_sc1",
        "coverage_qimd3_slider_sc1",
        "coverage_qimd4_slider_sc1",
        "coverage_qimd5_slider_sc1",
        "coverage_cost_qimd0_sc1",
        "coverage_cost_qimd1_sc1",
        "coverage_cost_qimd2_sc1",
        "coverage_cost_qimd3_sc1",
        "coverage_cost_qimd4_sc1",
        "coverage_cost_qimd5_sc1",
        "uptake_slider_sc1",
        "uptake_cost_sc1",
        "uptake_detailed_checkbox_sc1",
        "uptake_equalprob_checkbox_sc1",
        "uptake_table_help_sc1",
        "uptake_table_sc1",
        "statin_px_slider_sc1",
        "antihtn_px_slider_sc1",
        "px_detailed_checkbox_sc1",
        "statin_px_table_help_sc1",
        "statin_px_table_sc1",
        "antihtn_px_table_help_sc1",
        "antihtn_px_table_sc1",
        "smkcess_slider_sc1",
        "smkcess_cost_sc1",
        "smkcess_cost_ovrhd_sc1",
        "wghtpct_slider_sc1",
        "wghtreduc_slider_sc1",
        "wghtloss_cost_sc1",
        "wghtloss_cost_ovrhd_sc1",
        "papct_slider_sc1",
        "papincr_slider_sc1",
        "pa_cost_sc1",
        "pa_cost_ovrhd_sc1",
        "alcoholpct_slider_sc1",
        "alcoholreduc_slider_sc1",
        "alcoholreduc_cost_sc1",
        "alcoholreduc_cost_ovrhd_sc1",
        "lifestyle_attrition_slider_sc1",
        "ignore_cholesterol_checkbox_sc1",
        "ignore_sbp_checkbox_sc1",
        "ignore_bmi_checkbox_sc1"
      ),
      hideElement,
      anim = TRUE
    )
  } else {
    lapply(
      c(
        "frequency_eligibility_slider_sc1",
        "age_eligibility_slider_sc1",
        "invite_known_hypertensives_checkbox_sc1",
        "invite_known_diabetics_checkbox_sc1",
        "coverage_detailed_checkbox_sc1",
        "smkcess_slider_sc1",
        "smkcess_cost_sc1",
        "smkcess_cost_ovrhd_sc1",
        "wghtpct_slider_sc1",
        "wghtreduc_slider_sc1",
        "wghtloss_cost_sc1",
        "wghtloss_cost_ovrhd_sc1",
        "papct_slider_sc1",
        "papincr_slider_sc1",
        "pa_cost_sc1",
        "pa_cost_ovrhd_sc1",
        "alcoholpct_slider_sc1",
        "alcoholreduc_slider_sc1",
        "alcoholreduc_cost_sc1",
        "alcoholreduc_cost_ovrhd_sc1",
        "lifestyle_attrition_slider_sc1",
        "uptake_slider_sc1",
        "uptake_cost_sc1",
        "uptake_detailed_checkbox_sc1",
        "px_detailed_checkbox_sc1",
        "ignore_cholesterol_checkbox_sc1",
        "ignore_sbp_checkbox_sc1",
        "ignore_bmi_checkbox_sc1"
      ),
      showElement,
      anim = TRUE
    )
    if (!input$coverage_detailed_checkbox_sc1) {
      lapply(
        c(
          "coverage_qimd0_slider_sc1",
          "coverage_cost_qimd0_sc1"
        ),
        showElement,
        anim = TRUE
      )
    } else {
      lapply(
        c(
          "coverage_qimd1_slider_sc1",
          "coverage_qimd2_slider_sc1",
          "coverage_qimd3_slider_sc1",
          "coverage_qimd4_slider_sc1",
          "coverage_qimd5_slider_sc1",
          "coverage_cost_qimd1_sc1",
          "coverage_cost_qimd2_sc1",
          "coverage_cost_qimd3_sc1",
          "coverage_cost_qimd4_sc1",
          "coverage_cost_qimd5_sc1"
        ),
        showElement,
        anim = TRUE
      )
    }
    if (input$uptake_detailed_checkbox_sc1) {
      lapply(
        c(
          "uptake_table_sc1",
          "uptake_table_help_sc1",
          "uptake_structural_zeroes_sc1"
        ),
        showElement,
        anim = TRUE
      )
    } else {
      # showElement(id = "uptake_equalprob_checkbox_sc1",
      #             anim = TRUE)
    }

    if (!input$px_detailed_checkbox_sc1) {
      lapply(c("statin_px_slider_sc1", "antihtn_px_slider_sc1"),
             showElement,
             anim = TRUE)
    } else {
      lapply(
        c(
          "statin_px_slider_sc1",
          "antihtn_px_slider_sc1",
          "statin_px_table_help_sc1",
          "statin_px_table_sc1",
          "antihtn_px_table_help_sc1",
          "antihtn_px_table_sc1"
        ),
        showElement,
        anim = TRUE
      )
    }
  }
})

observeEvent(input$coverage_detailed_checkbox_sc1, {
  toggleElement(id = "coverage_qimd0_slider_sc1",
                condition = !input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd0_sc1",
                condition = !input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_qimd1_slider_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd1_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_qimd2_slider_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd2_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_qimd3_slider_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd3_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_qimd4_slider_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd4_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_qimd5_slider_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "coverage_cost_qimd5_sc1",
                condition = input$coverage_detailed_checkbox_sc1, anim = TRUE)
})

observeEvent(input$uptake_detailed_checkbox_sc1, {
  toggleElement(id = "uptake_table_sc1",
                condition = input$uptake_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "uptake_table_help_sc1",
                condition = input$uptake_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "uptake_structural_zeroes_sc1",
                condition = input$uptake_detailed_checkbox_sc1, anim = TRUE)
})

observeEvent(input$px_detailed_checkbox_sc1, {
  # toggleElement(id = "statin_px_slider_sc1",
  #               condition = !input$px_detailed_checkbox_sc1, anim = TRUE)
  # toggleElement(id = "antihtn_px_slider_sc1",
  #               condition = !input$px_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "statin_px_table_help_sc1",
                condition = input$px_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "statin_px_table_sc1",
                condition = input$px_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "antihtn_px_table_help_sc1",
                condition = input$px_detailed_checkbox_sc1, anim = TRUE)
  toggleElement(id = "antihtn_px_table_sc1",
                condition = input$px_detailed_checkbox_sc1, anim = TRUE)
})
