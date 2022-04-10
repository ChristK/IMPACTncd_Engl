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

output$out_year_slider <- renderUI({
  tagList(
    sliderInput(
      "inout_year_slider",
      "Year",
      input$simulation_period_slider[[1]], input$simulation_period_slider[[2]],
      input$simulation_period_slider[[2]], 1,
      sep = "",
      ticks = FALSE
    ) %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Please drag the slider to change the end date of the scenarios and see the differences.")
      )
  )
})

hlp_frienly_names <- reactive({
  hlp_frienly_names <- list()
  for (i in seq_len(input$scenarios_number_slider)) {
    if (!input[[paste0("baseline_sc", i)]]) {
      # hlp_frienly_names[link()$n] <- link()$n
      nam <- input[[paste0("friendly_name_sc", i)]]
      if (!nam %in% hlp_frienly_names)
        hlp_frienly_names[[nam]] <- nam
    }
  }
  hlp_frienly_names
})

observe({
  updatePickerInput(session, "inout_scenario_select_cypp", choices = hlp_frienly_names(), selected = hlp_frienly_names()[[1]])
  updatePickerInput(session, "inout_scenario_select_cpp", choices = hlp_frienly_names(), selected = hlp_frienly_names()[[1]])

  updatePickerInput(session, "inout_scenario_diseases_select_cypp", choices = hlp_frienly_names(), selected = hlp_frienly_names())
  updatePickerInput(session, "inout_scenario_diseases_select_cpp", choices = hlp_frienly_names(), selected = hlp_frienly_names())
})

output$out_scenario_select <- renderUI({
  tagList(
    pickerInput(inputId = "inout_scenario_select",
                label = "Scenarios",
                choices = hlp_frienly_names(),
                selected = hlp_frienly_names(),
                options = list(`actions-box` = TRUE, `live-search` = FALSE),
                multiple = TRUE)    %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Please select the scenario(s) you would like to remove from the graphs.")
      )
  )
})

# out_proc_raw ----

out_proc_raw <- reactive(
  # input$update_output,
  {
    agegroup_filter <-
      c("30-49", "50-69", "70-89")[c("30-49", "50-69", "70-89") %in% input$out_characteristics_select]
    sex_filter <-
      c("men", "women")[c("men", "women") %in% input$out_characteristics_select]
    qimd_filter <-
      c("1 most deprived", "2", "3", "4", "5 least deprived")[c("1 most deprived", "2", "3", "4", "5 least deprived") %in% input$out_characteristics_select]
    ethn_filter <- c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )[c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    ) %in% input$out_characteristics_select]

  if (input$produce_report > 0L) {
      dt <- out_report()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    } else {
      dt <- out()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    }


    dt[, `:=` (
      # TODO use lapply and grep _cost$, and net utility
      # discounting
      net_utility = deflate(net_utility, input$out_discount_qalys_slider, year, 2019L),
      invitation_cost = deflate(invitation_cost, input$out_discount_costs_slider, year, 2019L),
      attendees_cost = deflate(attendees_cost, input$out_discount_costs_slider, year, 2019L),
      active_days_cost = deflate(active_days_cost, input$out_discount_costs_slider, year, 2019L),
      bmi_cost = deflate(bmi_cost, input$out_discount_costs_slider, year, 2019L),
      alcohol_cost = deflate(alcohol_cost, input$out_discount_costs_slider, year, 2019L),
      smoking_cost = deflate(smoking_cost, input$out_discount_costs_slider, year, 2019L),
      healthcare_cost = deflate(healthcare_cost, input$out_discount_costs_slider, year, 2019L),
      socialcare_cost = deflate(socialcare_cost, input$out_discount_costs_slider, year, 2019L),
      productivity_cost = deflate(productivity_cost, input$out_discount_costs_slider, year, 2019L),
      informal_care_cost = deflate(informal_care_cost, input$out_discount_costs_slider, year, 2019L),
      smkcess_ovrhd_cost = deflate(smkcess_ovrhd_cost, input$out_discount_costs_slider, year, 2019L),
      wghtloss_ovrhd_cost = deflate(wghtloss_ovrhd_cost, input$out_discount_costs_slider, year, 2019L),
      pa_ovrhd_cost = deflate(pa_ovrhd_cost, input$out_discount_costs_slider, year, 2019L),
      alcoholreduc_ovrhd_cost = deflate(alcoholreduc_ovrhd_cost, input$out_discount_costs_slider, year, 2019L),
      policy_cost = deflate(policy_cost, input$out_discount_costs_slider, year, 2019L),
      net_policy_cost = deflate(net_policy_cost, input$out_discount_costs_slider, year, 2019L),
      net_healthcare_cost = deflate(net_healthcare_cost, input$out_discount_costs_slider, year, 2019L),
      net_socialcare_cost = deflate(net_socialcare_cost, input$out_discount_costs_slider, year, 2019L),
      net_informal_care_cost = deflate(net_informal_care_cost, input$out_discount_costs_slider, year, 2019L),
      net_productivity_cost = deflate(net_productivity_cost, input$out_discount_costs_slider, year, 2019L),
      total_hcp_cost = deflate(total_hcp_cost, input$out_discount_costs_slider, year, 2019L),
      total_hscp_cost = deflate(total_hscp_cost, input$out_discount_costs_slider, year, 2019L),
      societal_cost = deflate(societal_cost, input$out_discount_costs_slider, year, 2019L)
    )
    ][,
      `:=` (
        invitation_cost_cml = round(cumsum(invitation_cost)),
        attendees_cost_cml = round(cumsum(attendees_cost)),
        active_days_cost_cml = round(cumsum(active_days_cost)),
        bmi_cost_cml = round(cumsum(bmi_cost)),
        alcohol_cost_cml = round(cumsum(alcohol_cost)),
        smoking_cost_cml = round(cumsum(smoking_cost)),

        cvd_incd_cml = round(cumsum(cvd_incd)),
        chd_incd_cml = round(cumsum(chd_incd)),
        stroke_incd_cml = round(cumsum(stroke_incd)),
        poststroke_dementia_incd_cml = round(cumsum(poststroke_dementia_incd)),
        copd_incd_cml = round(cumsum(copd_incd)),
        t2dm_incd_cml = round(cumsum(t2dm_incd)),
        lung_ca_incd_cml = round(cumsum(lung_ca_incd)),
        colon_ca_incd_cml = round(cumsum(colon_ca_incd)),
        breast_ca_incd_cml = round(cumsum(breast_ca_incd)),
        # af_incd_cml = round(cumsum(af_incd)),
        # htn_incd_cml = round(cumsum(htn_incd)),
        pops_cml = round(cumsum(pops)),

        cvd_prvl_cml = round(cumsum(cvd_prvl)),
        chd_prvl_cml = round(cumsum(chd_prvl)),
        stroke_prvl_cml = round(cumsum(stroke_prvl)),
        poststroke_dementia_prvl_cml = round(cumsum(poststroke_dementia_prvl)),
        copd_prvl_cml = round(cumsum(copd_prvl)),
        t2dm_prvl_cml = round(cumsum(t2dm_prvl)),
        lung_ca_prvl_cml = round(cumsum(lung_ca_prvl)),
        colon_ca_prvl_cml = round(cumsum(colon_ca_prvl)),
        breast_ca_prvl_cml = round(cumsum(breast_ca_prvl)),

        chd_dgn_cml = round(cumsum(chd_dgn)),
        stroke_dgn_cml = round(cumsum(stroke_dgn)),
        poststroke_dementia_dgn_cml = round(cumsum(poststroke_dementia_dgn)),
        copd_dgn_cml = round(cumsum(copd_dgn)),
        t2dm_dgn_cml = round(cumsum(t2dm_dgn)),
        lung_ca_dgn_cml = round(cumsum(lung_ca_dgn)),
        colon_ca_dgn_cml = round(cumsum(colon_ca_dgn)),
        breast_ca_dgn_cml = round(cumsum(breast_ca_dgn)),

        chd_mrtl_cml = round(cumsum(chd_mrtl)),
        stroke_mrtl_cml = round(cumsum(stroke_mrtl)),
        copd_mrtl_cml = round(cumsum(copd_mrtl)),
        lung_ca_mrtl_cml = round(cumsum(lung_ca_mrtl)),
        colon_ca_mrtl_cml = round(cumsum(colon_ca_mrtl)),
        breast_ca_mrtl_cml = round(cumsum(breast_ca_mrtl)),
        nonmodelled_mrtl_cml = round(cumsum(nonmodelled_mrtl)),
        all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)),

        eq5d_cml = cumsum(eq5d),
        healthcare_cost_cml = round(cumsum(healthcare_cost)),
        socialcare_cost_cml = round(cumsum(socialcare_cost)),
        productivity_cost_cml = round(cumsum(productivity_cost)),
        informal_care_cost_cml = round(cumsum(informal_care_cost)),

        smkcess_ovrhd_cost_cml = round(cumsum(smkcess_ovrhd_cost)),
        wghtloss_ovrhd_cost_cml = round(cumsum(wghtloss_ovrhd_cost)),
        pa_ovrhd_cost_cml = round(cumsum(pa_ovrhd_cost)),
        alcoholreduc_ovrhd_cost_cml = round(cumsum(alcoholreduc_ovrhd_cost)),
        policy_cost_cml = round(cumsum(policy_cost)),
        net_utility_cml = signif(cumsum(net_utility), 2L),
        net_policy_cost_cml = round(cumsum(net_policy_cost)),
        net_healthcare_cost_cml = round(cumsum(net_healthcare_cost)),
        net_socialcare_cost_cml = round(cumsum(net_socialcare_cost)),
        net_informal_care_cost_cml = round(cumsum(net_informal_care_cost)),
        net_productivity_cost_cml = round(cumsum(net_productivity_cost)),

        cpp_cvd_cml = round(cumsum(cpp_cvd)),
        cpp_chd_cml = round(cumsum(cpp_chd)),
        cpp_stroke_cml = round(cumsum(cpp_stroke)),
        cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
        # cpp_htn_cml = round(cumsum(cpp_htn)),
        cpp_copd_cml = round(cumsum(cpp_copd)),
        # cpp_af_cml = round(cumsum(cpp_af)),
        cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
        cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
        cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
        cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),

        cypp_cvd_cml = round(cumsum(cypp_cvd)),
        cypp_chd_cml = round(cumsum(cypp_chd)),
        cypp_stroke_cml = round(cumsum(cypp_stroke)),
        cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
        # cypp_htn_cml = round(cumsum(cypp_htn)),
        # cypp_af_cml = round(cumsum(cypp_af)),
        cypp_copd_cml = round(cumsum(cypp_copd)),
        cypp_t2dm_cml = round(cumsum(cypp_t2dm)),
        cypp_lung_ca_cml = round(cumsum(cypp_lung_ca)),
        cypp_colon_ca_cml = round(cumsum(cypp_colon_ca)),
        cypp_breast_ca_cml = round(cumsum(cypp_breast_ca)),

        dpp_nonmodelled_cml = round(cumsum(dpp_nonmodelled)),
        dpp_chd_cml = round(cumsum(dpp_chd)),
        dpp_stroke_cml = round(cumsum(dpp_stroke)),
        dpp_copd_cml = round(cumsum(dpp_copd)),
        dpp_lung_ca_cml = round(cumsum(dpp_lung_ca)),
        dpp_colon_ca_cml = round(cumsum(dpp_colon_ca)),
        dpp_breast_ca_cml = round(cumsum(dpp_breast_ca)),
        dpp_all_cause_cml = round(cumsum(dpp_all_cause)),

        total_hcp_cost_cml = round(cumsum(total_hcp_cost)),
        total_hscp_cost_cml = round(cumsum(total_hscp_cost)),
        societal_cost_cml = round(cumsum(societal_cost))
      ),
      by = setdiff(strata_for_gui, "year") # strata_for_gui
    ][, nmb_cml := net_monetary_benefit_cml(.SD, input$health_econ_perspective_checkbox, input$out_wtp_box)
    ]
    # CBR cannot be summed because it is a ratio so it is inappropriate to put here
  })

# out_proc ----
out_proc <- reactive(
  # input$update_output,
  {
    agegroup_filter <-
      c("30-49", "50-69", "70-89")[c("30-49", "50-69", "70-89") %in% input$out_characteristics_select]
    sex_filter <-
      c("men", "women")[c("men", "women") %in% input$out_characteristics_select]
    qimd_filter <-
      c("1 most deprived", "2", "3", "4", "5 least deprived")[c("1 most deprived", "2", "3", "4", "5 least deprived") %in% input$out_characteristics_select]
    ethn_filter <- c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )[c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    ) %in% input$out_characteristics_select]

    if (input$produce_report > 0L) {
      dt <- out_report()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    } else {
      dt <- out()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    }
    dt <- sum_dt(dt)

    # Discounting costs and utility
    nam <- grep("_utility$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_qalys_slider, dt$year, 2019L))
    }
    nam <- grep("_cost$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_costs_slider, dt$year, 2019L))
    }

    setkey(dt, year, friendly_name, mc)
    dt[,
       `:=` (
         invitation_cost_cml = round(cumsum(invitation_cost)),
         attendees_cost_cml = round(cumsum(attendees_cost)),
         active_days_cost_cml = round(cumsum(active_days_cost)),
         bmi_cost_cml = round(cumsum(bmi_cost)),
         alcohol_cost_cml = round(cumsum(alcohol_cost)),
         smoking_cost_cml = round(cumsum(smoking_cost)),

         cvd_incd_cml = round(cumsum(cvd_incd)),
         chd_incd_cml = round(cumsum(chd_incd)),
         stroke_incd_cml = round(cumsum(stroke_incd)),
         poststroke_dementia_incd_cml = round(cumsum(poststroke_dementia_incd)),
         copd_incd_cml = round(cumsum(copd_incd)),
         t2dm_incd_cml = round(cumsum(t2dm_incd)),
         lung_ca_incd_cml = round(cumsum(lung_ca_incd)),
         colon_ca_incd_cml = round(cumsum(colon_ca_incd)),
         breast_ca_incd_cml = round(cumsum(breast_ca_incd)),
         # af_incd_cml = round(cumsum(af_incd)),
         # htn_incd_cml = round(cumsum(htn_incd)),
         pops_cml = round(cumsum(pops)),

         cvd_prvl_cml = round(cumsum(cvd_prvl)),
         chd_prvl_cml = round(cumsum(chd_prvl)),
         stroke_prvl_cml = round(cumsum(stroke_prvl)),
         poststroke_dementia_prvl_cml = round(cumsum(poststroke_dementia_prvl)),
         copd_prvl_cml = round(cumsum(copd_prvl)),
         t2dm_prvl_cml = round(cumsum(t2dm_prvl)),
         lung_ca_prvl_cml = round(cumsum(lung_ca_prvl)),
         colon_ca_prvl_cml = round(cumsum(colon_ca_prvl)),
         breast_ca_prvl_cml = round(cumsum(breast_ca_prvl)),

         chd_dgn_cml = round(cumsum(chd_dgn)),
         stroke_dgn_cml = round(cumsum(stroke_dgn)),
         poststroke_dementia_dgn_cml = round(cumsum(poststroke_dementia_dgn)),
         copd_dgn_cml = round(cumsum(copd_dgn)),
         t2dm_dgn_cml = round(cumsum(t2dm_dgn)),
         lung_ca_dgn_cml = round(cumsum(lung_ca_dgn)),
         colon_ca_dgn_cml = round(cumsum(colon_ca_dgn)),
         breast_ca_dgn_cml = round(cumsum(breast_ca_dgn)),

         chd_mrtl_cml = round(cumsum(chd_mrtl)),
         stroke_mrtl_cml = round(cumsum(stroke_mrtl)),
         copd_mrtl_cml = round(cumsum(copd_mrtl)),
         lung_ca_mrtl_cml = round(cumsum(lung_ca_mrtl)),
         colon_ca_mrtl_cml = round(cumsum(colon_ca_mrtl)),
         breast_ca_mrtl_cml = round(cumsum(breast_ca_mrtl)),
         nonmodelled_mrtl_cml = round(cumsum(nonmodelled_mrtl)),
         all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)),

         eq5d_cml = cumsum(eq5d),
         healthcare_cost_cml = round(cumsum(healthcare_cost)),
         socialcare_cost_cml = round(cumsum(socialcare_cost)),
         productivity_cost_cml = round(cumsum(productivity_cost)),
         informal_care_cost_cml = round(cumsum(informal_care_cost)),

         smkcess_ovrhd_cost_cml = round(cumsum(smkcess_ovrhd_cost)),
         wghtloss_ovrhd_cost_cml = round(cumsum(wghtloss_ovrhd_cost)),
         pa_ovrhd_cost_cml = round(cumsum(pa_ovrhd_cost)),
         alcoholreduc_ovrhd_cost_cml = round(cumsum(alcoholreduc_ovrhd_cost)),
         policy_cost_cml = round(cumsum(policy_cost)),
         net_utility_cml = signif(cumsum(net_utility), 2L),
         net_policy_cost_cml = round(cumsum(net_policy_cost)),
         net_healthcare_cost_cml = round(cumsum(net_healthcare_cost)),
         net_socialcare_cost_cml = round(cumsum(net_socialcare_cost)),
         net_informal_care_cost_cml = round(cumsum(net_informal_care_cost)),
         net_productivity_cost_cml = round(cumsum(net_productivity_cost)),

         cpp_cvd_cml = round(cumsum(cpp_cvd)),
         cpp_chd_cml = round(cumsum(cpp_chd)),
         cpp_stroke_cml = round(cumsum(cpp_stroke)),
         cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
         # cpp_htn_cml = round(cumsum(cpp_htn)),
         cpp_copd_cml = round(cumsum(cpp_copd)),
         # cpp_af_cml = round(cumsum(cpp_af)),
         cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
         cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
         cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
         cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),

         cypp_cvd_cml = round(cumsum(cypp_cvd)),
         cypp_chd_cml = round(cumsum(cypp_chd)),
         cypp_stroke_cml = round(cumsum(cypp_stroke)),
         cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
         # cypp_htn_cml = round(cumsum(cypp_htn)),
         # cypp_af_cml = round(cumsum(cypp_af)),
         cypp_copd_cml = round(cumsum(cypp_copd)),
         cypp_t2dm_cml = round(cumsum(cypp_t2dm)),
         cypp_lung_ca_cml = round(cumsum(cypp_lung_ca)),
         cypp_colon_ca_cml = round(cumsum(cypp_colon_ca)),
         cypp_breast_ca_cml = round(cumsum(cypp_breast_ca)),

         dpp_nonmodelled_cml = round(cumsum(dpp_nonmodelled)),
         dpp_chd_cml = round(cumsum(dpp_chd)),
         dpp_stroke_cml = round(cumsum(dpp_stroke)),
         dpp_copd_cml = round(cumsum(dpp_copd)),
         dpp_lung_ca_cml = round(cumsum(dpp_lung_ca)),
         dpp_colon_ca_cml = round(cumsum(dpp_colon_ca)),
         dpp_breast_ca_cml = round(cumsum(dpp_breast_ca)),
         dpp_all_cause_cml = round(cumsum(dpp_all_cause)),

         total_hcp_cost_cml = round(cumsum(total_hcp_cost)),
         total_hscp_cost_cml = round(cumsum(total_hscp_cost)),
         societal_cost_cml = round(cumsum(societal_cost))       ),
       by = c("friendly_name", "mc")
    ][, nmb_cml := net_monetary_benefit_cml(.SD, input$health_econ_perspective_checkbox, input$out_wtp_box)
    ]
    # CBR cannot be summed because it is a ratio so it is inappropriate to put here
  })

# out_proc_qimd ----
out_proc_qimd <- reactive(
  # input$update_output,
  {
    agegroup_filter <-
      c("30-49", "50-69", "70-89")[c("30-49", "50-69", "70-89") %in% input$out_characteristics_select]
    sex_filter <-
      c("men", "women")[c("men", "women") %in% input$out_characteristics_select]
    qimd_filter <-
      c("1 most deprived", "2", "3", "4", "5 least deprived")[c("1 most deprived", "2", "3", "4", "5 least deprived") %in% input$out_characteristics_select]
    ethn_filter <- c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )[c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    ) %in% input$out_characteristics_select]

    if (input$produce_report > 0L) {
      dt <- out_report()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    } else {
      dt <- out()[year <= input$inout_year_slider &
          friendly_name %in% input$inout_scenario_select &
          agegrp %in% agegroup_filter &
          sex %in% sex_filter &
          qimd %in% qimd_filter &
          ethnicity %in% ethn_filter , ]
    }

    dt <- sum_dt(dt, c("year", "friendly_name", "mc", "qimd"))


    # Discounting costs and utility
    nam <- grep("_utility$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_qalys_slider, dt$year, 2019L))
    }
    nam <- grep("_cost$", names(dt), value = TRUE)
    for (j in nam) {
      set(dt, NULL, j, deflate(dt[[j]], input$out_discount_costs_slider, dt$year, 2019L))
    }

    setkey(dt, year, friendly_name, qimd, mc)
    dt[,
       `:=` (
         invitation_cost_cml = round(cumsum(invitation_cost)),
         attendees_cost_cml = round(cumsum(attendees_cost)),
         active_days_cost_cml = round(cumsum(active_days_cost)),
         bmi_cost_cml = round(cumsum(bmi_cost)),
         alcohol_cost_cml = round(cumsum(alcohol_cost)),
         smoking_cost_cml = round(cumsum(smoking_cost)),

         cvd_incd_cml = round(cumsum(cvd_incd)),
         chd_incd_cml = round(cumsum(chd_incd)),
         stroke_incd_cml = round(cumsum(stroke_incd)),
         poststroke_dementia_incd_cml = round(cumsum(poststroke_dementia_incd)),
         copd_incd_cml = round(cumsum(copd_incd)),
         t2dm_incd_cml = round(cumsum(t2dm_incd)),
         lung_ca_incd_cml = round(cumsum(lung_ca_incd)),
         colon_ca_incd_cml = round(cumsum(colon_ca_incd)),
         breast_ca_incd_cml = round(cumsum(breast_ca_incd)),
         # af_incd_cml = round(cumsum(af_incd)),
         # htn_incd_cml = round(cumsum(htn_incd)),
         pops_cml = round(cumsum(pops)),

         cvd_prvl_cml = round(cumsum(cvd_prvl)),
         chd_prvl_cml = round(cumsum(chd_prvl)),
         stroke_prvl_cml = round(cumsum(stroke_prvl)),
         poststroke_dementia_prvl_cml = round(cumsum(poststroke_dementia_prvl)),
         copd_prvl_cml = round(cumsum(copd_prvl)),
         t2dm_prvl_cml = round(cumsum(t2dm_prvl)),
         lung_ca_prvl_cml = round(cumsum(lung_ca_prvl)),
         colon_ca_prvl_cml = round(cumsum(colon_ca_prvl)),
         breast_ca_prvl_cml = round(cumsum(breast_ca_prvl)),

         chd_dgn_cml = round(cumsum(chd_dgn)),
         stroke_dgn_cml = round(cumsum(stroke_dgn)),
         poststroke_dementia_dgn_cml = round(cumsum(poststroke_dementia_dgn)),
         copd_dgn_cml = round(cumsum(copd_dgn)),
         t2dm_dgn_cml = round(cumsum(t2dm_dgn)),
         lung_ca_dgn_cml = round(cumsum(lung_ca_dgn)),
         colon_ca_dgn_cml = round(cumsum(colon_ca_dgn)),
         breast_ca_dgn_cml = round(cumsum(breast_ca_dgn)),

         chd_mrtl_cml = round(cumsum(chd_mrtl)),
         stroke_mrtl_cml = round(cumsum(stroke_mrtl)),
         copd_mrtl_cml = round(cumsum(copd_mrtl)),
         lung_ca_mrtl_cml = round(cumsum(lung_ca_mrtl)),
         colon_ca_mrtl_cml = round(cumsum(colon_ca_mrtl)),
         breast_ca_mrtl_cml = round(cumsum(breast_ca_mrtl)),
         nonmodelled_mrtl_cml = round(cumsum(nonmodelled_mrtl)),
         all_cause_mrtl_cml = round(cumsum(all_cause_mrtl)),

         eq5d_cml = cumsum(eq5d),
         healthcare_cost_cml = round(cumsum(healthcare_cost)),
         socialcare_cost_cml = round(cumsum(socialcare_cost)),
         productivity_cost_cml = round(cumsum(productivity_cost)),
         informal_care_cost_cml = round(cumsum(informal_care_cost)),

         smkcess_ovrhd_cost_cml = round(cumsum(smkcess_ovrhd_cost)),
         wghtloss_ovrhd_cost_cml = round(cumsum(wghtloss_ovrhd_cost)),
         pa_ovrhd_cost_cml = round(cumsum(pa_ovrhd_cost)),
         alcoholreduc_ovrhd_cost_cml = round(cumsum(alcoholreduc_ovrhd_cost)),
         policy_cost_cml = round(cumsum(policy_cost)),
         net_utility_cml = signif(cumsum(net_utility), 2L),
         net_policy_cost_cml = round(cumsum(net_policy_cost)),
         net_healthcare_cost_cml = round(cumsum(net_healthcare_cost)),
         net_socialcare_cost_cml = round(cumsum(net_socialcare_cost)),
         net_informal_care_cost_cml = round(cumsum(net_informal_care_cost)),
         net_productivity_cost_cml = round(cumsum(net_productivity_cost)),

         cpp_cvd_cml = round(cumsum(cpp_cvd)),
         cpp_chd_cml = round(cumsum(cpp_chd)),
         cpp_stroke_cml = round(cumsum(cpp_stroke)),
         cpp_poststroke_dementia_cml = round(cumsum(cpp_poststroke_dementia)),
         # cpp_htn_cml = round(cumsum(cpp_htn)),
         cpp_copd_cml = round(cumsum(cpp_copd)),
         # cpp_af_cml = round(cumsum(cpp_af)),
         cpp_t2dm_cml = round(cumsum(cpp_t2dm)),
         cpp_lung_ca_cml = round(cumsum(cpp_lung_ca)),
         cpp_colon_ca_cml = round(cumsum(cpp_colon_ca)),
         cpp_breast_ca_cml = round(cumsum(cpp_breast_ca)),

         cypp_cvd_cml = round(cumsum(cypp_cvd)),
         cypp_chd_cml = round(cumsum(cypp_chd)),
         cypp_stroke_cml = round(cumsum(cypp_stroke)),
         cypp_poststroke_dementia_cml = round(cumsum(cypp_poststroke_dementia)),
         # cypp_htn_cml = round(cumsum(cypp_htn)),
         # cypp_af_cml = round(cumsum(cypp_af)),
         cypp_copd_cml = round(cumsum(cypp_copd)),
         cypp_t2dm_cml = round(cumsum(cypp_t2dm)),
         cypp_lung_ca_cml = round(cumsum(cypp_lung_ca)),
         cypp_colon_ca_cml = round(cumsum(cypp_colon_ca)),
         cypp_breast_ca_cml = round(cumsum(cypp_breast_ca)),

         dpp_nonmodelled_cml = round(cumsum(dpp_nonmodelled)),
         dpp_chd_cml = round(cumsum(dpp_chd)),
         dpp_stroke_cml = round(cumsum(dpp_stroke)),
         dpp_copd_cml = round(cumsum(dpp_copd)),
         dpp_lung_ca_cml = round(cumsum(dpp_lung_ca)),
         dpp_colon_ca_cml = round(cumsum(dpp_colon_ca)),
         dpp_breast_ca_cml = round(cumsum(dpp_breast_ca)),
         dpp_all_cause_cml = round(cumsum(dpp_all_cause)),

         total_hcp_cost_cml = round(cumsum(total_hcp_cost)),
         total_hscp_cost_cml = round(cumsum(total_hscp_cost)),
         societal_cost_cml = round(cumsum(societal_cost))
       ),
       by = c("friendly_name", "mc", "qimd")
    ][, nmb_cml := net_monetary_benefit_cml(.SD, input$health_econ_perspective_checkbox, input$out_wtp_box)
    ]
    # CBR cannot be summed because it is a ratio so it is inappropriate to put here
  })

out_summary <- reactive({
  strata <- c("year", "friendly_name", "mc")
  dt     <-  melt(out_proc(), strata)
  setkey(dt, variable)
  probabilities <- c(0.5, 0.025, 0.975, 0.1, 0.9)
  if (nrow(dt) > 0) { # for nrow == 0 fquantile crashes
  dt <- dt[, fquantile_byid(value, probabilities, variable, TRUE),
           by = c("year", "friendly_name")]
  setnames(dt, c("year", "friendly_name", paste0("V", 1:6)),
           c("Year", "Scenario", "Output", "Median", "2.5% UI", "97.5% UI",
             "10% UI", "90% UI"))
  }
  dt
})

# CE plane ----
output$cep1_1 <- renderPlotly({

  if (input$health_econ_perspective_checkbox == "Societal perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = societal_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Health and social care perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = total_hscp_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Healthcare perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = total_hcp_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  }

  # factor(friendly_name) otherwise levels include baseline scenario
  # [, sum_dt(.SD, c("mc", "friendly_name"), character(0))]

  max_x <- tt[, max(abs(net_utility_cml))] * 1.2
  wtp_thres <- reactive(max_x * input$out_wtp_box)
  max_y <-  max(tt[, max(abs(cost_cml))] * 1.2, wtp_thres())
  trng_path <- paste0("M 0 0 L ",  max_x, " ", wtp_thres(), " L ", max_x, " 0 Z")


  # TODO separate this code from this specific graph because it is universal. Move it somewhere it is obviously universal
  # TODO synchronise colours and graphs with the scenario selection on the left side bar

  if (input$res_display_cep1_1) tt <- median_dt(tt, "friendly_name", "mc", digits = 1)

  p <-
    plot_ly(
      tt,
      x = ~ net_utility_cml,
      y = ~ cost_cml,
      color = ~ friendly_name,
      colors =  colsymb()$colour,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )

  p <-
    layout(
      p,
      yaxis = list(title = "Incremental cumulative cost (£)"),
      xaxis = list(title = "Incremental cumulative effects (QALYs)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y"),
        list(type = "path",
             fillcolor = "blue", line = list(color = "blue"), opacity = 0.2,
             layer = "below",
             path = trng_path),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.2,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = max_y, y1 = -max_y, yref = "y"),
        list(type = "line",
          line = list(color = "black", width = 1, dash = "dash"),
          x0 = max_x, x1 = -max_x,
          y0 = wtp_thres(), y1 = -wtp_thres())
      ))

  # p <- animation_opts(p, frame = 1000, redraw = FALSE)
  # p <- animation_slider(p,
  #                       currentvalue = list(prefix = "Year: ",
  #                                           font = list(color = "red")))

})

output$cep1 <- renderPlotly({

  if (input$health_econ_perspective_checkbox == "Societal perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = societal_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Health and social care perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = total_hscp_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Healthcare perspective") {
    tt <- out_proc()[year == max(year),
      .(net_utility_cml, cost_cml = total_hcp_cost_cml, mc,
        friendly_name = factor(friendly_name))]
  }
  # [, sum_dt(.SD, c("mc", "friendly_name"), character(0))]
  max_x <- tt[, max(abs(net_utility_cml))] * 1.2
  wtp_thres <- reactive(max_x * input$out_wtp_box)
  max_y <-  max(tt[, max(abs(cost_cml))] * 1.2, wtp_thres())
  trng_path <- paste0("M 0 0 L ",  max_x, " ", wtp_thres(), " L ", max_x, " 0 Z")


  # TODO separate this code from this specific graph because it is universal.
  # Move it somewhere it is obviously universal
  # TODO synchronise colours and graphs with the scenario selection on the left
  # side bar
  if (input$res_display_cep1) tt <- median_dt(tt, "friendly_name", "mc", digits = 1)


  p <-
    plot_ly(
      tt,
      x = ~ net_utility_cml,
      y = ~ cost_cml,
      color = ~ friendly_name,
      colors =  colsymb()$colour,
      # frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )

  p <-
    layout(
      p,
      yaxis = list(title = "Incremental cumulative cost (£)"),
      xaxis = list(title = "Incremental cumulative effects (QALYs)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y"),
        list(type = "path",
             fillcolor = "blue", line = list(color = "blue"), opacity = 0.2,
             layer = "below",
             path = trng_path),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.2,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = max_y, y1 = -max_y, yref = "y"),
        list(type = "line",
          line = list(color = "black", width = 1, dash = "dash"),
          x0 = max_x, x1 = -max_x,
          y0 = wtp_thres(), y1 = -wtp_thres())
      ))

  # p <- animation_opts(p, frame = 1000, redraw = FALSE)
  # p <- animation_slider(p,
  #                       currentvalue = list(prefix = "Year: ",
  #                                           font = list(color = "red")))

})

output$cep_anim <- renderPlotly({

  if (input$health_econ_perspective_checkbox == "Societal perspective") {
    tt <- out_proc()[,
      .(net_utility_cml, cost_cml = societal_cost_cml, mc, year,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Health and social care perspective") {
    tt <- out_proc()[,
      .(net_utility_cml, cost_cml = total_hscp_cost_cml, mc, year,
        friendly_name = factor(friendly_name))]
  } else if (input$health_econ_perspective_checkbox == "Healthcare perspective") {
    tt <- out_proc()[,
      .(net_utility_cml, cost_cml = total_hcp_cost_cml, mc, year,
        friendly_name = factor(friendly_name))]
  }

  max_x <- tt[, max(abs(net_utility_cml))] * 1.2
  wtp_thres <- reactive(max_x * input$out_wtp_box)
  max_y <-  max(tt[, max(abs(cost_cml))] * 1.2, wtp_thres())
  trng_path <- paste0("M 0 0 L ",  max_x, " ", wtp_thres(), " L ", max_x, " 0 Z")



  p <-
    plot_ly(
      tt,
      x = ~ net_utility_cml,
      y = ~ cost_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Incremental cumulative cost (£)"),
      xaxis = list(title = "Incremental cumulative effects (QALYs)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y"),
        list(type = "path",
             fillcolor = "blue", line = list(color = "blue"), opacity = 0.2,
             layer = "below",
             path = trng_path),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.2,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = max_y, y1 = -max_y, yref = "y"),
        list(type = "line",
          line = list(color = "black", width = 1, dash = "dash"),
          x0 = max_x, x1 = -max_x,
          y0 = wtp_thres(), y1 = -wtp_thres())
      ))

  # NOTE gives warning when colours and symbols are defined. Known bug
  # https://github.com/ropensci/plotly/issues/1696

  # TODO the warning message concerns the colours and symbols that are not added
  # on every frames of the animation as they should be. They are only added on
  # the first. This needs to be done 3 times, for the 3 graphs with animations
  p <- animation_opts(p, frame = 1000, redraw = FALSE)
  p <- animation_slider(p,
                        currentvalue = list(prefix = "Year: ",
                                            font = list(color = "red")))
})

output$cep_p_ce <- renderPlotly({
  tt <-
    out_proc()[, .(nmb_cml, # nmb updates with healthecon perspective
      mc,
      friendly_name = factor(friendly_name),
      year)][, .(prop_if(nmb_cml > 0)), by = .(friendly_name, year)][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]

  # [, sum_dt(.SD, c("mc", "friendly_name", "year"), character(0))]
  plot_ly(
    tt,
    x = ~ year,
    y = ~ V2,
    type = "scatter",
    mode = "lines+markers",
    color = ~ friendly_name,
    colors = colsymb()$colour,
    symbol = ~ friendly_name,
    symbols = colsymb()$symbol,
    line = list(shape = "spline", smoothing = 1.3)
  ) %>%

    add_lines(
      x = ~ year,
      y = input$decision_aid_gui,
      name = "Decision aid",
      color = NULL,
      symbol = NULL,
      line = list(color = "black", dash = "dot")
    ) %>%
    layout(
      yaxis = list(
        title = "Probability of cost-effective policy",
        range = c(-0.05, 1.05),
        tickformat = ",.0%"
      ),
      xaxis = list(title = "Year")
    )
})

output$cep_p_cs <- renderPlotly({
  if (input$health_econ_perspective_checkbox == "Societal perspective") {
    tt <-
      out_proc()[, .(societal_cost_cml,
        mc,
        friendly_name = factor(friendly_name),
        year)][, .(prop_if(societal_cost_cml <= 0)), by = .(friendly_name, year)][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]
  } else if (input$health_econ_perspective_checkbox == "Health and social care perspective") {
    tt <-
      out_proc()[, .(total_hscp_cost_cml,
        mc,
        friendly_name = factor(friendly_name),
        year)][, .(prop_if(total_hscp_cost_cml <= 0)), by = .(friendly_name, year)][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]
  } else if (input$health_econ_perspective_checkbox == "Healthcare perspective") {
    tt <-
      out_proc()[, .(total_hcp_cost_cml,
        mc,
        friendly_name = factor(friendly_name),
        year)][, .(prop_if(total_hcp_cost_cml <= 0)), by = .(friendly_name, year)][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]
  }

  plot_ly(
    tt,
    x = ~ year,
    y = ~ V2,
    type = "scatter",
    mode = "lines+markers",
    color = ~ friendly_name,
    colors = colsymb()$colour,
    symbol = ~ friendly_name,
    symbols = colsymb()$symbol,
    line = list(shape = "spline", smoothing = 1.3)
  ) %>%
    add_lines(
      x = ~ year,
      y = input$decision_aid_gui,
      name = "Decision aid",
      color = NULL,
      symbol = NULL,
      line = list(color = "black", dash = "dot")
    ) %>%
    layout(
      yaxis = list(
        title = "Probability of cost-effective policy",
        range = c(-0.05, 1.05),
        tickformat = ",.0%"
      ),
      xaxis = list(title = "Year")
    )
})

# EQU plane ----
output$equ1_1 <- renderPlotly({
  tt <-
    out_proc_qimd()[year == max(year), .(
      net_utility_cml,
      nmb_cml,
      eq5d_cml,
      pops_cml,
      mc,
      friendly_name = factor(friendly_name),
      qimd
    )]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  tt <-
    tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name")]

  max_x <- tt[, max(abs(sei))] * 1.2
  max_y <- tt[, max(abs(nmb_cml))] * 1.2

  if (input$res_display_equ1_1)
    tt <- median_dt(tt, "friendly_name", "mc", digits = 1)

  p <-
    plot_ly(
      tt,
      x = ~ sei,
      y = ~ nmb_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      # frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Net monetary benefit (£)"),
      xaxis = list(title = "Absolute inequality reduction (SII)"),
      shapes = list(
        list(
          type = "rect",
          fillcolor = "green",
          line = list(color = "green"),
          opacity = 0.3,
          layer = "below",
          x0 = 0,
          x1 = max_x,
          xref = "x",
          y0 = 0,
          y1 = max_y,
          yref = "y"
        ),
        list(
          type = "rect",
          fillcolor = "red",
          line = list(color = "red"),
          opacity = 0.3,
          layer = "below",
          x0 = 0,
          x1 = -max_x,
          xref = "x",
          y0 = 0,
          y1 = -max_y,
          yref = "y"
        )
      )
    )

  # p <- animation_opts(p, frame = 1000, redraw = FALSE)
  # p <- animation_slider(p,
  #                       currentvalue = list(prefix = "Year: ",
  #                                           font = list(color = "red")))
})

output$equ1 <- renderPlotly({
  tt <- out_proc_qimd()[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml,
    pops_cml, mc, friendly_name = factor(friendly_name), qimd)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name")]

  max_x <- tt[, max(abs(sei))] * 1.2
  max_y <- tt[, max(abs(nmb_cml))] * 1.2

  if (input$res_display_equ1) tt <- median_dt(tt, "friendly_name", "mc", digits = 1)

  p <-
    plot_ly(
      tt,
      x = ~ sei,
      y = ~ nmb_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      # frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Net monetary benefit (£)"),
      xaxis = list(title = "Absolute inequality reduction (SII)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = max_y, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y")
      )
    )


  # p <- animation_opts(p, frame = 1000, redraw = FALSE)
  # p <- animation_slider(p,
  #                       currentvalue = list(prefix = "Year: ",
  #                                           font = list(color = "red")))
})

output$equ_rel <- renderPlotly({
  tt <- out_proc_qimd()[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml,
    pops_cml, mc, friendly_name = factor(friendly_name), qimd)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_rei(tt, c("mc", "friendly_name"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), rei = mean(rei)), by = c("mc", "friendly_name")]

  max_x <- tt[, max(abs(rei))] * 1.2
  max_y <- tt[, max(abs(nmb_cml))] * 1.2

  if (input$res_display_equ_rel) tt <- median_dt(tt, "friendly_name", "mc", digits = 6)


  p <-
    plot_ly(
      tt,
      x = ~ rei,
      y = ~ nmb_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      # frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Net monetary benefit (£)"),
      xaxis = list(title = "Relative inequality reduction (RII)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = max_y, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y")
      )
    )


  # p <- animation_opts(p, frame = 1000, redraw = FALSE)
  # p <- animation_slider(p,
  #                       currentvalue = list(prefix = "Year: ",
  #                                           font = list(color = "red")))
})

output$equ_anim_abs <- renderPlotly({
  tt <- out_proc_qimd()[, .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc,
    friendly_name = factor(friendly_name), qimd, year)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd", "year"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name", "year"))
  calc_sei(tt, c("mc", "friendly_name", "year"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name", "year")]


  max_x <- tt[, max(abs(sei))] * 1.2
  max_y <- tt[, max(abs(nmb_cml))] * 1.2

  p <-
    plot_ly(
      tt,
      x = ~ sei,
      y = ~ nmb_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Net monetary benefit (£)"),
      xaxis = list(title = "Absolute inequality reduction (SII)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = max_y, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y")
      ))


  p <- animation_opts(p, frame = 1000, redraw = FALSE)
  p <- animation_slider(p,
                        currentvalue = list(prefix = "Year: ",
                                            font = list(color = "red")))
})

output$equ_anim_rel <- renderPlotly({
  tt <- out_proc_qimd()[, .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc,
    friendly_name = factor(friendly_name), qimd, year)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd", "year"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name", "year"))
  calc_rei(tt, c("mc", "friendly_name", "year"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), rei = mean(rei)), by = c("mc", "friendly_name", "year")]

  max_x <- tt[, max(abs(rei))] * 1.2
  max_y <- tt[, max(abs(nmb_cml))] * 1.2

  p <-
    plot_ly(
      tt,
      x = ~ rei,
      y = ~ nmb_cml,
      color = ~ friendly_name,
      colors = colsymb()$colour,
      frame = ~ year,
      type = "scatter",
      mode = "markers",
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      showlegend = TRUE
    )
  p <-
    layout(
      p,
      yaxis = list(title = "Net monetary benefit (£)"),
      xaxis = list(title = "Relative inequality reduction (SII)"),
      shapes = list(
        list(type = "rect",
             fillcolor = "green", line = list(color = "green"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = max_x, xref = "x",
             y0 = 0, y1 = max_y, yref = "y"),
        list(type = "rect",
             fillcolor = "red", line = list(color = "red"), opacity = 0.3,
             layer = "below",
             x0 = 0, x1 = -max_x, xref = "x",
             y0 = 0, y1 = -max_y, yref = "y")
      ))


  p <- animation_opts(p, frame = 1000, redraw = FALSE)
  p <- animation_slider(p,
                        currentvalue = list(prefix = "Year: ",
                                            font = list(color = "red")))
})

output$equ_p_abs <- renderPlotly({
  tt <- out_proc_qimd()[, .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc,
    friendly_name = factor(friendly_name), qimd, year)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd", "year"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name", "year"))
  calc_sei(tt, c("mc", "friendly_name", "year"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei)), by = c("mc", "friendly_name", "year")
  ][, (prop_if(sei > 0)), by = .(friendly_name, year)
  ][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]


  plot_ly(
    tt,
    x = ~ year,
    y = ~ V2,
    type = "scatter",
    mode = "lines+markers",
    color = ~ friendly_name,
    colors = colsymb()$colour,
    symbol = ~ friendly_name,
    symbols = colsymb()$symbol,
    line = list(shape = "spline", smoothing = 1.3)
  ) %>%
    add_lines(
      x = ~ year,
      y = input$decision_aid_gui,
      name = "Decision aid",
      color = NULL,
      symbol = NULL,
      line = list(color = "black", dash = "dot")
    ) %>%
    layout(
      yaxis = list(
        title = "Probability of equitable policy",
        range = c(-0.05, 1.05),
        tickformat = ",.0%"
      ),
      xaxis = list(title = "Year")
    )

})

output$equ_p_rel <- renderPlotly({
  tt <- out_proc_qimd()[, .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc,
    friendly_name = factor(friendly_name), qimd, year)
  ]
  # [, sum_dt(.SD, c("mc", "friendly_name", "qimd", "year"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name", "year"))
  calc_sei(tt, c("mc", "friendly_name", "year"))
  calc_rei(tt, c("mc", "friendly_name", "year"))
  tt <- tt[, .(nmb_cml = sum(nmb_cml), sei = mean(sei), rei = mean(rei)), by = c("mc", "friendly_name", "year")
  ][, (prop_if(rei > 0 & sei > 0)), by = .(friendly_name, year)
  ][, V2 := clamp(predict(loess(V1 ~ year, span = 0.5))), by = friendly_name]

  plot_ly(
    tt,
    x = ~ year,
    y = ~ V2,
    type = "scatter",
    mode = "lines+markers",
    color = ~ friendly_name,
    colors = colsymb()$colour,
    symbol = ~ friendly_name,
    symbols = colsymb()$symbol,
    line = list(shape = "spline", smoothing = 1.3)
  ) %>%
    add_lines(
      x = ~ year,
      y = input$decision_aid_gui,
      name = "Decision aid",
      color = NULL,
      symbol = NULL,
      line = list(color = "black", dash = "dot")
    ) %>%
    layout(
      yaxis = list(
        title = "Probability of equitable policy",
        range = c(-0.05, 1.05),
        tickformat = ",.0%"
      ),
      xaxis = list(title = "Year")
    )
})


# EFFECTIVENESS ----
output$cypp_1 <- renderPlotly({
  tt <- out_proc()[year == max(year), .(cypp_chd_cml,
                                        cypp_stroke_cml,
                                        cypp_poststroke_dementia_cml,
                                        cypp_t2dm_cml,
                                        # cypp_af_cml,
                                        # cypp_htn_cml,
                                        cypp_copd_cml,
                                        cypp_lung_ca_cml,
                                        cypp_colon_ca_cml,
                                        cypp_breast_ca_cml,
                                        mc,
                                        friendly_name = factor(friendly_name))
  ][, median_dt(.SD, c("friendly_name"), "mc")]

  # [, sum_dt(.SD, c("mc", "friendly_name"), character(0))]

  p <-
    plot_ly(
      tt,
      x = ~ friendly_name,
      y = ~ cypp_chd_cml,
      name = '',
      text = 'CHD',
      textposition = 'auto',
      insidetextfont = list(
        size = 15,
        color = 'black',
        opacity = 1
      ),
      type = 'bar',
      marker = list(
        opacity = 0.6,
        line = list(color = colsymb()$colour, width = 5)
      )
    ) %>%
    add_trace(
      y = ~ cypp_stroke_cml,
      name = "",
      text = 'stroke',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cypp_poststroke_dementia_cml,
      name = "",
      text = 'poststroke dementia',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    # add_trace(y = ~cypp_af_cml, name = "", text = 'af', textposition = 'auto', insidetextfont = list(size=15, color = 'black')) %>%
    # add_trace(y = ~cypp_htn_cml, name = "", text = 'HTN', textposition = 'auto', insidetextfont = list(size=15, color = 'black')) %>%
    add_trace(
      y = ~ cypp_t2dm_cml,
      name = "",
      text = 'T2DM',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cypp_copd_cml,
      name = "",
      text = 'COPD',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cypp_lung_ca_cml,
      name = "",
      text = 'lung cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cypp_colon_ca_cml,
      name = "",
      text = 'colon cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cypp_breast_ca_cml,
      name = "",
      text = 'breast cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%

    layout(
      yaxis = list(title = 'Cases'),
      showlegend = FALSE,
      barmode = 'relative'
    )

})

output$cpp_1 <- renderPlotly({
  tt <- out_proc()[year == max(year), .(cpp_chd_cml,
                                        cpp_stroke_cml,
                                        cpp_poststroke_dementia_cml,
                                        cpp_copd_cml,
                                        cpp_t2dm_cml,
                                        # cypp_af_cml,
                                        # cpp_htn_cml,
                                        cpp_lung_ca_cml,
                                        cpp_colon_ca_cml,
                                        cpp_breast_ca_cml,
                                        mc,
                                        friendly_name = factor(friendly_name))
  ][, median_dt(.SD, c("friendly_name"), "mc")]

  # [, sum_dt(.SD, c("mc", "friendly_name"), character(0))]
  p <-
    plot_ly(
      tt,
      x = ~ friendly_name,
      y = ~ cpp_chd_cml,
      name = '',
      text = 'CHD',
      textposition = 'auto',
      insidetextfont = list(
        size = 15,
        color = 'black',
        opacity = 1
      ),
      type = 'bar',
      marker = list(
        opacity = 0.6,
        line = list(color = colsymb()$colour, width = 5)
      )
    ) %>%
    add_trace(
      y = ~ cpp_stroke_cml,
      name = "",
      text = 'stroke',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cpp_poststroke_dementia_cml,
      name = "",
      text = 'post-stroke dementia',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    # add_trace(y = ~cpp_af_cml, name = "", text = 'af', textposition = 'auto', insidetextfont = list(size=15, color = 'black')) %>%
    # add_trace(y = ~cpp_htn_cml, name = "", text = 'HTN', textposition = 'auto', insidetextfont = list(size=15, color = 'black')) %>%
    add_trace(
      y = ~ cpp_t2dm_cml,
      name = "",
      text = 'T2DM',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cpp_copd_cml,
      name = "",
      text = 'COPD',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cpp_lung_ca_cml,
      name = "",
      text = 'lung cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cpp_colon_ca_cml,
      name = "",
      text = 'colon cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%
    add_trace(
      y = ~ cpp_breast_ca_cml,
      name = "",
      text = 'breast cancer',
      textposition = 'auto',
      insidetextfont = list(size = 15, color = 'black')
    ) %>%

    layout(
      yaxis = list(title = 'Cases'),
      showlegend = FALSE,
      barmode = 'relative'
    )

})

output$dpp_1 <- renderPlotly({
  tt <-
    out_proc()[year == max(year), .(dpp_all_cause_cml, mc,
      friendly_name = factor(friendly_name))
    ][, median_dt(.SD, c("friendly_name"), "mc")]
  # [, sum_dt(.SD, c("mc", "friendly_name"), character(0))]
  plot_ly(
    tt,
    x = ~ friendly_name,
    y = ~ dpp_all_cause_cml,
    color = ~ friendly_name,
    colors = colsymb()$colour,
    text = 'DPP',
    textposition = 'auto',
    insidetextfont = list(
      size = 15,
      color = 'black',
      opacity = 1
    ),
    type = 'bar',
    marker = list(opacity = 0.7)
  ) %>%
    layout(
      yaxis = list(title = 'All-cause deaths'),
      showlegend = FALSE,
      barmode = 'relative'
    )

})

output$dppy_spline <- renderPlotly({
  tt <- out_proc()[, .(dpp_all_cause_cml,
                       mc,
                       friendly_name = factor(friendly_name), year)
  ][, median_dt(.SD, c("friendly_name", "year"), "mc")
  ][, V2 := round(predict(loess(dpp_all_cause_cml ~ year))), by = .(friendly_name)]
  # [, sum_dt(.SD, c("mc", "friendly_name", "year"), character(0))]

  plot_ly(
    tt,
    x = ~ year,
    y = ~ V2,
    type = "scatter",
    mode = "lines+markers",
    color = ~ friendly_name,
    colors = colsymb()$colour,
    symbol = ~ friendly_name,
    symbols = colsymb()$symbol,
    line = list(shape = "spline", smoothing = 1.3)
  ) %>%

    layout(
      yaxis = list(title = "Cumulative number of deaths<br>prevented or postponed",
                   title_standoff = 15),

      xaxis = list(title = "Years")
    )
})

output$cyppy_spline <- renderPlotly({
  if (is.null(input$out_diseases_select_cyppy)) {
    tt_n <- as.data.table(NULL)
    plot_ly(
      tt_n,
      x = 0,
      y = 0,
      type = "scatter",
      mode = "lines+markers",
      line = list(shape = "spline", smoothing = 1.3)
    )
  } else {
    tt <- out_proc()[, .(
      cypp_chd_cml,
      cypp_stroke_cml,
      cypp_poststroke_dementia_cml,
      cypp_t2dm_cml,
      # cypp_af_cml,
      # cypp_htn_cml,
      cypp_copd_cml,
      cypp_lung_ca_cml,
      cypp_colon_ca_cml,
      cypp_breast_ca_cml,
      mc,
      year,
      friendly_name = factor(friendly_name)
    )][, median_dt(.SD, c("friendly_name", "year"), "mc")
    ][, cpp_all := Reduce("+", mget(input$out_diseases_select_cyppy))
    ][, V2 := round(predict(loess(cpp_all ~ year, span = 0.5))),
      by = .(friendly_name)]
    # [, sum_dt(.SD, c("mc", "friendly_name", "year"), character(0))]
    plot_ly(
      tt,
      x = ~ year,
      y = ~ V2,
      type = "scatter",
      mode = "lines+markers",
      color = ~ friendly_name,
      colors = colsymb()$colour,
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      line = list(shape = "spline", smoothing = 1.3)
    ) %>%

      layout(
        yaxis = list(title = "Cumulative number of case-years<br>prevented or postponed",
                     title_standoff = 15),
        xaxis = list(title = "Years")
      )
  }


})

output$cppy_spline <- renderPlotly({
  if (is.null(input$out_diseases_select_cppy)) {
    tt_n <- as.data.table(NULL)
    plot_ly(
      tt_n,
      x = 0,
      y = 0,
      type = "scatter",
      mode = "lines+markers",
      line = list(shape = "spline", smoothing = 1.3)
    )
  } else {
    tt <- out_proc()[, .(
      cpp_chd_cml,
      cpp_stroke_cml,
      cpp_poststroke_dementia_cml,
      cpp_t2dm_cml,
      # cpp_af_cml,
      # cpp_htn_cml,
      cpp_copd_cml,
      cpp_lung_ca_cml,
      cpp_colon_ca_cml,
      cpp_breast_ca_cml,
      mc,
      year,
      friendly_name = factor(friendly_name)
    )][, median_dt(.SD, c("friendly_name", "year"), "mc")
    ][, cpp_all := Reduce("+", mget(input$out_diseases_select_cppy))
    ][, V2 := round(predict(loess(cpp_all ~ year, span = 0.5))),
      by = .(friendly_name)]
    # [, sum_dt(.SD, c("mc", "friendly_name", "year"), character(0))]


    plot_ly(
      tt,
      x = ~ year,
      y = ~ V2,
      type = "scatter",
      mode = "lines+markers",
      color = ~ friendly_name,
      colors = colsymb()$colour,
      symbol = ~ friendly_name,
      symbols = colsymb()$symbol,
      line = list(shape = "spline", smoothing = 1.3)
    ) %>%

      layout(
        yaxis = list(title = "Cumulative number of cases<br>prevented or postponed",
                     title_standoff = 15),
        xaxis = list(title = "Years")
      )
  }


})

scn_nams <- reactive(out_proc()[, as.character(unique(friendly_name))])

output$out_scenario_disease_select_cypp <- renderUI({
  tagList(
    pickerInput(inputId = "inout_scenario_diseases_select_cypp",
                label = "Scenario",
                choices = scn_nams(),
                selected = first(scn_nams()),
                options = list(`actions-box` = TRUE, `live-search` = FALSE),
                multiple = FALSE)    %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Please select the scenario that you would like to plot.")
      )
  )
})

output$out_scenario_disease_select_cpp <- renderUI({
  tagList(
    pickerInput(inputId = "inout_scenario_diseases_select_cpp",
                label = "Scenario",
                choices = scn_nams(),
                selected = first(scn_nams()),
                options = list(`actions-box` = TRUE, `live-search` = FALSE),
                multiple = FALSE)    %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Please select the scenario that you would like to plot.")
      )
  )
})



output$cyppy_stacked_area <- renderPlotly({
  scn_sel <-
    reactive({
      ifelse(
        length(input$inout_scenario_diseases_select_cypp) == 0,
        first(scn_nams()),
        input$inout_scenario_diseases_select_cypp
      )
    })
  tt <- out_proc()[, .(
    cypp_chd_cml,
    cypp_stroke_cml,
    cypp_poststroke_dementia_cml,
    cypp_t2dm_cml,
    # cypp_af_cml,
    # cypp_htn_cml,
    cypp_copd_cml,
    cypp_lung_ca_cml,
    cypp_colon_ca_cml,
    cypp_breast_ca_cml,
    mc,
    year,
    friendly_name = factor(friendly_name)
  )][, median_dt(.SD, c("friendly_name", "year"), "mc")
  ][friendly_name == scn_sel(),
    c("CHD", "Stroke", "Post-stroke dementia",
      "T2DM", "COPD", "Lung Cancer",
      "Colon cancer", "Breast cancer") := .(
        round(predict(loess(cypp_chd_cml ~ year))),
        round(predict(loess(cypp_stroke_cml ~ year))),
        round(predict(loess(cypp_poststroke_dementia_cml ~ year))),
        round(predict(loess(cypp_t2dm_cml ~ year))),
        # predict(loess(cypp_af_cml ~ year)),
        # predict(loess(cypp_htn_cml ~ year)),
        round(predict(loess(cypp_copd_cml ~ year))),
        round(predict(loess(cypp_lung_ca_cml ~ year))),
        round(predict(loess(cypp_colon_ca_cml ~ year))),
        round(predict(loess(cypp_breast_ca_cml ~ year)))
      ),
    by = .(friendly_name)][friendly_name == scn_sel(),]
  tt[, (grep("_cml$", names(tt))) := NULL]
  tt <- melt(tt, id.vars = 1:2)[, id := as.integer(factor(variable))]
  tt <- split(tt, by = "variable")
  p <-
    lapply(tt, function(d) {
      plot_ly(
        d,
        x = ~ year,
        y = ~ value,
        mode = "lines+markers",
        color = ~ variable,
        colors = viridis_pal(option = "C")(unique(d$id)) #,
        # colors = colsymb()$colour,
        # symbol = ~ friendly_name, symbols = colsymb()$symbol
      ) %>%
        add_lines() %>% add_annotations(
          text = unique(d$variable),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        ) %>% layout(
          xaxis = list(showgrid = FALSE),
          yaxis = list(rangemode = "tozero",
                       showgrid = FALSE)
        )
    })

  subplot(
    p,
    nrows = 4,
    shareX = TRUE,
    titleX = FALSE,
    titleY = FALSE,
    margin = 0.05
  ) %>%
    layout(
      title = scn_sel(),
      legend = list(
        orientation = 'h',
        xanchor = "center",
        x = 0.5,
        y = -0.18
      )
    ) %>%
    add_annotations(
      text = "Years",
      x = 0.5,
      y = 0,
      yref = "paper",
      xref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      yshift = -45,
      showarrow = FALSE,
      font = list(size = 15)
    ) %>%
    add_annotations(
      text = "Case-years",
      x = 0,
      y = 0.5,
      yref = "paper",
      xref = "paper",
      xanchor = "left",
      yanchor = "center",
      textangle = -90,
      xshift = -55,
      showarrow = FALSE,
      font = list(size = 15)
    )
})


output$cppy_stacked_area <- renderPlotly({
  scn_sel <-
    reactive({
      ifelse(
        length(input$inout_scenario_diseases_select_cpp) == 0,
        first(scn_nams()),
        input$inout_scenario_diseases_select_cpp
      )
    })
  tt <- out_proc()[, .(
    cpp_chd_cml,
    cpp_stroke_cml,
    cpp_poststroke_dementia_cml,
    cpp_t2dm_cml,
    # cpp_af_cml,
    # cpp_htn_cml,
    cpp_copd_cml,
    cpp_lung_ca_cml,
    cpp_colon_ca_cml,
    cpp_breast_ca_cml,
    mc,
    year,
    friendly_name = factor(friendly_name)
  )][, median_dt(.SD, c("friendly_name", "year"), "mc")
  ][friendly_name == scn_sel(),
    c("CHD", "Stroke", "Post-stroke dementia",
      "T2DM", "COPD", "Lung Cancer",
      "Colon cancer", "Breast cancer") := .(
        round(predict(loess(cpp_chd_cml ~ year))),
        round(predict(loess(cpp_stroke_cml ~ year))),
        round(predict(loess(cpp_poststroke_dementia_cml ~ year))),
        round(predict(loess(cpp_t2dm_cml ~ year))),
        # round(predict(loess(cpp_af_cml ~ year))),
        # round(predict(loess(cpp_htn_cml ~ year))),
        round(predict(loess(cpp_copd_cml ~ year))),
        round(predict(loess(cpp_lung_ca_cml ~ year))),
        round(predict(loess(cpp_colon_ca_cml ~ year))),
        round(predict(loess(cpp_breast_ca_cml ~ year)))
      ),
    by = .(friendly_name)][friendly_name == scn_sel(),]
  tt[, (grep("_cml$", names(tt))) := NULL]
  tt <- melt(tt, id.vars = 1:2)[, id := as.integer(factor(variable))]
  tt <- split(tt, by = "variable")
  p <-
    lapply(tt, function(d) {
      plot_ly(
        d,
        x = ~ year,
        y = ~ value,
        mode = "lines+markers",
        color = ~ variable,
        colors = viridis_pal(option = "C")(unique(d$id)) #,
        # colors = colsymb()$colour,
        # symbol = ~ friendly_name, symbols = colsymb()$symbol
      ) %>%
        add_lines() %>% add_annotations(
          text = unique(d$variable),
          xref = "paper",
          yref = "paper",
          yanchor = "bottom",
          xanchor = "center",
          align = "center",
          x = 0.5,
          y = 1,
          showarrow = FALSE
        ) %>% layout(
          xaxis = list(showgrid = FALSE),
          yaxis = list(rangemode = "tozero",
                       showgrid = FALSE)
        )
    })

  subplot(
    p,
    nrows = 4,
    shareX = TRUE,
    titleX = FALSE,
    titleY = FALSE,
    margin = 0.05
  ) %>%
    layout(
      title = scn_sel(),
      legend = list(
        orientation = 'h',
        xanchor = "center",
        x = 0.5,
        y = -0.18
      )
    ) %>%
    add_annotations(
      text = "Years",
      x = 0.5,
      y = 0,
      yref = "paper",
      xref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      yshift = -45,
      showarrow = FALSE,
      font = list(size = 15)
    ) %>%
    add_annotations(
      text = "Cases",
      x = 0,
      y = 0.5,
      yref = "paper",
      xref = "paper",
      xanchor = "left",
      yanchor = "center",
      textangle = -90,
      xshift = -55,
      showarrow = FALSE,
      font = list(size = 15)
    )
})


# myModal <- function() {
#   div(id = "test",
#       modalDialog(downloadButton("download","Download the tab as csv"),
#                   easyClose = TRUE, title = "Download Table")
#   )
# }


# output$tbl_raw <- DT::renderDT(
#       out_proc_raw()[1:10, ],
#       style = "default", # others mess download buttons
#       callback = JS("$('div.dwnld').append($('#download_raw'));"),
#       extensions = "Buttons",
#       options = list(scrollX = TRUE, scrollY = "300px",
#                      stateSave = FALSE,
#                      dom = 'B<"dwnld">frtip',
#                      buttons = c('copy', 'print')),
#       filter = "none",
#       caption = htmltools::tags$caption(
#         style = 'caption-side: top; text-align: left;',
#         'Table 2: ',
#         htmltools::em('Raw model outputs. For efficiency only a hundred rows are visible but you can download the full dataset as .csv'))
#     )

# output$out_columns_select_raw <- renderUI({
#   tagList(
#     pickerInput(inputId  = "inout_columns_select_raw",
#                 label    = "Select table columns",
#                 choices  = names(out_proc_raw()),
#                 selected = names(out_proc_raw()),
#                 options  = list(`actions-box` = TRUE, `live-search` = FALSE),
#                 multiple = TRUE)    %>%
#       shinyInput_label_embed(
#         icon("info") %>%
#           bs_embed_popover(title = "Please select the columns you want to see in the tab.")
#       )
#   )
# })

output$download_raw <- downloadHandler(
  filename = function() {
    paste("IMPACTncd_Engl_raw_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    fwrite(out_proc_raw(), file)
  }
)

observe({
  out_summary()

  # when it updates, save the search strings so they're not lost
  isolate({
    # update global search and column search strings
    default_search_tbl_summary <- input$dataTable_search
    default_search_columns_tbl_summary <-
      c("", input$dataTable_search_columns)

    # update the search terms on the proxy table (see below)
    proxy_tbl_summary %>% updateSearch(
      keywords =
        list(global = default_search_tbl_summary,
             columns = default_search_columns_tbl_summary)
    )
  })
})

# from https://dev.to/awwsmm/reactive-datatables-in-r-with-persistent-filters-l26
# and https://github.com/rstudio/DT/issues/267 for download all button
output$tbl_summary <- DT::renderDT(
  out_summary(),
  style = "default", # others mess up download buttons
  # callback = JS("$('div.dwnld').append($('#download_summary'));"),
  extensions = "Buttons",
  options = list(scrollX = TRUE, scrollY = "420 px",
                 stateSave = FALSE,
                 dom = 'B<"dwnld">frtip',
                 buttons = list(
                   'copy', 'print'
                 ),
                 # default column search strings and global search string
                 searchCols = default_search_columns_tbl_summary,
                 search = list(regex = FALSE, caseInsensitive = FALSE, search = default_search_tbl_summary)),
  filter = list(position = "top"),
  caption = htmltools::tags$caption(
    style = 'caption-side: top; text-align: left;',
    'Table 1: ',
    htmltools::em('Summarised model outputs and their uncertainty')))

# make a proxy of the data table so it can be edited after it's been rendered
proxy_tbl_summary <- dataTableProxy("tbl_summary")

output$download_summary <- downloadHandler(
  filename = function() {
    paste("IMPACTncd_Engl_summary_", Sys.Date(), ".csv", sep = "")
  },
  content = function(file) {
    fwrite(out_summary(), file)
  }
)



# output$out_columns_select_summary <- renderUI({
#   tagList(
#     pickerInput(inputId  = "inout_columns_select_summary",
#                 label    = "Select table columns",
#                 choices  = names(out_summary()),
#                 selected = names(out_summary()),
#                 options  = list(`actions-box` = TRUE, `live-search` = FALSE),
#                 multiple = TRUE)    %>%
#       shinyInput_label_embed(
#         icon("info") %>%
#           bs_embed_popover(title = "Please select the columns you want to see in the tab.")
#       )
#   )
# })


outputOptions(output, "out_scenario_select", suspendWhenHidden = FALSE)
outputOptions(output, "out_year_slider", suspendWhenHidden = FALSE)
# outputOptions(input, "out_discount_slider", suspendWhenHidden = FALSE)

