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

observeEvent(input$age_eligibility_slider_sc1, {
  # Dynamic age groups ----
  if (input$age_eligibility_slider_sc1[[2]] - input$age_eligibility_slider_sc1[[1]] < 11) {
    agegrp_nam <- paste0(input$age_eligibility_slider_sc1[[1]], "-", input$age_eligibility_slider_sc1[[2]])
  } else {
  agegrp_nam <- agegrp_name(floor(input$age_eligibility_slider_sc1[[1]]/10) * 10L,
                            floor(input$age_eligibility_slider_sc1[[2]]/10) * 10L, 10L)
  }
  l_agegrp_nam <- length(agegrp_nam)

  # Prepare uptake table ----------------------------------------------------
  output$uptake_table_sc1 <- renderTable({
    num.inputs.col1 <- paste0("<input id='uptake_qrisk_low_",  seq_len(l_agegrp_nam * 10), "_sc1", "' class='shiny-bound-input' type='number' value=`0.007`>")
    num.inputs.col2 <- paste0("<input id='uptake_qrisk_mid_",  seq_len(l_agegrp_nam * 10), "_sc1", "' class='shiny-bound-input' type='number' value=`0.007`>")
    num.inputs.col3 <- paste0("<input id='uptake_qrisk_high_", seq_len(l_agegrp_nam * 10), "_sc1", "' class='shiny-bound-input' type='number' value=`0.007`>")
    CJ(QIMD = c("1 (most deprived)", "2", "3", "4", "5 (least deprived)"),
       Sex = c("men", "women"),
       Ages = agegrp_nam)[,
                          c("Low risk (QRISK <10)", "Mid risk (QRISK 10-20)", "High risk (QRISK 20+)") := .(num.inputs.col1, num.inputs.col2, num.inputs.col3)
                          ]
  }, spacing = "xs",
  sanitize.text.function = function(x) x)
})

output$uptake_table_help_sc1 <- renderUI({
  helpText(
    p(
      'In the following table, please enter the number of participants over a period of time by sex, age, QIMD and QRISK score.'
    )
  )
})



# Prepare statins Px table -------------------------------
  output$statin_px_table_sc1 <- renderTable({
    num.inputs.col1 <- paste0("<input id='statin_px_qrisk_low_",  1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    num.inputs.col2 <- paste0("<input id='statin_px_qrisk_mid_",  1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    num.inputs.col3 <- paste0("<input id='statin_px_qrisk_high_", 1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    setnames(data.frame(QIMD =
                          c("1 (most deprived)", 2, 3, 4, "5 (least deprived)"),
                        num.inputs.col1, num.inputs.col2, num.inputs.col3),
             c("num.inputs.col1", "num.inputs.col2", "num.inputs.col3"),
             c("Low risk (QRISK <10)", "Mid risk (QRISK 10-20)", "High risk (QRISK 20+)")
    )
  }, spacing = "xs",
  sanitize.text.function = function(x) x)

output$statin_px_table_help_sc1 <- renderUI({
  helpText(
    p(
      'In the following table, please enter the number of participants that were prescribed a statin after a Health Check by QIMD, and QRISK score.'
    )
  )
})

# Prepare antihtn Px table -------------------------------
  output$antihtn_px_table_sc1 <- renderTable({
    num.inputs.col1 <- paste0("<input id='antihtn_px_qrisk_low_",  1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    num.inputs.col2 <- paste0("<input id='antihtn_px_qrisk_mid_",  1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    num.inputs.col3 <- paste0("<input id='antihtn_px_qrisk_high_", 1:5, "_sc1", "' class='shiny-bound-input' type='number' value=`0`>")
    setnames(data.frame(QIMD =
                          c("1 (most deprived)", 2, 3, 4, "5 (least deprived)"),
                        num.inputs.col1, num.inputs.col2, num.inputs.col3),
             c("num.inputs.col1", "num.inputs.col2", "num.inputs.col3"),
             c("Low risk (QRISK <10)", "Mid risk (QRISK 10-20)", "High risk (QRISK 20+)")
    )
  }, spacing = "xs",
  sanitize.text.function = function(x) x)

output$antihtn_px_table_help_sc1 <- renderUI({
  helpText(
    p(
      'In the following table, please enter the number of participants that were prescribed antihypertensive medication after a Health Check by QIMD, and QRISK score.'
    )
  )
})
