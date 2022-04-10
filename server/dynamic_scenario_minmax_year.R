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

output$init_year_slider_sc1 <- renderUI({
  tagList(
    sliderInput(
      "ininit_year_slider_sc1",
      "First year of implementation",
      input$simulation_period_slider[[1]], input$simulation_period_slider[[2]],
      input$simulation_period_slider[[1]], 1,
      sep = "",
      ticks = FALSE
    ) %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Please select the first year of implementation for this scenario.
                             If this scenario is part of serial scenarion ensemble,
                             then the last year of implementation for this scenario
                             will up to the first year of implementation for the next scenario in the ensemble.")
      )
  )
})
