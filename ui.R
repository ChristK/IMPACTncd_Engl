## workHORSE is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## workHORSE is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

# myenv <- environment()
ui = tagList(
  use_bs_tooltip(),
  use_bs_popover(),
  useShinyjs(), # for conditional deactivating inputs
  chooseSliderSkin("Flat", color = "#337ab7"),
  # 'Nice' is also good

  navbarPage(
    theme = shinytheme("paper"),
    #cosmo

    title = "workHORSE",

    id = "inTabset",
    # Welcome tab --------------------------
    source(file.path("ui", "welcome_tab.R"), local = TRUE)$value,

    # Simulation parameters tab --------------------------
    source(file.path("ui", "simulation_parameters_tab.R"), local = TRUE)$value,

    # Scenario parameters tab --------------------------
    source(file.path("ui", "scenario_parameters_tab.R"), local = TRUE)$value,

    # Output tab --------------------------
    source(file.path("ui", "output_tab.R"), local = TRUE)$value,
    # br(), # causes a warning from bslib

    navbarMenu("More",
      # Advance settings tab --------------------------
      source(
        file.path("ui", "advanced_settings_tab.R"), local = TRUE
      )$value,

      # Documentation tab --------------------------
      source(
        file.path("ui", "documentation_tab.R"), local = TRUE
      )$value,

      # About tab --------------------------
      source(file.path("ui", "about_tab.R"), local = TRUE)$value
      )
  ),


  # Add version number to the right
  HTML(
    "<script>var parent = document.getElementsByClassName('navbar-nav');
parent[0].insertAdjacentHTML( 'afterend', '<ul class=\"nav navbar-nav navbar-right\"><li class=\"disabled\"><a href=\"#\">v0.99b</a></li></ul>' );</script>"
  )
)

