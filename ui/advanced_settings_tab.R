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

tabPanel("Advanced settings",
  # Sysload ----

  wellPanel(
    h5("System load"),
    p("CPU and RAM utilisation"),
    verbatimTextOutput("sysload", placeholder = TRUE)
  ),


         # Theme selector ----
         wellPanel(
         h5("Theme selector"),
         p("Changes the colour theme of the interface. The default is 'paper'"),
         mythemeSelector()
         ),

         # sim parameters ----
         wellPanel(
           h5("Simulation parameters"),
           p("Please do not alter unless you know what you are doing!"),

           numericInput(
             "iteration_n_gui",
             "Number of Monte Carlo iterations for the interactive exploration (GUI)",
             design$sim_prm$iteration_n,
             1L,
             1e3,
             1L
           ),

           numericInput(
             "iteration_n_final_gui",
             "Number of Monte Carlo iterations for the final results (GUI)",
             design$sim_prm$iteration_n_final,
             1L,
             1e3,
             1L
           ),

           numericInput(
             "clusternumber_gui",
             "Number of cores to be used for explicit parallelisation",
             design$sim_prm$clusternumber,
             1L,
             1e3,
             1L
           ),

           numericInput(
             "n_cpus_gui",
             "Number of cores to be used for implicit parallelisation",
             design$sim_prm$n_cpus,
             1L,
             1e3,
             1L
           ),

           numericInput(
             "n_gui",
             "Size of synthetic population",
             design$sim_prm$n,
             1e5,
             1e6,
             1e5L
           ),

           numericInput(
             "n_synthpop_aggregation_gui",
             "Number of synthetic population files to combine together",
             design$sim_prm$n_synthpop_aggregation,
             1L,
             10L,
             1L
           ),

           numericInput(
             "n_primers_gui",
             "Number of synthetic population primer files to be produced. Better be a multiple of the setting above.",
             design$sim_prm$n_primers,
             1L,
             100L,
             1L
           ),

           numericInput(
             "cvd_lag_gui",
             "Median lag between exposure and CVD incidence",
             design$sim_prm$cvd_lag,
             2L,
             9L,
             1L
           ),

           numericInput(
             "copd_lag_gui",
             "Median lag between exposure and COPD incidence",
             design$sim_prm$copd_lag,
             2L,
             9L,
             1L
           ),

           numericInput(
             "cancer_lag_gui",
             "Median lag between exposure and cancer incidence",
             design$sim_prm$cancer_lag,
             2L,
             9L,
             1L
           ),

           numericInput(
             "nonmodelled_lag_gui",
             "Median lag between exposure and death from causes not explicitly modelled",
             design$sim_prm$nonmodelled_lag,
             2L,
             9L,
             1L
           ),

           numericInput(
             "cancer_cure_gui",
             "Number of years after which a prevalent cancer case is considered cured",
             design$sim_prm$cancer_cure,
             2L,
             10L,
             1L
           ),

           numericInput(
             "jumpiness_gui",
             "Increase for more erratic jumps in trajectories",
             design$sim_prm$jumpiness,
             0.1,
             10,
             0.1
           ),

           numericInput(
             "statin_adherence_gui",
             "Statin adherence. The mean of a beta distribution with shape2 = 0.2",
             design$sim_prm$statin_adherence,
             0.1,
             1,
             0.01
           ),

           numericInput(
             "bpmed_adherence_gui",
             "BP medication adherence. The mean of a beta distribution with shape2 = 0.2",
             design$sim_prm$bpmed_adherence,
             0.1,
             1,
             0.05
           ),

           numericInput(
             "decision_aid_gui",
             "Adjust the decision aid line used in some of the graphs",
             design$sim_prm$decision_aid,
             0,
             1,
             0.05
           )

         ), # wellPanel

  # Logging ----
  wellPanel(
    h5("Enable or disable logging"),
    p("Enable logging to save R warning and errors in text files"),

    switchInput(
      "logs_gui",
      "Enable logging",
      design$sim_prm$logs,
      onLabel = "Yes",
      offLabel = "No",
      labelWidth = "100%"
      # size = "large",
    ),

    downloadBttn(
      "download_logs_gui",
      "Download logs",
      "simple",
      color = "primary",
      size = "sm",
      block = TRUE
      # icon = icon("cloud-download-alt")
    ),

    br(),

    tags$head(tags$script(src = "message-handler.js")),
    actionButton(
      "delete_logs_gui",
      "Delete logs",
      icon = icon("trash-alt"),
      class = "btn-danger",
      width = "100%"
    )

  ), # wellPanel

  # Synthpop files ----
  wellPanel(
    h5("Manage synthpop files"),
    p("Actions to manipulate synthpop files"),

    actionButton(
      "locality_select_validator",
      "Produce synthpops for area selection",
      icon = icon("check-circle"),
      style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
      #class = "btn-info",
      width = "100%" # "auto"
    ),

    verbatimTextOutput("checksum", placeholder = TRUE),

    verbatimTextOutput("synthpop_info", placeholder = TRUE),

    actionButton(
      "check_synthpops_gui",
      "Check the integridy of the synthpop files and delete any orphan or malformed files",
      icon = icon("trash-alt"),
      class = "btn-warning",
      width = "100%"
    ),
    # . Make sure that you are the only user running the model when you click
    # this. Otherwise it may crash the session of other concurrent users

    actionButton(
      "check_synthpops_gui_check_checksum",
      "Check the integridy of the synthpop files, which produced from parameters same as in the GUI, and delete any orphan or malformed files",
      icon = icon("trash-alt"),
      class = "btn-warning",
      width = "100%"
    ),
    # Make sure that you are the only user running the model when you click
    # this. Otherwise it may crash the session of other concurrent users

    actionButton(
      "delete_synthpops_gui_spare_primer",
      "Delete synthpop files which produced from parameters same as in the GUI,  sparing the primer files",
      icon = icon("trash-alt"),
      class = "btn-danger",
      width = "100%"
    ),

    actionButton(
      "delete_synthpops_gui_check_checksum",
      "Delete synthpop files which produced from parameters same as in the GUI, including the primer files",
      icon = icon("trash-alt"),
      class = "btn-danger",
      width = "100%"
    ),

    actionButton(
      "delete_synthpops_gui_all",
      "Delete everything in the synthpop folder!",
      icon = icon("trash-alt"),
      class = "btn-danger",
      width = "100%"
    )

  ) # wellPanel
)
