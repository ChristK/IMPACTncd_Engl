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

# From https://stackoverflow.com/questions/32922190/saving-state-of-shiny-app-to-be-restored-later
output$save_sc1 <- downloadHandler(


  filename = function()
    paste0(input$friendly_name_sc1, ".yaml"),

  content = function(file) {
    updateCollapse(session, "collapse_panel_sc1",
    open = c("Eligibility Criteria",
             "Appointments Offered Yearly (%)",
             "Health Checks Received (%)",
             "Prescription Rate",
             "Impact on Lifestyle",
             "Advanced",
             "Notes"))
    Sys.sleep(1)
   tt_sc1 <- lapply(reactiveValuesToList(input), unclass)

    # tt_sc1 <- reactiveValuesToList(input)
    tt_sc1 <- tt_sc1[grep("^.*_sc1$", names(tt_sc1))]
    tt_sc1$save_sc1 <- NULL
      tt_sc1$load_sc1 <- NULL
      tt_sc1$collapse_panels_button_sc1 <- NULL
      tt_sc1$collapse_panel_sc1 <- NULL
      tt_sc1$run_simulation_sc1 <- NULL
      tt_sc1$next_sc1 <- NULL
      tt_sc1$previous_sc1 <- NULL
      tt_sc1$baseline_sc1 <- NULL #TODO crashes if two loaded scenarios have baseline = TRUE

    tt_sc1 <- setNames(tt_sc1, gsub("_sc1$", "", names(tt_sc1)))
    write_yaml(tt_sc1[sort(names(tt_sc1))], file = file)
  }
)


  observeEvent(input$load_sc1, {
    updateCollapse(session, "collapse_panel_sc1",
                   open = c("Eligibility Criteria",
                            "Appointments Offered Yearly (%)",
                            "Health Checks Received (%)",
                            "Prescription Rate",
                            "Impact on Lifestyle",
                            "Advanced",
                            "Notes"))

    inFile <- input$load_sc1

    if (is.null(inFile)) return(NULL)
    savedInputs   <- read_yaml(inFile$datapath)
    names(savedInputs)      <- paste0(names(savedInputs), "_sc1")

    if (savedInputs[["uptake_detailed_checkbox_sc1"]]) showElement(id = "uptake_table_sc1")
 if (savedInputs[["px_detailed_checkbox_sc1"]]) {
   showElement(id = "statin_px_table_sc1")
   showElement(id = "antihtn_px_table_sc1")
 }

# delay necessary so all elements are expanded on screen before load
    delay(2000, lapply(sort(names(savedInputs)),
           function(x) session$sendInputMessage(x, list(value = savedInputs[[x]]))
    ))
# TODO When age eligibility changes between scenarios the loading of the
# detailed uptake table does not work the first time. As a workaround until
# better, fix I load each scenario twice
    lapply(sort(names(savedInputs)),
      function(x) session$sendInputMessage(x, list(value = savedInputs[[x]]))
    )
  })

# Ensure uptake & px is open and close once, otherwise the tables are not produced and they cannot be saved or load
observeEvent(input$collapse_panels_button_sc1, ({
  if (length(input$collapse_panel_sc1) > 1) {
    updateCollapse(session, "collapse_panel_sc1",
                   close = setdiff(input$collapse_panel_sc1, "General Parameters"))
  } else {
    updateCollapse(session, "collapse_panel_sc1",
                   open = c("General Parameters",
                            "Eligibility Criteria",
                            "Appointments Offered Yearly (%)",
                            "Health Checks Received (%)",
                            "Prescription Rate",
                            "Impact on Lifestyle"#,
                           # "Advanced",
                           # "Notes"
                            ))
  }
}))
