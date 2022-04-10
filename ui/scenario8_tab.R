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

conditionalPanel(
  condition = "input.level == 8",

  bsCollapse(
    id = "collapse_panel_sc8",
    multiple = TRUE,
    open = c(
      "General Parameters"
      # "Eligibility Criteria",
      # "Appointments Offered Yearly (%)",
      # "Health Checks received (%)",
      # "Prescription Rate",
      # "Impact on Lifestyle"
    ),


# General panel ----------------------------------------------------------

bsCollapsePanel(
  "General Parameters",
  style = "default",
  wellPanel(fluidRow(
    ## col 1 ----
    column(
      2, align = "left", style='padding:0px;',
      div(
      textInput("friendly_name_sc8", "Friendly name", "Scenario 8")
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please pick a name for this scenario.
                             This will be used in all model outputs.")
        )),
      div(colourpicker::colourInput(
        inputId = "col_sc8",
        label = "Pick a color for the scenario",
        value = def_col_small[[8]],
        showColour = "background",
        palette = "limited",
        allowedCols = def_col,
        returnName = FALSE
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please pick a colour for this scenario. This will be displayed in all model outputs.")
        )),
      uiOutput("panel_col_sc8"),

      div(selectInput(
        inputId = "symbol_sc8",
        label = "Pick a symbol for the scenario",
        choices = def_sym,
        selected = def_sym[[8]]
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please pick a symbol for this scenario. This will be displayed in all model outputs to help for differentiation.")
        ))

    ),
    ## col 2 ----
    column(
      2, offset = 2, align = "left", style='padding:0px;',
      div(
      switchInput(
        "baseline_sc8",
        "Baseline scenario",
        FALSE,
        onLabel = "Yes",
        offLabel = "No",
        labelWidth = "100%"
        # size = "large",
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
             bs_embed_popover(title = "Please check if you would like this to be
                             the baseline scenario which will serve as point of
                             comparison with all other scenarios. In every run,
                             you have to have one and only one baseline scenario
                             (or ensemble scenario).",
                       placement = "top")
        ),
      div(uiOutput("init_year_slider_sc8"))

)),
    ## col 3 ----
    column(
      3, offset = 2, align = "left", style='padding:0px;',
      div(actionButton(
        "collapse_panels_button_sc8",
        "Expand/Collapse panels",
        icon = icon("bars")
      ), style='vertical-align:top;'),
      br(),
      br(),
      div(downloadButton(
        "save_sc8",
        "Save scenario         ",
        icon = icon("save"),
        # class = "btn-info",
        labelWidth = "100%"
      ), style='vertical-align:center;'),
      br(),
      div(fileInput(
        "load_sc8",
        "",
        multiple = FALSE,
        accept = ".yaml",
        placeholder = "",
        buttonLabel = "Load scenario"
      ), style='vertical-align:bottom;') # ,
      # tags$style("
      #        .btn-file {
      #        background-color:red;
      #        border-color: red;
      #        }
      #
      #        .progress-bar {
      #        background-color: red;
      #        }
      #
      #        ")
    ),
    ## col 4 ----
    column(
      1, align = "left",
      div(icon("info") %>%
        bs_embed_popover(title = "Expands all scenario parameter panels, if all
                         panels are collapsed. Collapses all scenario parameter
                         panels, if at least one is expanded.")),
      br(),
      br(),
      div(icon("info") %>%
        bs_embed_popover(title = "Saves the current scenario in a .rds file
                                  in your default 'download' folder for reuse or
                                  archiving purposes")),
      br(),
      br(),
      div(icon("info") %>%
        bs_embed_popover(title = "Loads a .rds file
                                  which contains a previously saved scenario
                                  specification. WARNING: This will overwrite
                                  all the parameters currently in this scenario."))
      )
  )
)
  ),


# Eligibility panel ----------------------------------------------------------
bsCollapsePanel("Eligibility Criteria",
                style = "success",
                wellPanel(fluidRow(
                  ## col 1 ----
                  column(
                    4,
                    sliderInput(
                      "age_eligibility_slider_sc8",
                      "Age eligibility",
                      30,
                      84,
                      c(40, 74),
                      1,
                      sep = "",
                      ticks = FALSE
                    )
                      %>%
                        shinyInput_label_embed(
                          icon("info") %>%
                            bs_embed_popover(title = "Please select the age range of people to be invited")
                        )
                    ,
                    sliderInput(
                      "frequency_eligibility_slider_sc8",
                      "Minimum allowed years between two Health Checks",
                      1,
                      25,
                      5,
                      1,
                      sep = "",
                      ticks = FALSE
                    )
                    %>%
                      shinyInput_label_embed(
                        icon("info") %>%
                          bs_embed_popover(title = "Please select the minimum number of years between two concecutive Health Checks")
                      )
                  ),
                  ## col 2 ----
                  column(
                    4,
                    switchInput(
                      "invite_known_hypertensives_checkbox_sc8",
                      "Invite people known to have  hypertension",
                      value = FALSE,
                      onLabel = "Yes",
                      offLabel = "No",
                      labelWidth = "100%"
                    )
                    %>%
                      shinyInput_label_embed(
                        icon("info") %>%
                          bs_embed_popover(title = "Please select whether known hypertensives should be eligible for a Health Check")
                    ),
                    br(),
                    br(),
                    switchInput(
                      "invite_known_diabetics_checkbox_sc8",
                      "Invite people known to have diabetes",
                      value = FALSE,
                      onLabel = "Yes",
                      offLabel = "No",
                      labelWidth = "100%"
                    )
                    %>%
                      shinyInput_label_embed(
                        icon("info") %>%
                          bs_embed_popover(title ="Please select whether known diabetics should be eligible for a Health Check")
                      )
                  ),
                  ## col 3 ----
                  column(
                    4,
                    switchInput(
                      "cancel_program_checkbox_sc8",
                      "Nobody eligible",
                      value = FALSE,
                      onLabel = "Yes",
                      offLabel = "No",
                      labelWidth = "100%"
                    )
                    %>%
                      shinyInput_label_embed(
                        icon("info") %>%
                          bs_embed_popover(title = "Please select to simulate a no Health Check scenario. It hides all relevant parameters from the user interface")
                      )
                  )
                ))),

# Coverage panel ----------------------------------------------------------
bsCollapsePanel(
  "Appointments Offered Yearly (%)",
  style = "success",
  wellPanel(
    fluidRow(
      ## col 1 ----
      column(
        4,
        sliderInput(
          "coverage_qimd0_slider_sc8",
          "Invitations (percentage of eligible population)",
          0,
          100,
          20,
          0.5,
          sep = "",
          ticks = FALSE,
          post  = " %"
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please select the percentage of eligible population that is invited every year")
          )
      ),
      ## col 2 ----
      column(
        4,
        numericInput(
          "coverage_cost_qimd0_sc8",
          "Cost per invitation",
          20,
          0,
          500,
          1
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please enter the cost per invitation")
          )
      ),
      ## col 3 ----
      column(
        4,
        br(),
        switchInput(
          "coverage_detailed_checkbox_sc8",
          "Detailed input",
          value = FALSE,
          onLabel = "Show",
          offLabel = "Hide",
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Reveals more granular inputs")
          )
      )
    ),
    ## row 2 ----
    fluidRow(
      ## col 1 ----
      column(
      4,
      sliderInput(
        "coverage_qimd1_slider_sc8",
        "Invitations as percentage of eligible population in QIMD 1 (most deprived)",
        0,
        100,
        20,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please select the percentage of eligible population in QIMD 1 (most deprived) areas that is invited every year")
        )
    ),
    ## col 2 ----
    column(
      4,
      numericInput(
        "coverage_cost_qimd1_sc8",
        "Cost per invitation in QIMD 1 (most deprived)",
        20,
        0,
        500,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per invitation in QIMD 1 (most deprived) areas")
        )
    )),
    ## row 3 ----
    fluidRow(
      ## col 1 ----
      column(
      4,
      sliderInput(
        "coverage_qimd2_slider_sc8",
        "Invitations as percentage of eligible population in QIMD 2",
        0,
        100,
        20,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please select the percentage of eligible population in QIMD 2 areas that is invited every year")
        )
    ),
    ## col 2 ----
    column(
      4,
      numericInput(
        "coverage_cost_qimd2_sc8",
        "Cost per invitation in QIMD 2",
        20,
        0,
        500,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per invitation in QIMD 2 areas")
        )
    )),
    ## row 4 ----
    fluidRow(
      ## col 1 ----
      column(
      4,
      sliderInput(
        "coverage_qimd3_slider_sc8",
        "Invitations as percentage of eligible population in QIMD 3",
        0,
        100,
        20,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please select the percentage of eligible population in QIMD 3 areas that is invited every year")
        )
    ),
    ## col 2 ----
    column(
      4,
      numericInput(
        "coverage_cost_qimd3_sc8",
        "Cost per invitation in QIMD 3",
        20,
        0,
        500,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per invitation in QIMD 3 areas")
        )
    )),
    fluidRow(column(
      4,
      sliderInput(
        "coverage_qimd4_slider_sc8",
        "Invitations as percentage of eligible population in QIMD 4",
        0,
        100,
        20,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please select the percentage of eligible population in QIMD 4 areas that is invited every year")
        )
    ),
    column(
      4,
      numericInput(
        "coverage_cost_qimd4_sc8",
        "Cost per invitation in QIMD 4",
        20,
        0,
        500,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per invitation in QIMD 4 areas")
        )
    )),
    fluidRow(column(
      4,
      sliderInput(
        "coverage_qimd5_slider_sc8",
        "Invitations as percentage of eligible population in QIMD 5 (least deprived)",
        0,
        100,
        20,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please select the percentage of eligible population in QIMD 5 (least deprived) areas that is invited every year")
        )
    ),
    column(
      4,
      numericInput(
        "coverage_cost_qimd5_sc8",
        "Cost per invitation in QIMD 5 (least deprived)",
        20,
        0,
        500,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per invitation in QIMD 5 (least deprived) areas")
        )
    ))
  )
),

# Uptake panel ----------------------------------------------------------

bsCollapsePanel("Health Checks Received (%)",
                style = "success",
                wellPanel(
                  fluidRow(
                    ## col 1 ----
                    column(
                      4,
                      sliderInput(
                        "uptake_slider_sc8",
                        "Proportion of invitees attending a Health Check",
                        0,
                        100,
                        66,
                        0.5,
                        sep = "",
                        ticks = FALSE,
                        post  = " %"
                      )
                      %>%
                        shinyInput_label_embed(
                          icon("info") %>%
                            bs_embed_popover(title = "Please select the percentage of people who complete a Health Check after an invitation")
                        )
                    ),
                    ## col 2 ----
                    column(
                      4,
                      numericInput("uptake_cost_sc8",
                                   "Cost per completed Health Check", 20, 0, 500, 1)
                      %>%
                        shinyInput_label_embed(
                          icon("info") %>%
                            bs_embed_popover(title = "Please enter the cost per individual completing a Health Check")
                        )
                    ),
                    ## col 3 ----
                    column(
                      4,
                      br(),
                      switchInput(
                        "uptake_detailed_checkbox_sc8",
                        "Detailed input",
                        value = FALSE,
                        onLabel = "Show",
                        offLabel = "Hide",
                        labelWidth = "100%"
                      )
                      %>%
                        shinyInput_label_embed(
                          icon("info") %>%
                            bs_embed_popover(title = "Reveal more granular inputs. It is highly recommended to inform the model with this additional information.")
                        )
                    )
                  ),
                  fluidRow(
                    br(),
                    uiOutput("uptake_table_help_sc8"),
                    tableOutput("uptake_table_sc8"),
                    switchInput(
                      "uptake_structural_zeroes_sc8",
                      "Assume structural 0s",
                      value = FALSE,
                      onLabel = "Yes",
                      offLabel = "No",
                      labelWidth = "100%"
                    )
                    %>%
                      shinyInput_label_embed(
                        icon("info") %>%
                          bs_embed_popover(title = "If yes, 0s in the cells at the table above mean that the respective groups are excluded from participating in a Health Check. This option may be useful in parallel ensemble scenarios. The default is 'No', which means that groups with 0 in their respective cells, still have a small probability to complete a Health Check. This probability increases as the proportion of invitees attending a Health Check is increasing.")
                      )
                  )
                )),


# Prescription rate panel  ------------------------------------------------
bsCollapsePanel(
  "Prescription Rate",
  style = "success",
  wellPanel(
    fluidRow(
      ## col 1 ----
      column(
        4,
        sliderInput(
          "statin_px_slider_sc8",
          "Proportion of participants that have been prescribed a statin",
          0,
          100,
          10,
          0.5,
          sep = "",
          ticks = FALSE,
          post  = " %"
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please select the percentage of attendees who have completed a Health Check and have been prescribed a statin as a consequence. Note that the denomiator here are all attendees, not just those that a statin would be recommended.")
          )
      ),
      ## col 2 ----
      column(
        4,
        sliderInput(
          "antihtn_px_slider_sc8",
          "Proportion of participants that have been prescribed antihypertensives",
          0,
          100,
          10,
          0.5,
          sep = "",
          ticks = FALSE,
          post  = " %"
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please select the percentage of attendees who have completed a Health Check and have been prescribed antihypertensive medication as a consequence. Note that the denomiator here are all attendees, not just those that treatment for hypertension would be recommended." )
          )
      ),
      ## col 3 ----
      column(
        4,
        br(),
        switchInput(
          "px_detailed_checkbox_sc8",
          "Detailed input",
          value = FALSE,
          onLabel = "Show",
          offLabel = "Hide",
          labelWidth = "100%"
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Reveal more granular options. It is highly recommended to inform the model with this additional information.")
          )
      )
    )),

  wellPanel(
    fluidRow(
    uiOutput("statin_px_table_help_sc8"),
    tableOutput("statin_px_table_sc8")
)),

wellPanel(
  fluidRow(
    uiOutput("antihtn_px_table_help_sc8"),
    tableOutput("antihtn_px_table_sc8")
  ))
),

# Lifestyle panel  ------------------------------------------------
bsCollapsePanel(
  "Impact on Lifestyle",
  style = "success",
  wellPanel(helpText(h5('Smoking cessation')),
            fluidRow(
              ## col 1 ----
              column(
                4,
                sliderInput(
                  "smkcess_slider_sc8",
                  "Percentage of smoker participants successfully quit smoking for at least a year",
                  0,
                  100,
                  0,
                  0.5,
                  sep = "",
                  ticks = FALSE,
                  post  = " %"
                )
                  %>%
                    shinyInput_label_embed(
                      icon("info") %>%
                        bs_embed_popover(title = "Please enter the percentage of people quitting smoking in the last year")
                    )
              ),
              ## col 2 ----
              column(
                4,
                numericInput(
                  "smkcess_cost_sc8",
                  "Smoking cessation cost per successful quit",
                  0,
                  0,
                  5e6,
                  1
                )
                %>%
                  shinyInput_label_embed(
                    icon("info") %>%
                      bs_embed_popover(title = "Please enter the cost per participant in the smoking cessation program")
                  )
              ),
              ## col 3 ----
              column(
                4,
                numericInput(
                  "smkcess_cost_ovrhd_sc8",
                  "Smoking cessation annual overhead costs",
                  0,
                 -5e6,
                  5e6,
                  1
                )
                %>%
                  shinyInput_label_embed(
                    icon("info") %>%
                      bs_embed_popover(title = "Please enter the annual fixed cost of the smoking cessation program")
                  )
              )
            )),
  wellPanel(
  helpText(h5('Weight management')),
  fluidRow(
    column(
      4,
      sliderInput(
        "wghtpct_slider_sc8",
        "Percentage of obese participants losing weight",
        0,
        100,
        0,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )  %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the percentage of obese participants losing weigth")
      ),
      sliderInput(
        "wghtreduc_slider_sc8",
        "Average weight loss as percentage weight", # equivalent to percentage of BMI
        0,
        20,
        0,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter average percentage reduction in BMI per person in a year, e.g. 5% reduction in BMI")
        )
    ),
    column(
      4,
      numericInput(
        "wghtloss_cost_sc8",
        "Weight management annual cost per participant losing weight",
        0,
        0,
        5e6,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per participant in the weight management program")
        )
    ),
    column(
      4,
      numericInput(
        "wghtloss_cost_ovrhd_sc8",
        "Weight management annual overhead costs",
        0,
       -5e6,
        5e6,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the annual fixed cost of the weight management program")
        )
    )
  )
  ),
  wellPanel(
  helpText(h5('Physical activity')),
  fluidRow(
    column(
      4,
      sliderInput(
        "papct_slider_sc8",
        "Percentage of participants increasing their physical activity",
        0,
        100,
        0,
        0.5,
        sep = "",
        ticks = FALSE,
        post  = " %"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the percentage of participants increasing their physical activity levels")
      ),
      sliderInput(
        "papincr_slider_sc8",
        "Average increase in active days per week",
        0,
        7,
        0,
        1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the average increase in days of physical activity, e.g. 2 days")
        )
    ),
    column(
      4,
      numericInput(
        "pa_cost_sc8",
        "Physical activity annual cost per more active participant",
        0,
        0,
        5e6,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the cost per participant in the the physical activity program")
        )
    ),
    column(
      4,
      numericInput(
        "pa_cost_ovrhd_sc8",
        "Physical activity actions annual overhead costs",
        0,
       -5e6,
        5e6,
        1
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Please enter the annual fixed cost of the physical activity program")
        )
    )
  )),
  wellPanel(
    helpText(h5('Alcohol')),
    fluidRow(
      column(
        4,
        sliderInput(
          "alcoholpct_slider_sc8",
          "Percentage of participants consuming more than 14 units of alcohol per week, cutting down",
          # 1 unit is 8g of alcohol
          0,
          100,
          0,
          0.5,
          sep = "",
          ticks = FALSE,
          post  = " %"
        )  %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please enter the percentage of participants consuming more than 14 units of alcohol per week that cut down")
          ),
        sliderInput(
          "alcoholreduc_slider_sc8",
          "Alcohol intake percentage reduction",
          0,
          50,
          0,
          0.5,
          sep = "",
          ticks = FALSE,
          post  = " %"
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please enter average percentage reduction in alcohol intake, e.g. 5% reduction in units of alcohol consumed per week")
          )
      ),
      column(
        4,
        numericInput(
          "alcoholreduc_cost_sc8",
          "Alcohol service annual cost per participant losing weight",
          0,
          0,
          5e6,
          1
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please enter the cost per participant in the alcohol reduction service")
          )
      ),
      column(
        4,
        numericInput(
          "alcoholreduc_cost_ovrhd_sc8",
          "Alcohol service annual overhead costs",
          0,
          -5e6,
          5e6,
          1
        )
        %>%
          shinyInput_label_embed(
            icon("info") %>%
              bs_embed_popover(title = "Please enter the annual fixed cost of the alcohol reduction service")
          )
      )
    )
  ),
  wellPanel(helpText(h5('Attrition rate')),
            fluidRow(
              ## col 1 ----
              column(
                4,
                sliderInput(
                  "lifestyle_attrition_slider_sc8",
                  "Percentage of participants reverting back to pre health check lifestyle every year",
                  0,
                  100,
                  20,
                  0.5,
                  sep = "",
                  ticks = FALSE,
                  post  = " %"
                )
                %>%
                  shinyInput_label_embed(
                    icon("info") %>%
                      bs_embed_popover(title = "Please enter the percentage of people that revert to their previous lifestyle every year, after succesfully improving their lifestyle because of a health check. For example, 20% attrition rate per year means that after 5 years only (1-0.2)^5 = 33% of those successfuly improved their lifestyle initially are still observe the lifestyle change. The same attrition rate applies to all of the risk factor above, except for smoking.")
                  )
              )
            )
          )
),

# Advanced panel  ------------------------------------------------
bsCollapsePanel(
  "Advanced",
  style = "default",
  wellPanel(fluidRow(
    ## Scenario ensembles ----
    h5("Scenario ensembles"),
    h6("In this section you can specify the above scenario to be part of a serial
       or parallel ensemble. Scenarios in the same ensemble need to share the same
       friendly name!"),
    HTML("<br> <br/> "),
    column(6,
  switchInput(
    "serial_ensemble_checkbox_sc8",
    "This scenario is part of a serial ensemble",
    value = FALSE,
    onLabel = "Yes",
    offLabel = "No",
    labelWidth = "100%"
  )
  %>%
    shinyInput_label_embed(
      icon("info") %>%
        bs_embed_popover(title = "This allows you to run simultaneously different scenarios for different time spans, e.g. scenario 1 from years 1-5 and scenario 2 from years 6-10. ")
    )
  )),
  fluidRow(div(style="display: inline-block;vertical-align:top; width: 10px; height: 30px",HTML("<br>"))),
fluidRow(
    column(6,
           switchInput(
             "parallel_ensemble_checkbox_sc8",
             "This scenario is part of a parallel ensemble",
             value = FALSE,
             onLabel = "Yes",
             offLabel = "No",
             labelWidth = "100%"
           )
           %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "This allows you to run simultaneously different scenarios for different groups of the population, e.g lifestyle interventions for the most deprived and more appointments offered for the individuals at risk")
             )
    ),
    column(6,
    sliderInput(
      "parallel_ensemble_slider_sc8",
      "Percentage of population this scenario applies to",
      0,
      100,
      0,
      0.5,
      sep = "",
      ticks = FALSE,
      post  = " %"
    )
    %>%
      shinyInput_label_embed(
        icon("info") %>%
          bs_embed_popover(title = "Percentage of the population this scenario applies to if it is part of a parallel or sequential ensemble")
      )
    ))),
  wellPanel(fluidRow(
    ## Digital Health Checks ----
    h5("Digital Health Checks"),
    HTML("<br> <br/> "),
    column(3,
           switchInput(
             "ignore_cholesterol_checkbox_sc8",
             "Ignore cholesterol in QRISK calculation",
             value = FALSE,
             onLabel = "Yes",
             offLabel = "No",
             labelWidth = "50%",
             width = "100%"
           )
           %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "When calculating the Qrisk function this option allows you to exclude cholesterol estimates")
             )
    ),
    column(4, offset = 1,
           switchInput(
             "ignore_sbp_checkbox_sc8",
             "Ignore systolic blood pressure in QRISK calculation",
             value = FALSE,
             onLabel = "Yes",
             offLabel = "No",
             labelWidth = "100%",
             width = "100%"
           )
           %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "When calculating the Qrisk function this option allows you to exclude systolic blood pressure estimates")
             )
    ),
    column(3, offset = 1,
           switchInput(
             "ignore_bmi_checkbox_sc8",
             "Ignore BMI in QRISK calculation",
             value = FALSE,
             onLabel = "Yes",
             offLabel = "No",
             labelWidth = "50%",
             width = "100%"
           )
           %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "When calculating the Qrisk function this option allows you to exclude BMI estimates")
             )
    )
    )),
  wellPanel(fluidRow(
    ## Structural Policies ----
    h5("Structural Policies / Calibration"),
    h6("In this section the effect of implementing structural/population level
       strategies  can be modelled, in addition",
       "to the health checks programm specified above",
       "On occasion, this section can be used to calibrate the synthetic
        population exposure to risk factors."),
    HTML("<br> <br/> "),
    column(6,
           sliderInput(
             "structural_smk_slider_sc8",
             "Smoking prevalence relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )
           %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative (%) in smoking
                                  prevalence as result of a structural policy.
                                  It does not affect the prevalence of never
                                  smokers. When negative, it forces more smokers
                                  to quit (i.e. it increases the probability of
                                  cessation). When positive, smokers are coming
                                  from the pool of ex-smokers (i.e. it increases
                                  the probability of relapse).")
           ),
           # sliderInput(
           #   "structural_neversmk_slider_sc8",
           #   "Never-smoking prevalence absolute percentage change across
           #   population",
           #   -20,
           #   20,
           #   0,
           #   0.1,
           #   sep = "",
           #   ticks = FALSE,
           #   post  = " %"
           # )
           # %>%
           #   shinyInput_label_embed(
           #     icon("info") %>%
           #       bs_embed_popover(title = "Absolute change (%) in the
           # cprevalence of smoking as result of a structural policy")
           #   ),
           sliderInput(
             "structural_fv_slider_sc8",
             "Fruit & veg relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative change (%) in the prevalence
                                  of fruit and vegetable consumption as result
                                  of a structural policy")
           ),
           sliderInput(
             "structural_alcohol_slider_sc8",
             "Alcohol relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative change (%) in the population
                                  mean of alcohol consumption as result of a
                                  structural policy")
           ),
           sliderInput(
             "structural_bmi_slider_sc8",
             "BMI relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative change (%) in the population mean of BMI as result of a structural policy")
           ),
           sliderInput(
             "structural_sbp_slider_sc8",
             "Systolic blood pressure relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative change (%) in the population mean of systolic blood pressure as result of a structural policy")
           ),
           sliderInput(
             "structural_chol_slider_sc8",
             "Cholesterol relative percentage change across population",
             -20,
             20,
             0,
             0.1,
             sep = "",
             ticks = FALSE,
             post  = " %"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Relative change (%) in the population mean of cholesterol as result of a structural policy")
           ),
           sliderInput(
             "structural_pa_slider_sc8",
             "Physical activity change in active days per week across
             population",
             -7,
             7,
             0,
             1,
             sep = "",
             ticks = TRUE,
             post  = " active days/week"
           )  %>%
             shinyInput_label_embed(
               icon("info") %>%
                 bs_embed_popover(title = "Change in the active days per week
                 as result of a structural policy")
             )
    ),
    column(3

    ),
    column(3

    )
  )),
  wellPanel(fluidRow(
    ## Social Policies ----
    h5("Social Policies"),
    h6("Note: these policies may override some of the changes set above"),
    HTML("<br> <br/> "),
    column(6,
      sliderInput(
        "qimd_one_slider_sc8",
        "Reassign QIMD 1 (most deprived) to QIMD...",
        1, 5, 1, 1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Assign QIMD 1 to a different one")
        ),
      sliderInput(
        "qimd_two_slider_sc8",
        "Reassign QIMD 2 to QIMD...",
        1, 5, 2, 1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Assign QIMD 2 to a different one")
        ),
      sliderInput(
        "qimd_three_slider_sc8",
        "Reassign QIMD 3 to QIMD...",
        1, 5, 3, 1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Assign QIMD 3 to a different one")
        ),
      sliderInput(
        "qimd_four_slider_sc8",
        "Reassign QIMD 4 to QIMD...",
        1, 5, 4, 1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Assign QIMD 4 to a different one")
        ),
      sliderInput(
        "qimd_five_slider_sc8",
        "Reassign QIMD 5 (least deprived) to QIMD...",
        1, 5, 5, 1,
        sep = "",
        ticks = FALSE,
        post  = ""
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Assign QIMD 5 to a different one")
        )

    ),
    column(6,
      switchInput(
        "qimd_disease_fatalities_checkbox_sc8",
        "Should QIMD reassignment affect diseases case fatalities?",
        value = FALSE,
        onLabel = "Yes",
        offLabel = "No",
        labelWidth = "100%"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Should QIMD reassignment affect disease case fatalities?")
        ),
      switchInput(
        "qimd_nonmodelled_fatality_checkbox_sc8",
        "Should QIMD reassignment affect nonmodelled case fatality?",
        value = FALSE,
        onLabel = "Yes",
        offLabel = "No",
        labelWidth = "100%"
      )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Should QIMD reassignment affect nonmodelled diseases case fatalitys?")
        ),
      checkboxGroupButtons("qimd_risk_factors_checkbox_sc8",
        "Which risk factors should be affected by the changes in QIMD?",
        c("Smoking" = "smok",
          "Environmental Tobacco Smoking" = "ets",
          "Fruit & Veg" = "fv",
          "Alcohol" = "alc",
          "Physical Activity" = "pa",
          "BMI" = "bmi",
          "SBP" = "sbp",
          "Cholesterol" = "tchol"),
        selected = c("smok",
          "ets",
          "fv",
          "alc",
          "pa",
          "bmi" ,
          "sbp",
          "tchol"),
        direction = "vertical",
        status = "default",
        checkIcon = list(
          yes = icon("check-square"),
          no = icon("square")
        ),
        justified = TRUE,
        individual = FALSE
        )
      %>%
        shinyInput_label_embed(
          icon("info") %>%
            bs_embed_popover(title = "Which risk factors should be affected by the changes in QIMD?")
        )
    ) # end column
  ))
),

# Notes  ------------------------------------------------
bsCollapsePanel("Notes",
                style = "default",
                wellPanel(fluidRow(
                  textAreaInput(
                    "notes_sc8",
                    "",
                    "",
                    width = "100%",
                    rows = 6,
                    resize = "both"
                  )
                  %>%
                    shinyInput_label_embed(
                      icon("info") %>%
                        bs_embed_popover(title = "Add a description of this scenario or any important notes to remember")
                    )
                )))
  ),

conditionalPanel(condition = "input.level <= input.scenarios_number_slider",
                 fluidRow(
                 column(6,
                        actionButton(
                          "previous_sc8",
                          "Go to previous scenario",
                          icon = icon("step-backward"),
                          #style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                          class = "btn btn-primary",
                          width = "100%"
                        )),
                 column(
                   6,
                   actionButton(
                     "next_sc8",
                     "Go to next scenario",
                     icon = icon("step-forward"),
                     #style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                     class = "btn btn-primary",
                     width = "100%"
                   )
                 ))),

conditionalPanel(condition = "input.level == input.scenarios_number_slider",
                 fluidRow(
                   br(),
                   bsAlert("alert_locality_sc8"),
                   br(),
                   bsAlert("alert_baseline_sc8"),
                   br(),
                   actionButton(
                     "run_simulation_sc8",
                     "Run simulation (all scenarios)",
                     icon = icon("paper-plane"),
                     #style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
                     class = "btn-info",
                     width = "100%"
                   )
                 )))



