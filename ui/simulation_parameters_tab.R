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

tabPanel(
  "Simulation parameters",
  wellPanel(fluidRow(
    column(
      6,
      pickerInput(inputId = "locality_select",
                  label = "Area to simulate",
                  choices = localitities_list,
                  options = list(`actions-box` = TRUE, `live-search` = TRUE,
                    liveSearchNormalize = TRUE, header = "Close list"),
                  multiple = TRUE)
      %>%
        shinyInput_label_embed(
          shiny_iconlink("info") %>%
            bs_embed_popover(title = "Please select the locality for your simulation",
                             content = "The model will use a synthetic population that closely resembles that of the selected area. If you select multiple areas, the overlapping ones will only be included once. For example, if you select both 'England' and 'London', London will be ignored because it is already included in England. The first time you select an area combination, workHORSE will generate and save the synthetic population for that area combination which may take up to 30 min. Subsequent selections of the same area combination will be much faster as the model will reuse the already existing synthetic population.

      The 'Produce synthpops for area selection' button in the advanced settings tab initiates the production of the synthetic populations in the background while you design the scenarios, to save time. WorkHORSE will proceed as normal even if you do not press this button, although it will take longer to produce results the first time you simulate an area.
                             ",
                             placement = "bottom")
        ) #,

      # actionButton(
      #   "locality_select_validator",
      #   "Confirm area selection",
      #   icon = icon("check-circle"),
      #   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
      #   #class = "btn-info",
      #   width = "auto"
      # )
      ),
    column(
      6,
      sliderInput(
        "simulation_period_slider",
        "Period to simulate",
        2013,
        2041,
        c(2021, 2041),
        1,
        sep = "",
        ticks = FALSE
      ) %>%
        shinyInput_label_embed(
           shiny_iconlink("info") %>%
            bs_embed_popover(title = "Please select the start and end year of the simulation")
    ))
  )),

  wellPanel(fluidRow(
    column(
      6,
      sliderInput(
        "scenarios_number_slider",
        "Number of scenarios to be simulated",
        1,
        9,
        value = 2,
        1,
        sep = "",
        ticks = FALSE
      ) %>%
        shinyInput_label_embed(
          # shiny_iconlink("info") %>%
           shiny_iconlink("info") %>%
            bs_embed_popover(title = "Please choose how many scenarios you want to simulate")
        )
    ),
    column(
      5,
      radioButtons(
        "health_econ_perspective_checkbox",
        "Health economic analysis perspective",
        choices = c("Societal perspective",
          "Health and social care perspective",
          "Healthcare perspective"
        ),
        selected = c("Societal perspective"))
    ),
    column( # This only needed for the icon of the tooltip
      1, align = 'right',
      # shiny_iconlink("info") %>%
       shiny_iconlink("info") %>%
        bs_embed_popover(title = "Please choose the perspective of your cost:benefit analysis")
    )
  )),

  wellPanel(fluidRow(column(
    6,
    switchInput(
      "national_qimd_checkbox",
      # "Produce results by local IMD",
      value = TRUE,
      onLabel = "National IMD",
      offLabel = "Local IMD",
      offStatus = "danger",
      width = "100%"
    )
    %>%
      shinyInput_label_embed(
         shiny_iconlink("info") %>%
          bs_embed_popover(title = "Please choose to use either the local or the national quantile groups of the Index of Multideprivation (QIMD). Your choice here affects every aspect of workHORSE that uses QIMD including scenario inputs and model results.")
      ))
  # ,
  # column(
  #   6,
  #   switchInput(
  #     "ward_output_checkbox",
  #     "Report results by Ward",
  #     value = FALSE,
  #     onLabel = "Yes",
  #     offLabel = "No",
  #     labelWidth = "100%"
  #   )
  #   %>%
  #     shinyInput_label_embed(
  #        shiny_iconlink("info") %>%
  #         bs_embed_popover(title = "Please choose if you would like the results to be presented at the Ward level")
  #     ))
  )),
  icon = icon("user-cog")
)
