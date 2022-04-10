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
  title = "Output",
  value = "output_panel",
  dashboardPage(
    dashboardHeader(disable = TRUE),
    # Misbehaves when framed in a tabPanel
    dashboardSidebar(
      tags$head(tags$style(
        HTML(
          "right-side,
                        .sidebar-toggle {height: 100vh; overflow: hidden;}"
        )
      )),
      sidebarMenu(
        id = "out_tabs",
        menuItem(
          "Dashboard",
          tabName = "dashboard",
          icon = icon("tachometer-alt")
        ),
        menuItem(
          "Cost-effectiveness",
          icon = icon("pound-sign"),
          tabName = "cost_effectiveness"
        ),
        menuItem(
          "Health Inequalities",
          icon = icon("balance-scale"),
          tabName = "equity"
        ),
        menuItem(
          "Effectiveness",
          icon = icon("heart"),
          tabName = "effectiveness"
        ),
        # menuItem("Wider Social Benefits", icon = icon("users"), tabName = "wider_social_benef"
        # ),
        menuItem(
          "Tables",
          icon = icon("table"),
          menuSubItem("Summarised results", tabName = "summarised_results"),
          menuSubItem("Raw model output", tabName = "raw_output")

        ),
        menuItem(
          "Filters",
          uiOutput("out_year_slider"),
          uiOutput("out_scenario_select"),
          pickerInput(
            inputId = "out_characteristics_select",
            label = "Sub-population",
            choices = list(
              "Age groups" = c("30-49", "50-69", "70-89"),
              "Sex"        = c("men", "women"),
              "QIMD"       = c("1 most deprived", "2", "3", "4", "5 least deprived"),
              "Ethnicity"  = c(
                "white",
                "indian",
                "pakistani",
                "bangladeshi",
                "other asian",
                "black caribbean",
                "black african",
                "chinese",
                "other"
              )
            ),
            selected = c(
              "30-49",
              "50-69",
              "70-89",
              "men",
              "women",
              "1 most deprived",
              "2",
              "3",
              "4",
              "5 least deprived",
              "white",
              "indian",
              "pakistani",
              "bangladeshi",
              "other asian",
              "black caribbean",
              "black african",
              "chinese",
              "other"
            ),
            options = list(`actions-box` = TRUE, `live-search` = TRUE),
            multiple = TRUE
          )  %>%
            shinyInput_label_embed(
              icon("info") %>%
                bs_embed_popover(title = "Please select the parameters you want to apply to the scenarios and see the differences in the tab")
            ),
          icon = icon("filter")
        ),
        menuItem(
          "Health Economics",
          sliderInput(
            "out_discount_costs_slider",
            "Annual discount rate for costs",
            0,
            10,
            3.5,
            0.5,
            sep = "",
            ticks = FALSE,
            post  = " %"
          ),
          sliderInput(
            "out_discount_qalys_slider",
            "Annual discount rate for QALYs",
            0,
            10,
            1.5,
            0.5,
            sep = "",
            ticks = FALSE,
            post  = " %"
          ),
          numericInput(
            "out_wtp_box",
            "Willingness to pay (cost per QALY)",
            2e4L,
            0,
            1e5L,
            1e3L
          ),
          icon = icon("toolbox")
        ),
        br(),
        actionButton(
          "produce_report",
          "Produce report",
          icon = icon("file-contract"),
          style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
          #class = "btn-info",
          width = "auto"
        ) # ,
        # br(),
        # br(),
        # downloadButton(
        #   "save_analysis",
        #   "Archive analysis",
        #   icon("save"),
        #   style = "margin-left: 15px;"
        #   # class = "btn-info",
        #   # labelWidth = "100%"
        # ),
        # fileInput(
        #   "load_analysis",
        #   "",
        #   multiple = FALSE,
        #   accept = ".qs",
        #   placeholder = "",
        #   buttonLabel = "Load archived analysis"
        # ),
        # actionButton(
        #   "view_analysis",
        #   "View archived analysis",
        #   icon = icon("file-contract"),
        #   # style = "color: #fff; background-color: #337ab7; border-color: #2e6da4",
        #   #class = "btn-info",
        #   width = "auto"
        # ),
        # bookmarkButton(id = "bookmarkBtn"),
        # tags$style(type = 'text/css', "#bookmarkBtn { width:100%; margin-top: 25px;}")
      )
    ),
    dashboardBody(
      # main body ----
      tabItems(
        tabItem(
          # dashboard tab ----
          "dashboard", fluidPage(
          fluidRow(
            infoBoxOutput("most_cost_effective_box"),
            infoBoxOutput("most_equitable_box"),
            infoBoxOutput("most_effective_box")
          ),
          fluidRow(
            tags$head(tags$style("html, body {overflow: visible; }")),
            box(
              title = "Cost-effectiveness plane",
              solidHeader = TRUE,
              collapsible = FALSE,
              fluidRow(
                column(
                  1,
                  bs_modal(
                    id = "modal_cep",
                    title = "Explanation of this cost-effectiveness plane",
                    body = textOutput("info_ce_plane")
                  ),
                  shiny_iconlink() %>%
                    bs_attach_modal(id_modal = "modal_cep")
                ),
                #var_out_proc_use_ui("plot1"),
                #var_tt_creation_ui("plot1_tt"),

                column(10, offset = 1,
                       plotlyOutput("cep1")),
                column(
                  1,
                  offset = 10,
                  dropdown(
                    "Result display",
                    switchInput(
                      "res_display_cep1",
                      "Mean values",
                      FALSE,
                      onLabel = "Yes",
                      offLabel = "No"
                    ),

                    style = "unite",
                    status = "primary",
                    size = "sm",
                    icon = icon("cog", class = "fa-xs"),
                    width = "200px"
                    #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                  )
                )
              )
            ) ,
            box(
              title = "Equity plane",
              solidHeader = TRUE,
              collapsible = FALSE,
              fluidRow(
                column(
                  1,
                  bs_modal(
                    id = "modal_equ",
                    title = "Explications of this equity plane",
                    body = textOutput("info_equ_plane")
                  ),
                  shiny_iconlink() %>%
                    bs_attach_modal(id_modal = "modal_equ")
                ),
                column(10, offset = 1,
                       plotlyOutput("equ1")),
                column(
                  1,
                  offset = 10,
                  dropdown(
                    "Result display",
                    switchInput(
                      "res_display_equ1",
                      "Mean values",
                      FALSE,
                      onLabel = "Yes",
                      offLabel = "No"
                    ),
                    style = "unite",
                    status = "primary",
                    size = "sm",
                    right = TRUE,
                    icon = icon("cog", class = "fa-xs"),
                    width = "200px"
                    #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                  )
                )
              )
            ),
            fluidRow(
              box(
                title = "Explications",
                solidHeader = TRUE,
                collapsible = FALSE,
                width = 12,
                uiOutput("automated_text_descr")
              )
            )
          )
        )),
        tabItem(
          # cost-effectiveness tab ----
          "cost_effectiveness",
          fluidRow(
            box(
              title = "Cost-effectiveness plane",
              solidHeader = TRUE,
              collapsible = FALSE,
              fluidRow(
                column(
                  1,
                  bs_modal(
                    id = "modal_cep1_1",
                    title = "Explications of this Cost-effectiveness plane",
                    body = textOutput("info_ce1_plane")
                  ),
                  shiny_iconlink() %>%
                    bs_attach_modal(id_modal = "modal_cep1_1")
                ),
                column(10, offset = 1,
                       plotlyOutput("cep1_1")),
                column(
                  1,
                  offset = 10,
                  dropdown(
                    "Result display",
                    switchInput(
                      "res_display_cep1_1",
                      "Mean values",
                      FALSE,
                      onLabel = "Yes",
                      offLabel = "No"
                    ),

                    style = "unite",
                    status = "primary",
                    size = "sm",
                    icon = icon("cog", class = "fa-xs"),
                    width = "200px"
                    #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                  )
                )
              )
            ),
            box(
              title = "Notes",
              solidHeader = TRUE,
              collapsible = FALSE,
              height = 505,
              uiOutput("note_cost_eff")
            )
          ),
          fluidRow(
            tabBox(
              title = "Probability of",
              side = "right",
              selected = "cost-effective policy",
              id = "out_cep_p",
              tabPanel(
                "cost saving policy",
                bs_modal(
                  id = "modal_cep_p_cs",
                  title = "Explications of this Cost saving policy plane",
                  body = textOutput("info_cepcs_plane")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_cep_p_cs"),
                plotlyOutput("cep_p_cs")
              ),
              tabPanel(
                "cost-effective policy",
                bs_modal(
                  id = "modal_cep_p_ce",
                  title = "Explications of this Cost-effective policy plane",
                  body = textOutput("info_cepce_plane")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_cep_p_ce"),
                plotlyOutput("cep_p_ce")
              )
            ),
            box(
              title = "Cost-effectiveness over time",
              solidHeader = TRUE,
              collapsible = FALSE,
              bs_modal(
                id = "modal_cep_anim",
                title = "Explications of this Cost-effectiveness over time plane",
                body = textOutput("info_cep_anim_plane")
              ),
              shiny_iconlink() %>%
                bs_attach_modal(id_modal = "modal_cep_anim"),
              plotlyOutput("cep_anim")
            )
          )
        ),
        tabItem(
          # equity tab ----
          "equity",
                fluidRow(
                  tabBox(
                    title = "Equity planes",
                    side = "right",
                    selected = "absolute health inequalities",
                    id = "out_equ_plane",
                    tabPanel(
                      "relative health inequalities",
                      fluidRow(
                        column(
                          1,
                          bs_modal(
                            id = "modal_equ_rel",
                            title = "Explications of this relative health inequalities plane concerning equity",
                            body = textOutput("info_equ_rel_plane")
                          ),
                          shiny_iconlink() %>%
                            bs_attach_modal(id_modal = "modal_equ_rel")
                        ),
                        column(10, offset = 1,
                               plotlyOutput("equ_rel")),
                        column(
                          1,
                          offset = 10,
                          dropdown(
                            "Result display",
                            switchInput(
                              "res_display_equ_rel",
                              "Mean values",
                              FALSE,
                              onLabel = "Yes",
                              offLabel = "No"
                            ),
                            style = "unite",
                            status = "primary",
                            size = "sm",
                            icon = icon("cog", class = "fa-xs"),
                            width = "200px"
                            #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                          )
                        )
                      )
                    ),
                    tabPanel(
                      "absolute health inequalities",
                      fluidRow(
                        column(
                          1,
                          bs_modal(
                            id = "modal_equ1_1",
                            title = "Explications of this absolute health inequalities plane concerning equity",
                            body = textOutput("info_equ1_1_plane")
                          ),
                          shiny_iconlink() %>%
                            bs_attach_modal(id_modal = "modal_equ1_1")
                        ),
                        column(10, offset = 1,
                               plotlyOutput("equ1_1")),
                        column(
                          1,
                          offset = 10,
                          dropdown(
                            "Result display",
                            switchInput(
                              "res_display_equ1_1",
                              "Mean values",
                              FALSE,
                              onLabel = "Yes",
                              offLabel = "No"
                            ),
                            style = "unite",
                            status = "primary",
                            size = "sm",
                            icon = icon("cog", class = "fa-xs"),
                            width = "200px"
                            #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                          )
                        )
                      )
                    )
                  ),
                  box(
                    title = "Notes",
                    solidHeader = TRUE,
                    collapsible = FALSE,
                    height = 540,
                    uiOutput("note_health_ineq")
                  )
                ),
                fluidRow(
                  tabBox(
                    title = "Probability of equitable policy",
                    side = "right",
                    selected = "absolute health inequalities",
                    id = "out_equ_p",
                    tabPanel(
                      "relative health inequalities",
                      bs_modal(
                        id = "modal_equ_p_rel",
                        title = "Explications of this relative health inequalities plane concerning equitable policy",
                        body = textOutput("info_equ_p_rel_plane")
                      ),
                      shiny_iconlink() %>%
                        bs_attach_modal(id_modal = "modal_equ_p_rel"),
                      plotlyOutput("equ_p_rel")
                    ),
                    tabPanel(
                      "absolute health inequalities",
                      bs_modal(
                        id = "modal_equ_p_abs",
                        title = "Explications of this absolute health inequalities plane concerning equitable policy",
                        body = textOutput("info_equ_p_abs_plane")
                      ),
                      shiny_iconlink() %>%
                        bs_attach_modal(id_modal = "modal_equ_p_abs"),
                      plotlyOutput("equ_p_abs")
                    )
                  ),
                  tabBox(
                    title = "Equity over time",
                    side = "right",
                    selected = "absolute health inequalities",
                    id = "out_equ_anim",
                    tabPanel(
                      "relative health inequalities",
                      bs_modal(
                        id = "modal_equ_anim_rel",
                        title = "Explications of this relative health inequalities plane concerning equity over time",
                        body = textOutput("info_equ_anim_rel_plane")
                      ),
                      shiny_iconlink() %>%
                        bs_attach_modal(id_modal = "modal_equ_anim_rel"),
                      plotlyOutput("equ_anim_rel")
                    ),
                    tabPanel(
                      "absolute health inequalities",
                      bs_modal(
                        id = "modal_equ_anim_abs",
                        title = "Explications of this absolute health inequalities plane concerning equity over time",
                        body = textOutput("info_equ_anim_abs_plane")
                      ),
                      shiny_iconlink() %>%
                        bs_attach_modal(id_modal = "modal_equ_anim_abs"),
                      plotlyOutput("equ_anim_abs")
                    )
                  )
                )),
        tabItem(
          # effectiveness tab ----
          "effectiveness",
          fluidRow(
            tabBox(
              title = "Prevented or postponed",
              side = "right",
              selected = "Case-years",
              id = "out_pp_plane",
              tabPanel(
                "Deaths",
                bs_modal(
                  id = "modal_dpp_1",
                  title = "Explications",
                  body = textOutput("info_dpp_1_chart")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_dpp_1"),
                plotlyOutput("dpp_1")
              ),
              tabPanel(
                "Case-years",
                bs_modal(
                  id = "modal_cypp_1",
                  title = "Explications",
                  body = textOutput("info_cypp_1_chart")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_cppy_1"),
                plotlyOutput("cypp_1")
              ),
              tabPanel(
                "Cases",
                bs_modal(
                  id = "modal_cpp_1",
                  title = "Explications",
                  body = textOutput("info_cpp_1_chart")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_cpp_1"),
                plotlyOutput("cpp_1")
              )
            ),
            box(
              title = "Notes",
              solidHeader = TRUE,
              collapsible = FALSE,
              height = 525,
              uiOutput("note_effvnss")
            )
          ),
          fluidRow(
            tabBox(
              title = "Prevented or postponed over time",
              side = "right",
              selected = "Case-years",
              id = "out_ppy_plane",
              tabPanel(
                "Deaths",
                bs_modal(
                  id = "modal_dppy_spline",
                  title = "Explications",
                  body = textOutput("info_dppy_spline")
                ),
                shiny_iconlink() %>%
                  bs_attach_modal(id_modal = "modal_dppy_spline"),
                plotlyOutput("dppy_spline")
              ),
              tabPanel("Case-years",
                       fluidRow(
                         column(
                           1,
                           bs_modal(
                             id = "modal_cyppy_spline",
                             title = "Explications",
                             body = textOutput("info_cyppy_spline")
                           ),
                           shiny_iconlink() %>%
                             bs_attach_modal(id_modal = "modal_cyppy_spline")
                         ),
                         column(10, offset = 1,
                                plotlyOutput("cyppy_spline")),
                         column(
                           1,
                           offset = 10,
                           dropdown(
                             checkboxGroupButtons(
                               inputId = "out_diseases_select_cyppy",
                               label   = "Diseases",
                               choices = list(
                                 "Coronary heart disease" = "cypp_chd_cml",
                                 "Stroke" = "cypp_stroke_cml",
                                 "Post-stroke dementia" = "cypp_poststroke_dementia_cml",
                                 # "Hypertension" = "cypp_htn_cml",
                                 "Diabetes" = "cypp_t2dm_cml",
                                 "COPD" = "cypp_copd_cml",
                                 "Lung cancer" = "cypp_lung_ca_cml",
                                 "Colon cancer" = "cypp_colon_ca_cml",
                                 "Breast cancer" = "cypp_breast_ca_cml"
                               ),
                               selected = list(
                                 "cypp_chd_cml",
                                 "cypp_stroke_cml",
                                 "cypp_poststroke_dementia_cml",
                                 "cypp_t2dm_cml",
                                 # "cypp_htn_cml",
                                 "cypp_copd_cml",
                                 "cypp_lung_ca_cml",
                                 "cypp_colon_ca_cml",
                                 "cypp_breast_ca_cml"
                               ),
                               size = "xs",
                               direction = "vertical",
                               width = "100px"
                             ),
                             style = "unite",
                             status = "primary",
                             size = "sm",
                             up = TRUE,
                             icon = icon("cog", class = "fa-xs"),
                             width = "200px"
                             #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                           )
                         )
                       )),
              tabPanel(
                "Cases",
                fluidRow(
                  column(
                    1,
                    bs_modal(
                      id = "modal_cppy_spline",
                      title = "Explications",
                      body = textOutput("info_cppy_spline")
                    ),
                    shiny_iconlink() %>%
                      bs_attach_modal(id_modal = "modal_cppy_spline")
                  ),
                  column(10, offset = 1,
                         plotlyOutput("cppy_spline")),
                  column(
                    1,
                    offset = 10,
                    dropdown(
                      checkboxGroupButtons(
                        inputId = "out_diseases_select_cppy",
                        label   = "Diseases",
                        choices = list(
                          "Coronary heart disease" = "cpp_chd_cml",
                          "Stroke" = "cpp_stroke_cml",
                          "Post-stroke dementia" = "cpp_poststroke_dementia_cml",
                          # "Hypertension" = "cpp_htn_cml",
                          "Diabetes" = "cpp_t2dm_cml",
                          "COPD" = "cpp_copd_cml",
                          "Lung cancer" = "cpp_lung_ca_cml",
                          "Colon cancer" = "cpp_colon_ca_cml",
                          "Breast cancer" = "cpp_breast_ca_cml"
                        ),
                        selected = list(
                          "cpp_chd_cml",
                          "cpp_stroke_cml",
                          "cpp_poststroke_dementia_cml",
                          "cpp_t2dm_cml",
                          # "cpp_htn_cml",
                          "cpp_copd_cml",
                          "cpp_lung_ca_cml",
                          "cpp_colon_ca_cml",
                          "cpp_breast_ca_cml"
                        ),
                        size = "xs",
                        direction = "vertical",
                        width = "100px"
                      ),
                      style = "unite",
                      status = "primary",
                      size = "sm",
                      up = TRUE,
                      icon = icon("cog", class = "fa-xs"),
                      width = "200px"
                      #tooltip = tooltipOptions(title = "Click to chose the way you want to see the results !")
                    )
                  )
                ))
            ),
            tabBox(
              title = "Prevented or postponed over time",
              side = "right",
              selected = "Case-years",
              id = "case_years_staked_graph",
              tabPanel("Case-years",
                       fluidRow(
                         column(
                           1,
                           bs_modal(
                             id = "modal_case_years",
                             title = "Explications",
                             body = textOutput("info_case_years_stack")
                           ),
                           shiny_iconlink() %>%
                             bs_attach_modal(id_modal = "modal_case_year")
                         ),
                         column(10, offset = 1,
                                plotlyOutput("cyppy_stacked_area")),
                         column(
                           1,
                           offset = 10,
                           dropdown(
                             uiOutput("out_scenario_disease_select_cypp"),
                             style = "unite",
                             status = "primary",
                             size = "sm",
                             right = TRUE,
                             up = TRUE,
                             icon = icon("cog", class = "fa-xs"),
                             width = "200px"
                           )
                         )
                       )),
              tabPanel("Cases",
                       fluidRow(
                         column(
                           1,
                           bs_modal(
                             id = "modal_cases",
                             title = "Explications",
                             body = textOutput("info_cases_stack")
                           ),
                           shiny_iconlink() %>%
                             bs_attach_modal(id_modal = "modal_cases")
                         ),
                         column(10, offset = 1,
                                plotlyOutput("cppy_stacked_area")),
                         column(
                           1,
                           offset = 10,
                           dropdown(
                             uiOutput("out_scenario_disease_select_cpp"),
                             style = "unite",
                             status = "primary",
                             size = "sm",
                             right = TRUE,
                             up = TRUE,
                             icon = icon("cog", class = "fa-xs"),
                             width = "200px"
                           )
                         )
                       ))
              # , bs_modal(id = "modal_cpp_1", title = "Explications", body = textOutput("info_cpp_1_chart")),
              # shiny_iconlink() %>%
              #   bs_attach_modal(id_modal = "modal_cpp_1"), plotlyOutput("cpp_1"))

            )
          )
        ),
        tabItem(
          # wider_social_benefit ----
          "wider_social_benef",
          fluidRow(
            tabBox(
              title = "1st plane",
              side = "right",
              selected = "absolute health inequalities",
              id = "out_equ_plane",
              tabPanel("A 1st tab"),
              # , bs_modal(id = "modal_equ_rel", title = "Explications of this relative health inequalities plane concerning equity", body = textOutput("info_equ_rel_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ_rel"), plotlyOutput("equ_rel")),
              tabPanel("A 2nd tab")
              # , bs_modal(id = "modal_equ1_1", title = "Explications of this absolute health inequalities plane concerning equity", body = textOutput("info_equ1_1_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ1_1"), plotlyOutput("equ1_1"))
            ),
            box(
              title = "Notes",
              solidHeader = TRUE,
              collapsible = FALSE,
              uiOutput("note_social_benef")
              #p("Please do not use these results for any real-life application")
            )
          ),
          fluidRow(
            tabBox(
              title = "2nd plane",
              side = "right",
              selected = "absolute health inequalities",
              id = "out_equ_p",
              tabPanel("A 1st tab"),
              # , bs_modal(id = "modal_equ_p_rel", title = "Explications of this relative health inequalities plane concerning equitable policy", body = textOutput("info_equ_p_rel_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ_p_rel"), plotlyOutput("equ_p_rel")),
              tabPanel("A 2nd tab")
              # , bs_modal(id = "modal_equ_p_abs", title = "Explications of this absolute health inequalities plane concerning equitable policy", body = textOutput("info_equ_p_abs_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ_p_abs"), plotlyOutput("equ_p_abs"))
            ),
            tabBox(
              title = "3rd plane",
              side = "right",
              selected = "absolute health inequalities",
              id = "out_equ_anim",
              tabPanel("A 1st tab"),
              # , bs_modal(id = "modal_equ_anim_rel", title = "Explications of this relative health inequalities plane concerning equity over time", body = textOutput("info_equ_anim_rel_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ_anim_rel"), plotlyOutput("equ_anim_rel")),
              tabPanel("A 2nd tab")
              # , bs_modal(id = "modal_equ_anim_abs", title = "Explications of this absolute health inequalities plane concerning equity over time", body = textOutput("info_equ_anim_abs_plane")),
              #          shiny_iconlink() %>%
              #            bs_attach_modal(id_modal = "modal_equ_anim_abs"), plotlyOutput("equ_anim_abs"))
            )
          )
          #div(p("Dashboard tab content"))
        ),
        tabItem(
          # tables ----
          "summarised_results",
                box(
                  title = "Summarised results",
                  width = 12,
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  # column(4, offset = 8, uiOutput("out_columns_select_summary")),
                  column(12,
                         downloadButton("download_summary","Download summary"),
                         DT::dataTableOutput("tbl_summary"))
                )
                ),
        tabItem(
          "raw_output",
          box(
            title = "Raw results",
            width = 12,
            solidHeader = TRUE,
            collapsible = TRUE,
            # column(4, offset = 8, uiOutput("out_columns_select_raw")),
            column(12,
                   downloadButton("download_raw",
                                  "Download raw outputs")
                   # ,
                   # DT::dataTableOutput("tbl_raw")
                   )
          )
        )
      )
    )
  ),
  icon = icon("poll")
)
