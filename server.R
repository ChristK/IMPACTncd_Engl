## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
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

server = function(input, output, session) {
  # enableBookmarking(store = "server")

  # restrict friendly names to 32 characters at max
  shinyjs::runjs("$('#friendly_name_sc1').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc2').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc3').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc4').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc5').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc6').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc7').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc8').attr('maxlength',32)")
  shinyjs::runjs("$('#friendly_name_sc9').attr('maxlength',32)")

# Dynamic scenario tab addition -------------------------------------------
  output$scenario_tabs <- renderUI({
    Tabs <- as.list(rep(0, input$scenarios_number_slider))
    for (i in 0:length(Tabs)) {
      Tabs[i] <- lapply(paste0("Scenario ", i), tabPanel, value = i)
    }

    #Tabs <- lapply(paste("Layer",0:input$subClust,sep=" "), tabPanel)
    do.call(tabsetPanel, c(Tabs, id = "level"))
  })


# Reactive function to choose scenario colours and shapes, and use their value
# according to their name

  colsymb <- reactive({
      setnames(transpose(as.data.table(lapply(seq_len(input$scenarios_number_slider), function(i) {
        c(input[[paste0("col_sc", i)]],
          input[[paste0("friendly_name_sc", i)]],
          input[[paste0("symbol_sc", i)]])
      }))), c("colour", "names", "symbol")
        )[names %in% input$inout_scenario_select, ]
  })

  # Reactive function to create an empty synthpop

  empty_synthpop <- reactive({
    parameters <- fromGUI_prune(reactiveValuesToList(input))
    design$update_fromGUI(parameters)
    SynthPop$new(0, design)
  })


  # firstFrame <- vapplyseq_len(input$scenarios_number_slider), function(i) isTRUE(tr[["frame"]] %in% frameNames[[1]]), logical(1))
  # input$scenarios_number_slider[firstFrame] <- p$x$frames[[1]]$data

  #p$x$data <- c(c(p$x$data, colours), NULL)

  # callModule(module = createPlot, id = "cep1")




  observeEvent(input$scenarios_number_slider, {
# Render uptake/Px tables -------------------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "render_scenario_tables.R")))), local = TRUE
             )$value
    }

    # Modify ,min/max scenario init year ==------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "dynamic_scenario_minmax_year.R")))), local = TRUE
             )$value
    }

    # Load/save scenario  -----------------------------------------------------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                      readLines(file.path("server", "save_load_scenarios.R")))), local = TRUE
             )$value
    }


    # Conditionally disable inputs -------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                       readLines(file.path("server", "show_hide_logic.R")))), local = TRUE
      )$value
    }

    # Prev/Next buttons -------
    for (i in seq_len(input$scenarios_number_slider)) {
      source(exprs = str2expression(gsub("_sc1",  paste0("_sc", i),
                                       readLines(file.path("server", "prev_next_buttons.R")))), local = TRUE
      )$value
    }
  },
  ignoreNULL = FALSE, ignoreInit = FALSE)


  # Check if synthpop exists and generate if not ----------------------------
   observeEvent(input$locality_select_validator, {
  # See https://stackoverflow.com/questions/50165443/async-process-blocking-r-shiny-app
  # and https://github.com/rstudio/promises/issues/23
     plan(list(
       # 1 worker for GUI and one for parallel population write.
       tweak(multiprocess, workers = 2L), # outer loop
       tweak(multiprocess, workers = design$sim_prm$clusternumber) # inner loop
       # Because from the outer loop 1 worker never parallelises the inner loop,
       # the inner loop need not be cpu/2. But RAM usage is excess.
     ))
    design$update_fromGUI(reactiveValuesToList(input))
    # if no are selected from gui, defaults to England
    future({
        SynthPop$
        new(0, design)$
        write_synthpop(1:(design$sim_prm$iteration_n * design$sim_prm$n_synthpop_aggregation))
    })
    NULL # This is necessary!! after future when there is no promise to return
  })


# Run simulation ----
  out <-
    eventReactive(input[[paste0("run_simulation_sc", input$scenarios_number_slider)]],
      {
        plan(multiprocess, workers = design$sim_prm$clusternumber)
        parameters <- fromGUI_prune(reactiveValuesToList(input))
        design$update_fromGUI(parameters)

        qsave(reactiveValuesToList(input, all.names = TRUE),
          file.path(design$sim_prm$output_dir, "input.qs"))
        qsave(parameters, file.path(design$sim_prm$output_dir, "parameters.qs"))

        withProgress(message = 'Running IMPACTncd_Engl model.',
          detail = 'This may take a couple of minutes...',
          value = 0,
          {
            # TODO remove before release

            if (file.exists(file.path(design$sim_prm$output_dir, "results.fst"))) {
              # out <- read_fst(file.path(design$sim_prm$output_dir, "results.fst"), as.data.table = TRUE)
              # out
              file.remove(file.path(design$sim_prm$output_dir, "results.fst"))
              run_simulation(parameters, design, FALSE)
            } else {
              run_simulation(parameters, design, FALSE)
            }
          })
      },
      ignoreNULL = TRUE,
      ignoreInit = TRUE)

  # Report button ----
  out_report <- eventReactive(input$produce_report, {
    plan(multiprocess, workers = design$sim_prm$clusternumber)

    parameters <- fromGUI_prune(reactiveValuesToList(input))
    design$update_fromGUI(parameters)

    qsave(reactiveValuesToList(input, all.names = TRUE),
      file.path(design$sim_prm$output_dir, "input.qs"))
    qsave(parameters, file.path(design$sim_prm$output_dir, "parameters.qs"))

    # progress$inc(1/n, detail = paste("Doing part", i))
    withProgress(message = 'Running IMPACTncd_Engl model.',
      detail = 'This may take a couple of hours...',
      value = 0,
      {
        if (file.exists(file.path(design$sim_prm$output_dir, "results.fst")))
          file.remove(file.path(design$sim_prm$output_dir, "results.fst"))
        run_simulation(parameters, design, TRUE)
      })
  },
    ignoreNULL = FALSE, # Otherwise require double-click
    ignoreInit = TRUE)

  # Export logs ----
  output$download_logs_gui <- downloadHandler(


    filename = "logs.yaml",

    content = function(file) {

      log_files <- list.files(file.path(design$sim_prm$output_dir, "logs"), full.names = TRUE)
      logs <- lapply(log_files, readLines)
      names(logs) <-
        gsub(
          paste0(file.path(design$sim_prm$output_dir, "logs/"), "|.txt$"),
          "", log_files)
      if (file.exists(file.path(design$sim_prm$output_dir, "input.qs"))) {
        logs$input <-
          qread(file.path(design$sim_prm$output_dir, "input.qs"))
      }
      if (file.exists(file.path(design$sim_prm$output_dir, "parameters.qs"))) {
        logs$parameters <-
          qread(file.path(design$sim_prm$output_dir, "parameters.qs"))
      }
      if (file.exists(file.path(design$sim_prm$output_dir, "times.txt"))) {
        logs$times<-
          readLines(file.path(design$sim_prm$output_dir, "times.txt"))
      }

      write_yaml(logs, file = file)
    }
  )

  # Delete logs
  observeEvent(input$delete_logs_gui, {
    session$sendCustomMessage(type = 'testmessage',
      message = 'Log files will be permanently deleted!')
    unlink(file.path(design$sim_prm$output_dir, "logs"), recursive = TRUE)
  })

  # Add sysload info -----
  output$sysload <- renderPrint({
    # Re-execute this reactive expression after 2000 milliseconds
    invalidateLater(2000)
    isolate(sysLoad())
  })

  # Manipulate synthpops ----
  output$checksum <- renderPrint({
    empty_synthpop()$get_checksum()
  })

  output$synthpop_info <- renderPrint({
    # Re-execute this reactive expression after 10 seconds
    invalidateLater(5e3)
    empty_synthpop()$count_synthpop()
  })

  # Delete synthpops
  observeEvent(input$check_synthpops_gui, {
    empty_synthpop()$
    check_integridy(remove_malformed = TRUE, check_checksum = FALSE)$
      delete_incomplete_synthpop(check_checksum = FALSE)
  })

  observeEvent(input$check_synthpops_gui_check_checksum, {
    empty_synthpop()$
      check_integridy(remove_malformed = TRUE, check_checksum = TRUE)$
      delete_incomplete_synthpop(check_checksum = TRUE)
  })

  observeEvent(input$delete_synthpops_gui_spare_primer, {
    empty_synthpop()$delete_synthpop(NULL, spare_primer = TRUE, check_checksum = TRUE)
  })

  observeEvent(input$delete_synthpops_gui_check_checksum, {
    empty_synthpop()$delete_synthpop(NULL, spare_primer = FALSE, check_checksum = TRUE)
  })

  observeEvent(input$delete_synthpops_gui_all, {
    empty_synthpop()$delete_synthpop(NULL, spare_primer = FALSE, check_checksum = FALSE)
  })


  # Save archived analysis ----
  # myBookmarks <- reactiveValues(urlDF = NULL)
  # observeEvent(input$bookmarkBtn, {
  #   session$doBookmark()
  # })
  #
  # output$save_analysis <- downloadHandler(
  #
  #   filename = function() "Analysis.qs",
  #
  #   content = function(file) {
  #     parameters <- fromGUI_prune(reactiveValuesToList(input))
  #     input <- reactiveValuesToList(input, all.names = TRUE)
  #     results <- out()
  #     qsavem(input, parameters, results,
  #       file = file, nthreads = design$sim_prm$clusternumber)
  #   }
  # )

  # Load archived analysis ----
  # observeEvent(input$load_analysis, {
  #
  #   inFile <- input$load_analysis
  #
  #   if (is.null(inFile)) return(NULL)
  #   qload(inFile$datapath)
  # })
  #
  # # View archived analysis ----
  # out_archived <- eventReactive(input$view_analysis, {
  #   lapply(names(input),
  #     function(x) session$sendInputMessage(x, list(value = input[[x]]))
  #   )
  #   write_fst(results, "./output/results2.fst")
  #   setDT(results)
  #   },
  #   ignoreNULL = FALSE, # Otherwise require double-click
  #   ignoreInit = TRUE)

# change to output panel
  observeEvent(input[[paste0("run_simulation_sc", input$scenarios_number_slider)]], {
    check_locality <- length(input$locality_select) > 0
    check_baseline <-
      any(
        input$baseline_sc1,
        input$baseline_sc2,
        input$baseline_sc3,
        input$baseline_sc4,
        input$baseline_sc5,
        input$baseline_sc6,
        input$baseline_sc7,
        input$baseline_sc8,
        input$baseline_sc9
      )
    # TODO What if more than one baseline scenario is selected. It needs to take
    # into account ensemble scenarios.

    if (!check_locality) {

      createAlert(
        session,
        paste0("alert_locality_sc", input$scenarios_number_slider),
        "alert_locality",
        title = "Oops... area is missing!",
        content = "Please select an area to simulate from the 'Simulation parameters' tab.",
        append = TRUE,
        style = "warning")

      } else closeAlert(session, "alert_locality")

     if (!check_baseline) {

      createAlert(
        session,
        paste0("alert_baseline_sc", input$scenarios_number_slider),
        "alert_baseline",
        title = "Oops...no baseline scenario!",
        content = "Please select a scenario as the baseline one. The baseline button is at top of each scenario page.",
        append =  TRUE,
        style = "warning")

    } else closeAlert(session, "alert_baseline")

      if (check_locality && check_baseline) {

        updateTabsetPanel(session, "inTabset",
          selected = "output_panel")

      }
  })

# Close alerts when error is corrected
  observeEvent(input$locality_select, {
    check_locality <- length(input$locality_select) > 0
    if (check_locality) closeAlert(session, "alert_locality")
  }
    )

  observeEvent({
    input$baseline_sc1
    input$baseline_sc2
    input$baseline_sc3
    input$baseline_sc4
    input$baseline_sc5
    input$baseline_sc6
    input$baseline_sc7
    input$baseline_sc8
    input$baseline_sc9
    1}, {

      check_baseline <-
        any(
          input$baseline_sc1,
          input$baseline_sc2,
          input$baseline_sc3,
          input$baseline_sc4,
          input$baseline_sc5,
          input$baseline_sc6,
          input$baseline_sc7,
          input$baseline_sc8,
          input$baseline_sc9
        )

      if (check_baseline) closeAlert(session, "alert_baseline")
    })



# Output tab --------------------------------------------------------
  output$most_cost_effective_box <- renderInfoBox({
    infoBox(
      "Most Cost-Effective", most_cost_effective(out_proc()), icon = icon("pound-sign"),
      color = "aqua", fill = FALSE
    )
  })

  output$most_equitable_box <- renderInfoBox({
    infoBox(
      "Most Equitable", most_equitable(out_proc_qimd()), icon = icon("balance-scale"),
      color = "purple"
    )
  })

  output$most_effective_box <- renderInfoBox({
    infoBox(
      "Most Effective", most_effective(out_proc()), icon = icon("heart"),
      color = "yellow"
    )
  })


  diff_year <- reactive({
    input$inout_year_slider - input$simulation_period_slider[1]
  })


# Automated graph explanation text  ------------------------------------------------------


  output$automated_text_descr <- renderUI({

    order_bcr_sc_nam <- # bcr = benefit:cost ratio
      reactive({
        most_benefit_cost_ratio(out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          10L
        )})

    most_bcr_sc_val <-
      reactive({
        round(benefit_cost_ratio_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          order_bcr_sc_nam()[1L]
        ), 2L)
      })

    most_nmb_sc_val <-
      reactive({
        net_monetary_benefit_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          order_bcr_sc_nam()[1L]
        )
      })

    icer_1st_vs_2nd <-
      reactive({
        round(
          scn_icer_cml(
            out_proc(),
            order_bcr_sc_nam()[1],
            input$health_econ_perspective_checkbox
          ) -
            ifelse(is.na(scn_icer_cml(
              out_proc(),
              order_bcr_sc_nam()[2],
              input$health_econ_perspective_checkbox
            )), 0, scn_icer_cml(
              out_proc(),
              order_bcr_sc_nam()[2],
              input$health_econ_perspective_checkbox
            ))
        )
      })

    HTML(
      paste0(
        "With time horizon of ",
        diff_year(),
        " years, starting in year ",
        input$simulation_period_slider[1],
        ", the most cost-effective scenario was ",
        order_bcr_sc_nam()[1L],
      " scenario. ",
      "It had a benefit:cost ratio (including the value of health gains) of £",
        most_bcr_sc_val(),
      " for every £1 spent.",
    br(), br(),


        "This scenario had a ", tolower(input$health_econ_perspective_checkbox),
        " incremental cost-effectiveness ratio of £",
        icer_1st_vs_2nd(),
        " per QALY when compared to the next best scenario. ", br(), br(),


       "The order of the scenarios, starting from the most cost-effective was ",
       head(order_bcr_sc_nam(), -1),
      " and finally, ",
      last(order_bcr_sc_nam()),
      "."
      ))
  })


  output$note_cost_eff <- renderUI({
    # TODO find a way to avoid repeating these function as they are also used in
    # output$automated_text_descr
    order_bcr_sc_nam <- # bcr = benefit:cost ratio
      reactive({
        most_benefit_cost_ratio(out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          10L
        )})

    most_bcr_sc_val <-
      reactive({
        round(benefit_cost_ratio_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          order_bcr_sc_nam()[1L]
        ), 2L)
      })

    most_nmb_sc_val <-
      reactive({
        net_monetary_benefit_value(
          out_proc(),
          input$health_econ_perspective_checkbox,
          input$out_wtp_box,
          order_bcr_sc_nam()[1L]
        )
      })

    icer_1st_vs_2nd <-
      reactive({
        round(
          scn_icer_cml(
            out_proc(),
            order_bcr_sc_nam()[1],
            input$health_econ_perspective_checkbox
          ) -
            ifelse(is.na(scn_icer_cml(
              out_proc(),
              order_bcr_sc_nam()[2],
              input$health_econ_perspective_checkbox
            )), 0, scn_icer_cml(
              out_proc(),
              order_bcr_sc_nam()[2],
              input$health_econ_perspective_checkbox
            ))
        )
      })

    HTML(
      paste0(
        "With time horizon of ",
        diff_year(),
        " years, starting in year ",
        input$simulation_period_slider[1],
        ": and based on valuing QALYs at ",
        input$out_wtp_box,
        ", the scenario with the greatest ",
        tolower(input$health_econ_perspective_checkbox),
        " benefit:cost ratio was ",
        order_bcr_sc_nam()[1L],
        " scenario which had a benefit:cost ratio of ",
        most_bcr_sc_val(),
        " for every £1 spent, and a total net monetary benefit of £",
        most_nmb_sc_val(),
        " from a ", tolower(input$health_econ_perspective_checkbox), ".",
        br(),
        br(),
        "This scenario had a ", tolower(input$health_econ_perspective_checkbox),
        " incremental cost-effectiveness ratio of £",
        icer_1st_vs_2nd(),
        " per QALY when compared to the next best scenario (",
        order_bcr_sc_nam()[2],
        " ).",
        br(),
        br(),
        "Ranking the scenarios from most to least cost-effective, the order was: ",
        paste0(paste(head(order_bcr_sc_nam(), -1), collapse = ", "), " and ",
          tail(order_bcr_sc_nam(), 1), ".")
      )
    )
  })

   #TODO: update with new functions for absolute and relative equitable

   output$note_health_ineq <- renderUI({
     HTML(paste0(most_equitable_rel(out_proc_qimd()), " scenario had the biggest impact in terms of reducing relative inequalities, reducing the relative index of inequalities by ", reduce_rel_index_ineq(), ".", br(), br(),
                 last(most_equitable_abs(out_proc_qimd())), " scenario did the least harm to absolute inequalities, increasing the absolute index of inequalities by ", increase_abs_index_ineq(), "."))

   })


   output$note_effvnss <- renderUI({
     HTML(paste0("The most effective scenario in terms of disease cases prevented and postponed was ", most_effective(out_proc()), " scenario which prevented '...' total disease cases,
                 made up of '...' cases of CVD, and '...' cases of other diseases"))
   })


   output$note_social_benef <- renderUI({
     HTML(paste0("As well as healthcare benefits, ", most_equitable(out_proc_qimd())," scenario would produce additional social care cost savings of ",
                 social_care_cost_sav(out_proc())[most_equitable(out_proc_qimd())]$V1, ", productivity benefits of ", prod_benef(out_proc())[most_equitable(out_proc_qimd())]$V1, " (including earnings
and the value of household productivity) and informal care cost savings of", inform_care_cost_sav(out_proc())[most_equitable(out_proc_qimd())]$V1, br()))
   })

 output$info_ce_plane <- renderText({
   "Needs to be implemented"
 })

 output$info_equ_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_ce1_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cepcs_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cepce_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_cep_anim_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_rel_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ1_1_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_p_rel_plane <- renderText({
   "The future explications for this specific graph"
 })

 output$info_equ_p_abs_plane <- renderText({
   "The future explications for this specific graph"
 })

output$info_equ_anim_rel_plane <- renderText({
  "The future explications for this specific graph"
})

output$info_equ_anim_abs_plane <- renderText({
  "The future explications for this specific graph"
})

 # user_inputs <- reactive(reactiveValuesToList(input))


  source(file.path("server", "output_inputs.R"),  local = TRUE)$value


# Tooltips ----------------------------------------------------------------
  # source(file.path("server", "tooltips.R"),  local = TRUE)$value


  # outputOptions(output, "nrows", suspendWhenHidden = FALSE)


# last_rank <- function(){
#   lr <- rank_cost_effective(out_proc())[length(rank_cost_effective(out_proc()))]
#   lr <- paste0(" and ", last(rank_cost_effective(out_proc())))
# }

}

