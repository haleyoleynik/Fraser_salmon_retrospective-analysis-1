# Sockeye retrospective shiny app 
# Haley Oleynik 
# February 2026

library(shiny)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)

# ---------------- CONFIG ----------------
FIT_YEARS <- 1952:2019
LAG_YEARS <- 4

# ---------------- LOAD DATA ----------------
obs <- read_csv("Walters_model_all-stocks.csv") %>%
  rename(
    Year            = `Year`,
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    DBE             = `DBE`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  arrange(Stock, Year)

# ---------------- RETRO FUNCTION ----------------
run_retro_model <- function(dat,
                            retroU,
                            useretro,
                            yrretro) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs  = pmin(0.999, Catch / RunJacks),
      ENS     = pmin(1, pmax(0.0001,
                             AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S       = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
  
  ra <- coef(fit)[1]
  rb <- -coef(fit)[2]
  
  obs3 <- obs2 %>%
    mutate(wt = lnR_S - (ra - rb * AdultEscapement))
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retroU_vec <- if (useretro)
    ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
  else
    obs3$Ut_obs
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retroR[i] <- retroS[j] *
      exp(ra - rb * retroS[j] + obs3$wt[j])
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR = retroR,
      retroU = retroU_vec,
      retroS = retroS,
      retroC = retroC
    )
}

# ---------------- UI ----------------
ui <- fluidPage(
  titlePanel("Fraser Sockeye Multi-Stock Retrospective Model"),
  
  sidebarLayout(
    sidebarPanel(
      
      sliderInput("retroU",
                  "Retrospective harvest rate",
                  min = 0,
                  max = 1,
                  value = 0.3,
                  step = 0.05),
      
      checkboxInput("useretro", "Use retrospective harvest rate?",
                    value = TRUE),
      
      numericInput("yrretro", "Retrospective start year",
                   value = 1990),
      
      hr(),
      
      selectInput("stock_choice",
                  "Select stock",
                  choices = c("All Stocks", sort(unique(obs$Stock))),
                  selected = "All Stocks")
    ),
    
    mainPanel(
      h4("Total catch"),
      verbatimTextOutput("total_summary"),
      
      plotOutput("catch_plot", height = 300),
      plotOutput("return_plot", height = 300),
      plotOutput("esc_plot", height = 300),
      
      h4("Time Series"),
      tableOutput("retro_ts_table")
      
      #h4("Stock Summary"),
      #tableOutput("stock_table")
    )
  )
)

# ---------------- SERVER ----------------
server <- function(input, output, session) {
  
  model_all <- reactive({
    obs %>%
      group_by(Stock) %>%
      group_modify(~ run_retro_model(.x,
                                     retroU = input$retroU,
                                     useretro = input$useretro,
                                     yrretro = input$yrretro)) %>%
      ungroup()
  })
  
  filtered_model <- reactive({
    out <- model_all()
    
    if (input$stock_choice == "All Stocks") {
      out
    } else {
      out %>% filter(Stock == input$stock_choice)
    }
  })
  
  # ---- Stock Summary ----
  output$stock_table <- renderTable({
    
    out <- filtered_model()
    
    if (input$stock_choice == "All Stocks") {
      out %>%
        group_by(Stock) %>%
        summarise(
          hist_catch  = sum(Catch,  na.rm = TRUE),
          retro_catch = sum(retroC, na.rm = TRUE),
          lost_catch  = retro_catch - hist_catch,
          .groups = "drop"
        )
    } else {
      out %>%
        summarise(
          hist_catch  = sum(Catch,  na.rm = TRUE),
          retro_catch = sum(retroC, na.rm = TRUE),
          lost_catch  = retro_catch - hist_catch
        )
    }
  })
  
  # ---- Retro Time Series Table ----
  output$retro_ts_table <- renderTable({
    
    df <- summed_ts() %>%
      select(Year, retroR, retroC, Run = RunJacks) %>%
      arrange(Year) %>%
      mutate(
        retroR = round(retroR, 0),
        retroC = round(retroC, 0),
        Run = round(Run, 0)
      )
    
    df
  })
  
  
  # ---- Total Summary ----
  output$total_summary <- renderPrint({
    
    out <- filtered_model()
    
    hist_total  <- sum(out$Catch,  na.rm = TRUE)
    retro_total <- sum(out$retroC, na.rm = TRUE)
    
    cat("Historical Catch: ",
        format(round(hist_total,0), big.mark=","), "\n")
    
    cat("Retrospective Catch: ",
        format(round(retro_total,0), big.mark=","), "\n")
    
    cat("Difference: ",
        format(round(retro_total - hist_total,0), big.mark=","), "\n")
  })
  
  
  # ---- Aggregate Timeseries (summed across stocks) ----
# summed_ts <- reactive({
# 
#   out <- filtered_model()
# 
#   out %>%
#     group_by(Year) %>%
#     summarise(
#       Catch = sum(Catch, na.rm=TRUE),
#       retroC = sum(retroC, na.rm=TRUE),
#       RunJacks = sum(RunJacks, na.rm=TRUE),
#       retroR = sum(retroR, na.rm=TRUE),
#       AdultEscapement = sum(AdultEscapement, na.rm=TRUE),
#       retroS = sum(retroS, na.rm=TRUE),
#       .groups="drop"
#     ) %>%
#     arrange(Year) %>%
#     mutate(
#       AdultReturn = lead(RunJacks, LAG_YEARS)
#     )
# })
  
  summed_ts <- reactive({
    
    out <- filtered_model()
    
    out %>%
      group_by(Year) %>%
      summarise(
        Catch = sum(Catch, na.rm=TRUE),
        retroC = sum(retroC, na.rm=TRUE),
        RunJacks = sum(RunJacks, na.rm=TRUE),
        retroR = sum(retroR, na.rm=TRUE),
        AdultEscapement = sum(AdultEscapement, na.rm=TRUE),
        retroS = sum(retroS, na.rm=TRUE),
        .groups="drop"
      ) %>%
      arrange(Year) %>%
      mutate(
        AdultReturn = lead(RunJacks, LAG_YEARS)
      )
  })
  
  
  # Catch plot
  output$catch_plot <- renderPlot({
    df <- summed_ts()
    
    plot_df <- df %>%
      select(Year, Catch, retroC) %>%
      pivot_longer(-Year)
    
    ggplot(plot_df,
           aes(Year, value, color=name, linetype=name)) +
      geom_line(linewidth=1.2) +
      labs(title="Catch",
           y="Catch",
           color = "") +
      guides(linetype = "none") +
      scale_color_manual(values = c("#4682B4","#FF4500")) +
      theme_minimal()
  })
  
  # Adult return plot
  output$return_plot <- renderPlot({
    df <- summed_ts()
    
    plot_df <- df %>%
      select(Year, AdultReturns = RunJacks, retroR) %>%
      pivot_longer(-Year)
    
    ggplot(plot_df,
           aes(Year, value, color=name,linetype=name)) +
      geom_line(linewidth=1.2) +
      labs(title="Return",
           y="Returns",
           color = "") +
      guides(linetype = "none") +
      scale_color_manual(values = c("#4682B4","#FF4500")) +
      theme_minimal()
  })
  
  # Escapement plot
  output$esc_plot <- renderPlot({
    df <- summed_ts()
    
    plot_df <- df %>%
      select(Year, AdultEscapement, retroS) %>%
      pivot_longer(-Year)
    
    ggplot(plot_df,
           aes(Year, value, color=name,linetype=name)) +
      geom_line(linewidth=1.2) +
      labs(title="Spawners",
           y="Spawners",
           color = "") +
      guides(linetype = "none") +
      scale_color_manual(values = c("#4682B4","#FF4500")) +
      theme_minimal()
  })
}

shinyApp(ui, server)
