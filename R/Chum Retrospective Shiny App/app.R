# Chum retrospective policy analysis shiny app 
# Haley Oleynik
# Feb 2026

library(shiny)
library(tidyverse)
library(readr)
library(readxl)

# ==============================
# ---- DATA + MODEL SETUP -----
# ==============================

fp <- "Fraser chum Oleynik retrospective_v2.xlsx"

chum <- read_csv("walters_chum_retrospective_data.csv")

dat <- read_excel(fp, sheet = "Orig data", skip = 3) %>%
  transmute(
    Year        = `Brood year`,
    total_stock = `total Stock, hr from Grant then Brittany`,
    U_hist      = `Historial Ut`,
    wt          = wt,
    p3          = rp3,
    p4          = rp4,
    p5          = rp5
  ) %>%
  filter(!is.na(Year)) %>%
  arrange(Year)

fit <- lm(lnrs ~ Spawners, data = chum)
aa <- unname(coef(fit)[[1]])
bb <- unname(-coef(fit)[[2]])

chum <- chum %>%
  mutate(wt = lnrs - (aa - bb * Spawners))


# ==============================
# ---- SIMULATOR --------------
# ==============================

simulate_policy <- function(dat, umethod = c("hist","constU","hcr"),
                            retroU = 0.5,
                            cslope = 0.8,
                            nlim   = 500000,
                            aa, bb,
                            util_exp = 0.6) {
  
  umethod <- match.arg(umethod)
  n <- nrow(dat)
  
  Run      <- numeric(n)
  U        <- numeric(n)
  Spawners <- numeric(n)
  Recruits <- numeric(n)
  Yield    <- numeric(n)
  
  for (i in seq_len(n)) {
    
    if (i <= 4) {
      Run[i] <- dat$total_stock[i]
    } else {
      run_i <- 0
      if (i - 3 >= 1) run_i <- run_i + dat$p3[i - 3] * Recruits[i - 3]
      if (i - 4 >= 1) run_i <- run_i + dat$p4[i - 4] * Recruits[i - 4]
      if (i - 5 >= 1) run_i <- run_i + dat$p5[i - 5] * Recruits[i - 5]
      Run[i] <- run_i
    }
    
    if (umethod == "hist") {
      U[i] <- dat$U_hist[i]
    } else if (umethod == "constU") {
      U[i] <- retroU
    } else {
      U[i] <- max(0, cslope * (Run[i] - nlim) / Run[i])
    }
    
    Spawners[i] <- Run[i] * (1 - U[i])
    Recruits[i] <- Spawners[i] * exp(aa - bb * Spawners[i] + dat$wt[i])
    Yield[i] <- Run[i] - Spawners[i]
  }
  
  tibble(
    Year = dat$Year,
    Run = Run,
    U = U,
    Spawners = Spawners,
    Recruits = Recruits,
    Yield = Yield,
    Utility = Yield^util_exp
  )
}

summarize_policy <- function(sim_tbl) {
  sim_tbl %>%
    summarise(
      total_yield = sum(Yield, na.rm = TRUE),
      utility     = sum(Utility, na.rm = TRUE)
    )
}

# ==============================
# ---- UI ----------------------
# ==============================

ui <- fluidPage(
  
  titlePanel("Fraser Chum Retrospective Policy Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      sliderInput("retroU", "Constant U",
                  min = 0, max = 1, value = 0.5, step = 0.1),
      
      numericInput("nlim", "Conservation minimum",
                   value = 500000),
      helpText("Escapement threshold below which harvest rate is zero"),
      
      sliderInput("cslope", "HCR slope",
                  min = 0, max = 1, value = 0.8, step = 0.05),
      helpText("Controls how aggressively harvest rate increases once run size exceeds the conservation minimum
           Also sets the maximum U at very large run sizes."),
      
      sliderInput("util_exp", "Utility exponent",
                  min = 0.1, max = 1.5, value = 0.6, step = 0.05),
      helpText(
        HTML("Determines risk preference: &gamma; &lt; 1 penalizes variability and favors stable harvests; &gamma; = 1 maximizes total yield.")
      )
    ),
    
    mainPanel(
      verbatimTextOutput("policySummary"),
      verbatimTextOutput("lossText"),
      verbatimTextOutput("lossHCRText"),
      
      hr(),
      
      plotOutput("yieldTS", height = 250),
      plotOutput("Uts", height = 250),
      plotOutput("spawnTS", height = 250),
      
      hr(),
      
      plotOutput("utilityPlot", height = 300),
      plotOutput("Uplot", height = 350)
      
      
    )
    
  )
)

# ==============================
# ---- SERVER ------------------
# ==============================

server <- function(input, output) {
  
  sims <- reactive({
    
    hist <- simulate_policy(dat, "hist",
                            aa = aa, bb = bb,
                            util_exp = input$util_exp)
    
    const <- simulate_policy(dat, "constU",
                             retroU = input$retroU,
                             aa = aa, bb = bb,
                             util_exp = input$util_exp)
    
    maxY <- simulate_policy(dat, "hcr",
                            cslope = input$cslope,
                            nlim   = input$nlim,
                            aa = aa, bb = bb,
                            util_exp = input$util_exp)
    
    maxU <- simulate_policy(
      dat,
      umethod = "hcr",
      cslope = 0.7763355348064025,
      nlim   = 487692.7340110076,
      aa = aa,
      bb = bb,
      util_exp = input$util_exp
    )
    
    list(hist = hist, const = const, maxY = maxY, maxU = maxU)
  })
  
  policy_summary <- reactive({
    
    pol <- bind_rows(
      Historical = summarize_policy(sims()$hist),
      ConstU     = summarize_policy(sims()$const),
      MaxY       = summarize_policy(sims()$maxY),
      MaxU       = summarize_policy(sims()$maxU),
      .id = "Policy"
    )
    pol <- pol %>%
      mutate(
        rel_yield   = total_yield / max(total_yield, na.rm = TRUE),
        rel_utility = utility     / max(utility, na.rm = TRUE)
      )
    
    pol
  })
  
  
  
  
  historical_series <- reactive({
    
    dat %>%
      mutate(
        Hist_Spawners = total_stock * (1 - U_hist),
        Hist_Yield    = total_stock - Hist_Spawners
      ) %>%
      select(Year, Hist_Yield, U_hist, Hist_Spawners)
  })
  
  output$lossText <- renderText({
    
    pol <- policy_summary()
    
    hist_yield  <- pol$total_yield[pol$Policy == "Historical"]
    const_yield <- pol$total_yield[pol$Policy == "ConstU"]
    
    loss <- const_yield - hist_yield
    
    paste0("Loss (ConstU − Historical): ",
           format(round(loss, 0), big.mark = ","))
  })
  
  output$lossHCRText <- renderText({
    
    pol <- policy_summary()
    
    hist_yield <- pol$total_yield[pol$Policy == "Historical"]
    hcr_yield  <- pol$total_yield[pol$Policy == "MaxY"]
    
    loss <- hcr_yield - hist_yield
    
    paste0("Loss (HCR − Historical): ",
           format(round(loss, 0), big.mark = ","))
  })
  
  
  
  # ---- Utility vs Yield plot ----
  
  
  output$utilityPlot <- renderPlot({
    
    pol <- policy_summary()
    
    ggplot(pol,
           aes(rel_yield, rel_utility,
               colour = Policy,
               label = Policy)) +
      
      geom_point(size = 4) +
      geom_text(nudge_y = 0.02, show.legend = FALSE) +
      
      scale_colour_manual(values = c(
        Historical = "steelblue",
        ConstU     = "goldenrod",
        MaxY       = "black",
        MaxU       = "red"
      )) +
      
      labs(
        x = "Relative Yield",
        y = "Relative Utility",
        colour = "Policy",
        title = "Utility vs Yield Tradeoff"
      ) +
      
      theme_minimal(base_size = 14) +
      theme(legend.position = "right")
  })
  
  
  
  
  # ---- U vs Run size plot ----
  
  output$Uplot <- renderPlot({
    
    hist  <- sims()$hist  %>% mutate(Policy = "Historical")
    const <- sims()$const %>% mutate(Policy = "ConstU")
    maxY  <- sims()$maxY  %>% mutate(Policy = "HCR")
    
    combined <- bind_rows(hist, const, maxY)
    
    ggplot(combined,
           aes(Run, U, colour = Policy)) +
      
      geom_point(
        data = filter(combined, Policy == "Historical"),
        alpha = 0.6
      ) +
      
      geom_line(
        data = filter(combined, Policy != "Historical"),
        size = 1.2
      ) +
      
      scale_colour_manual(
        values = c(
          Historical = "steelblue",
          ConstU     = "goldenrod",
          HCR        = "black"
        )
      ) +
      
      labs(
        x = "Run size",
        y = "Exploitation rate",
        colour = "Policy",
        title = "U vs Run Size"
      ) +
      
      theme_minimal(base_size = 14) +
      theme(legend.position = "top")
  })
  
  # yield plot 
  output$yieldTS <- renderPlot({
    
    hist  <- sims()$hist  %>% mutate(Policy = "Historical")
    const <- sims()$const %>% mutate(Policy = "ConstU")
    hcr   <- sims()$maxY  %>% mutate(Policy = "HCR")
    
    combined <- bind_rows(hist, const, hcr)
    
    ggplot(combined, aes(Year, Yield, colour = Policy)) +
      geom_line(size = 1.2) +
      scale_colour_manual(values = c(
        Historical = "steelblue",
        ConstU     = "goldenrod",
        HCR        = "black"
      )) +
      labs(y = "Yield",
           title = "Yield by Policy",
           colour = "Policy") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  
  # harvest rates 
  output$Uts <- renderPlot({
    
    hist  <- sims()$hist  %>% mutate(Policy = "Historical")
    const <- sims()$const %>% mutate(Policy = "ConstU")
    hcr   <- sims()$maxY  %>% mutate(Policy = "HCR")
    
    combined <- bind_rows(hist, const, hcr)
    
    ggplot(combined, aes(Year, U, colour = Policy)) +
      geom_line(size = 1.2) +
      scale_colour_manual(values = c(
        Historical = "steelblue",
        ConstU     = "goldenrod",
        HCR        = "black"
      )) +
      labs(y = "Harvest rate (U)",
           title = "Harvest Rate by Policy",
           colour = "Policy") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  
  # spawner timeseries 
  output$spawnTS <- renderPlot({
    
    hist  <- sims()$hist  %>% mutate(Policy = "Historical")
    const <- sims()$const %>% mutate(Policy = "ConstU")
    hcr   <- sims()$maxY  %>% mutate(Policy = "HCR")
    
    combined <- bind_rows(hist, const, hcr)
    
    ggplot(combined, aes(Year, Spawners, colour = Policy)) +
      geom_line(size = 1.2) +
      scale_colour_manual(values = c(
        Historical = "steelblue",
        ConstU     = "goldenrod",
        HCR        = "black"
      )) +
      labs(y = "Spawners",
           title = "Spawner Abundance by Policy",
           colour = "Policy") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top")
  })
  
  
  output$policySummary <- renderPrint({
    policy_summary()
  })
  
  # as a table 
  # output$policyTable <- renderTable({
  #   
  #   policy_summary() %>%
  #     mutate(
  #       total_yield = round(total_yield, 0),
  #       utility     = round(utility, 0)
  #     )
  #   
  # }, striped = TRUE, bordered = TRUE)
  
  
}

shinyApp(ui, server)

