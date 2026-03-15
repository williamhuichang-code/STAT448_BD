library(shiny)
library(ggplot2)
library(MASS)  # for mvrnorm, optional

ui <- fluidPage(
  titlePanel("Ridge vs OLS: Constraint Geometry"),
  sidebarLayout(
    sidebarPanel(
      h4("OLS estimate (β̂)"),
      sliderInput("b1_ols", "β̂₁ (OLS)", min = -3, max = 3, value = 2.0, step = 0.1),
      sliderInput("b2_ols", "β̂₂ (OLS)", min = -3, max = 3, value = 1.5, step = 0.1),
      hr(),
      h4("XᵀX matrix (predictor structure)"),
      sliderInput("rho", "Correlation ρ between X₁, X₂",
                  min = -0.95, max = 0.95, value = 0.7, step = 0.05),
      sliderInput("scale", "Relative scale of X₂ vs X₁",
                  min = 0.3, max = 3.0, value = 1.0, step = 0.1),
      hr(),
      h4("Regularization"),
      sliderInput("lambda", "λ (ridge penalty)",
                  min = 0.01, max = 5, value = 1.0, step = 0.05),
      hr(),
      h4("Display"),
      checkboxInput("show_lasso", "Also show Lasso constraint (diamond)", value = TRUE),
      checkboxInput("show_path", "Show ridge path as λ varies", value = FALSE)
    ),
    mainPanel(
      plotOutput("ridge_plot", height = "550px"),
      hr(),
      fluidRow(
        column(4, wellPanel(
          h5("Ridge estimate"),
          verbatimTextOutput("ridge_vals")
        )),
        column(4, wellPanel(
          h5("Shrinkage factor"),
          verbatimTextOutput("shrinkage")
        )),
        column(4, wellPanel(
          h5("Constraint radius t"),
          verbatimTextOutput("radius_info")
        ))
      )
    )
  )
)

server <- function(input, output) {
  
  # Build XᵀX from user inputs
  XtX <- reactive({
    s <- input$scale
    r <- input$rho
    # XᵀX = [[1, r*s], [r*s, s^2]]  (scaled correlation structure)
    matrix(c(1, r * s, r * s, s^2), nrow = 2)
  })
  
  b_ols <- reactive({
    c(input$b1_ols, input$b2_ols)
  })
  
  # Ridge solution: β_ridge = (XᵀX + λI)⁻¹ XᵀX β_ols
  # (since XᵀX β_ols = Xᵀy by definition)
  b_ridge <- reactive({
    A <- XtX()
    lam <- input$lambda
    solve(A + lam * diag(2)) %*% A %*% b_ols()
  })
  
  # Constraint radius: ||β_ridge||₂
  t_radius <- reactive({
    sqrt(sum(b_ridge()^2))
  })
  
  # Generate RSS ellipse contours
  # RSS(β) = (β - β_ols)ᵀ XᵀX (β - β_ols) = c
  ellipse_contour <- function(c_val, n = 200) {
    A <- XtX()
    b0 <- b_ols()
    eig <- eigen(A)
    # Parametric ellipse: β = β_ols + U D^{-1/2} [cos,sin] * sqrt(c)
    U <- eig$vectors
    D <- diag(1 / sqrt(eig$values))
    theta <- seq(0, 2 * pi, length.out = n)
    pts <- t(b0 + sqrt(c_val) * U %*% D %*% rbind(cos(theta), sin(theta)))
    data.frame(x = pts[, 1], y = pts[, 2], level = c_val)
  }
  
  # RSS value at the ridge solution (innermost touching ellipse)
  c_ridge <- reactive({
    br <- b_ridge()
    b0 <- b_ols()
    d <- br - b0
    as.numeric(t(d) %*% XtX() %*% d)
  })
  
  output$ridge_plot <- renderPlot({
    b0 <- b_ols()
    br <- as.vector(b_ridge())
    t_r <- t_radius()
    c_r <- c_ridge()
    lam <- input$lambda
    
    # Ellipse levels: touching one + 3 outer ones
    levels <- c_r * c(1, 2.5, 5, 9)
    ell_data <- do.call(rbind, lapply(levels, ellipse_contour))
    ell_data$level <- factor(ell_data$level)
    
    # Ridge constraint circle
    theta <- seq(0, 2 * pi, length.out = 300)
    circle_df <- data.frame(
      x = t_r * cos(theta),
      y = t_r * sin(theta)
    )
    
    # Lasso diamond
    lasso_df <- data.frame(
      x = c(t_r, 0, -t_r, 0, t_r),
      y = c(0, t_r, 0, -t_r, 0)
    )
    
    # Ridge path as lambda varies
    path_df <- NULL
    if (input$show_path) {
      lambdas <- seq(0.01, 8, length.out = 80)
      A <- XtX()
      path_pts <- t(sapply(lambdas, function(l) {
        as.vector(solve(A + l * diag(2)) %*% A %*% b0)
      }))
      path_df <- data.frame(x = path_pts[, 1], y = path_pts[, 2])
    }
    
    p <- ggplot() +
      # RSS ellipses
      geom_path(data = ell_data,
                aes(x = x, y = y, group = level, alpha = level),
                color = "#C0392B", linewidth = 0.8) +
      scale_alpha_manual(values = c(1, 0.7, 0.5, 0.3), guide = "none") +
      
      # Constraint circle fill
      geom_polygon(data = circle_df, aes(x = x, y = y),
                   fill = "#2980B9", alpha = 0.15, color = "#1A5276", linewidth = 1) +
      
      # Lasso diamond (optional)
      {if (input$show_lasso)
        geom_path(data = lasso_df, aes(x = x, y = y),
                  color = "#1A5276", linewidth = 1, linetype = "dashed")
      } +
      
      # Ridge path (optional)
      {if (!is.null(path_df))
        geom_path(data = path_df, aes(x = x, y = y),
                  color = "#8E44AD", linewidth = 0.8, linetype = "dotted")
      } +
      
      # OLS estimate
      geom_point(aes(x = b0[1], y = b0[2]),
                 color = "#C0392B", size = 4, shape = 16) +
      annotate("text", x = b0[1] + 0.12, y = b0[2] + 0.12,
               label = "OLS estimate", color = "#C0392B", size = 3.8, hjust = 0) +
      
      # Ridge estimate
      geom_point(aes(x = br[1], y = br[2]),
                 color = "#1A5276", size = 4, shape = 16) +
      annotate("text", x = br[1] - 0.12, y = br[2] - 0.18,
               label = "Ridge estimate", color = "#1A5276", size = 3.8, hjust = 1) +
      
      # Line from origin to ridge estimate
      annotate("segment", x = 0, y = 0, xend = br[1], yend = br[2],
               color = "#1A5276", linewidth = 0.5, linetype = "dotted") +
      
      # Axes
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
      
      labs(
        title = paste0("Ridge regression geometry  (λ = ", lam, ",  ρ = ", input$rho, ")"),
        x = expression(beta[1]),
        y = expression(beta[2])
      ) +
      coord_fixed() +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14)
      )
    
    # Annotations for constraint regions
    p <- p +
      annotate("text", x = -t_r * 0.6, y = -t_r * 0.6,
               label = paste0("Ridge ball\n||β||₂ ≤ ", round(t_r, 2)),
               color = "#1A5276", size = 3.2, hjust = 0.5) +
      {if (input$show_lasso)
        annotate("text", x = t_r * 0.75, y = t_r * 0.75,
                 label = "Lasso diamond\n||β||₁ ≤ t",
                 color = "#1A5276", size = 3.2, hjust = 0)
      }
    
    print(p)
  })
  
  output$ridge_vals <- renderText({
    br <- as.vector(b_ridge())
    paste0("β₁ = ", round(br[1], 4), "\nβ₂ = ", round(br[2], 4))
  })
  
  output$shrinkage <- renderText({
    br <- as.vector(b_ridge())
    b0 <- b_ols()
    s1 <- round(br[1] / b0[1], 3)
    s2 <- round(br[2] / b0[2], 3)
    paste0("β₁: shrunk to ", round(s1 * 100, 1), "% of OLS\n",
           "β₂: shrunk to ", round(s2 * 100, 1), "% of OLS")
  })
  
  output$radius_info <- renderText({
    paste0("||β_ridge||₂ = ", round(t_radius(), 4),
           "\nRSS at tangency = ", round(c_ridge(), 4))
  })
}

shinyApp(ui, server)