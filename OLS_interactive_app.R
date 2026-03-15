library(shiny)
library(ggplot2)

# ── helpers ──────────────────────────────────────────────────────────────────

build_XtX <- function(rho, scale) {
  matrix(c(1, rho * scale, rho * scale, scale^2), nrow = 2)
}

ridge_solution <- function(XtX, b_ols, lambda) {
  as.vector(solve(XtX + lambda * diag(2)) %*% XtX %*% b_ols)
}

lasso_solution <- function(b_ols, t_val) {
  # Soft-thresholding for orthogonal X (exact for uncorrelated predictors)
  # For general X: iterative coordinate descent
  sign(b_ols) * pmax(abs(b_ols) - t_val, 0)
}

elasticnet_solution <- function(XtX, b_ols, lambda, alpha) {
  # alpha=1 → Lasso, alpha=0 → Ridge
  # Elastic net: (XtX + lambda*(1-alpha)*I)^{-1} * (XtX * b_ols - lambda*alpha*sign(b))
  # Iterative coordinate descent (2 coords, closed form per step)
  b <- b_ols
  for (iter in seq_len(200)) {
    b_old <- b
    for (j in 1:2) {
      r <- XtX[j, -j] %*% b[-j]
      z <- (XtX[j, j] * b_ols[j] - r)
      denom <- XtX[j, j] + lambda * (1 - alpha)
      b[j] <- sign(z) * max(abs(z) - lambda * alpha / 2, 0) / denom
    }
    if (max(abs(b - b_old)) < 1e-8) break
  }
  b
}

ellipse_points <- function(XtX, b_ols, c_val, n = 300) {
  eig <- eigen(XtX)
  U   <- eig$vectors
  D   <- diag(1 / sqrt(eig$values))
  th  <- seq(0, 2 * pi, length.out = n)
  pts <- t(b_ols + sqrt(c_val) * U %*% D %*% rbind(cos(th), sin(th)))
  data.frame(x = pts[, 1], y = pts[, 2])
}

rss_at <- function(XtX, b_ols, b) {
  d <- b - b_ols
  as.numeric(t(d) %*% XtX %*% d)
}

constraint_circle <- function(r, n = 300) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = r * cos(th), y = r * sin(th))
}

constraint_diamond <- function(t, n = 200) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = t * cos(th) * abs(cos(th))^0 * sign(cos(th)),
    y = t * sin(th) * abs(sin(th))^0 * sign(sin(th))
  )
  # exact Lasso diamond corners
  data.frame(
    x = c(t, 0, -t, 0, t),
    y = c(0, t, 0, -t, 0)
  )
}

constraint_elasticnet <- function(t, alpha, n = 400) {
  th <- seq(0, 2 * pi, length.out = n)
  # |x|^p + |y|^p = t  where p is between 1 (lasso) and 2 (ridge)
  # For elastic net constraint shape approximation:
  # alpha controls mix: alpha=1 → diamond, alpha=0 → circle
  # We use the actual elastic net ball: alpha*||b||_1 + (1-alpha)*||b||_2^2/2 <= t
  pts <- matrix(NA, n, 2)
  for (i in seq_len(n)) {
    angle <- th[i]
    # Find r such that alpha*r*(|cos|+|sin|) + (1-alpha)*r^2/2 = t
    # Solve: (1-alpha)/2 * r^2 + alpha*(|cos|+|sin|)*r - t = 0
    ac <- abs(cos(angle)) + abs(sin(angle))
    a_coef <- (1 - alpha) / 2
    b_coef <- alpha * ac
    c_coef <- -t
    if (a_coef < 1e-9) {
      r <- t / b_coef
    } else {
      disc <- b_coef^2 - 4 * a_coef * c_coef
      r <- (-b_coef + sqrt(disc)) / (2 * a_coef)
    }
    pts[i, ] <- c(r * cos(angle), r * sin(angle))
  }
  data.frame(x = pts[, 1], y = pts[, 2])
}

shrinkage_path <- function(XtX, b_ols, method, alpha = 0.5, n = 60) {
  lambdas <- exp(seq(log(0.01), log(10), length.out = n))
  pts <- t(sapply(lambdas, function(lam) {
    switch(method,
           Ridge       = ridge_solution(XtX, b_ols, lam),
           Lasso       = elasticnet_solution(XtX, b_ols, lam, alpha = 1),
           ElasticNet  = elasticnet_solution(XtX, b_ols, lam, alpha = alpha)
    )
  }))
  data.frame(x = pts[, 1], y = pts[, 2])
}

# ── colour palette ────────────────────────────────────────────────────────────
COL <- list(
  ols    = "#C0392B",
  ridge  = "#1D6FB8",
  lasso  = "#1A7A3A",
  enet   = "#7B2FBE",
  ell    = "#C0392B",
  bg     = "#E8E8E8",
  grid   = "#D0D0D0",
  text   = "#555555",
  white  = "#1A1A1A"
)

METHODS <- c("Ridge", "Lasso", "ElasticNet")

METHOD_COL <- c(Ridge = COL$ridge, Lasso = COL$lasso, ElasticNet = COL$enet)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML(paste0("
    @import url('https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=Inter:wght@300;400;500;600&display=swap');
    html, body { overflow-x:hidden; } .container-fluid { padding:0 !important; max-width:100% !important; } body { background:#E8E8E8; color:#2A2A2A; font-family:'Inter',sans-serif; margin:0; }
    .well { background:#F2F2F2; border:1px solid #C8C8C8; border-radius:12px; padding:16px; }
    .shiny-input-container label { color:#555555; font-size:12px; font-weight:500; letter-spacing:.05em; text-transform:uppercase; }
    .irs--shiny .irs-bar { background:", COL$ridge, "; }
    .irs--shiny .irs-handle { background:#FFFFFF; border-color:#AAAAAA; }
    .irs--shiny .irs-line { background:#C8C8C8; }
    .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background:#F2F2F2; color:#555555; font-size:10px; }
    h1,h2,h3,h4,h5 { font-family:'Space Mono',monospace; color:#1A1A1A; }
    .card { background:#F2F2F2; border:1px solid #C8C8C8; border-radius:12px; padding:14px 18px; margin-bottom:12px; }
    .metric-label { font-size:10px; text-transform:uppercase; letter-spacing:.08em; color:#777777; margin-bottom:2px; }
    .metric-val   { font-family:'Space Mono',monospace; font-size:13px; color:#1A1A1A; }
    .method-tag { display:inline-block; padding:3px 10px; border-radius:20px; font-size:11px; font-weight:600; margin-right:6px; }
    .tag-ridge { background:rgba(29,111,184,.15); color:#1254A0; border:1px solid rgba(29,111,184,.5); }
    .tag-lasso { background:rgba(26,122,58,.15);  color:#145C2A; border:1px solid rgba(26,122,58,.5); }
    .tag-enet  { background:rgba(123,47,190,.15); color:#5A1E9A; border:1px solid rgba(123,47,190,.5); }
    .tag-ols   { background:rgba(192,57,43,.15);  color:#8B1A10; border:1px solid rgba(192,57,43,.5); }
    .section-title { font-size:10px; text-transform:uppercase; letter-spacing:.1em; color:#888888; margin-bottom:8px; border-bottom:1px solid #D8D8D8; padding-bottom:6px; }
    input[type=checkbox] { accent-color:", COL$ridge, "; }
    .shiny-plot-output { border-radius:12px; overflow:hidden; }
    #drag_hint { font-size:11px; color:#888888; text-align:center; margin-top:4px; font-style:italic; }
  ")))),
  div(style = "width:100%; padding:12px 16px; box-sizing:border-box;",
      
      div(style = "margin-bottom:12px;",
          h2(style = "margin:0; font-size:20px; letter-spacing:-.02em;",
             "Regularization Geometry Explorer"),
          p(style = "color:#777777; font-size:12px; margin:3px 0 0;",
            "Ridge \u00b7 Lasso \u00b7 Elastic Net \u2014 constraint shapes, RSS contours & shrinkage paths")
      ),
      
      fluidRow(style="margin:0;",
               
               column(2, style="padding:0 6px 0 0;",
                      div(class = "card",
                          div(class = "section-title", "OLS Estimate"),
                          p(style="font-size:10px;color:#888888;margin:0 0 6px;",
                            "Click the plot to move the OLS estimate"),
                          sliderInput("b1", "Beta1 (OLS)", min=-3, max=3, value=2.0, step=0.05, width="100%"),
                          sliderInput("b2", "Beta2 (OLS)", min=-3, max=3, value=1.5, step=0.05, width="100%")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Predictor Structure (XtX)"),
                          sliderInput("rho",   "Correlation rho",    min=-.95, max=.95, value=.65, step=.05, width="100%"),
                          sliderInput("scale", "Scale ratio s2/s1",  min=.3,  max=3,   value=1.0, step=.05, width="100%")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Penalty"),
                          sliderInput("lambda", "Lambda (strength)", min=.00, max=5, value=1.0, step=.05, width="100%"),
                          sliderInput("alpha",  "Alpha (L1 ratio - Elastic Net)",
                                      min=0, max=1, value=.5, step=.05, width="100%"),
                          p(style="font-size:10px;color:#888888;margin:3px 0 0;",
                            "alpha=1: Lasso  |  alpha=0: Ridge")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Display"),
                          checkboxInput("show_ridge",    span(class="method-tag tag-ridge","Ridge"), TRUE),
                          checkboxInput("show_lasso",    span(class="method-tag tag-lasso","Lasso"), TRUE),
                          checkboxInput("show_enet",     span(class="method-tag tag-enet", "Elastic Net"), TRUE),
                          checkboxInput("show_path",     "Shrinkage paths", FALSE),
                          checkboxInput("show_ellipses", "RSS contours", TRUE),
                          checkboxInput("show_axes_zero","Zero axes", TRUE)
                      )
               ),
               
               column(8, style="padding:0 6px;",
                      div(style="position:relative;",
                          plotOutput("main_plot",
                                     height = "calc(100vh - 120px)",
                                     click  = "plot_click"),
                          div(id="drag_hint", "Click anywhere on the plot to move the OLS estimate")
                      )
               ),
               
               column(2, style="padding:0 0 0 6px;",
                      div(class = "card",
                          div(class = "section-title", "Estimates"),
                          uiOutput("estimates_ui")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Shrinkage (% of OLS)"),
                          uiOutput("shrinkage_ui")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Constraint Radii"),
                          uiOutput("radii_ui")
                      ),
                      div(class = "card",
                          div(class = "section-title", "Key Insight"),
                          uiOutput("insight_ui")
                      )
               )
      )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # Allow clicking plot to move OLS estimate
  observeEvent(input$plot_click, {
    req(input$plot_click)
    x <- round(max(min(input$plot_click$x, 3), -3), 2)
    y <- round(max(min(input$plot_click$y, 3), -3), 2)
    updateSliderInput(session, "b1", value = x)
    updateSliderInput(session, "b2", value = y)
  })
  
  # Reactive computations
  XtX_r    <- reactive(build_XtX(input$rho, input$scale))
  b_ols_r  <- reactive(c(input$b1, input$b2))
  lam      <- reactive(input$lambda)
  alp      <- reactive(input$alpha)
  
  b_ridge_r <- reactive(ridge_solution(XtX_r(), b_ols_r(), lam()))
  b_lasso_r <- reactive(elasticnet_solution(XtX_r(), b_ols_r(), lam(), alpha = 1))
  b_enet_r  <- reactive(elasticnet_solution(XtX_r(), b_ols_r(), lam(), alpha = alp()))
  
  t_ridge_r <- reactive(sqrt(sum(b_ridge_r()^2)))
  t_lasso_r <- reactive(sum(abs(b_lasso_r())))
  t_enet_r  <- reactive({
    b <- b_enet_r(); alp() * sum(abs(b)) + (1 - alp()) / 2 * sum(b^2)
  })
  
  # Main plot
  output$main_plot <- renderPlot({
    b0  <- b_ols_r()
    XtX <- XtX_r()
    br  <- b_ridge_r()
    bl  <- b_lasso_r()
    be  <- b_enet_r()
    
    # Axis limits
    all_pts <- rbind(b0, br, bl, be, c(-3,-3), c(3,3))
    lim <- c(-3.2, 3.2)
    
    p <- ggplot() +
      theme_void() +
      theme(
        plot.background  = element_rect(fill = COL$bg, color = NA),
        panel.background = element_rect(fill = COL$bg, color = NA),
        panel.grid.major = element_line(color = "#D0D0D0", linewidth = .4),
        panel.grid.minor = element_line(color = "#DCDCDC", linewidth = .2),
        axis.text  = element_text(color = "#555555", size = 10, family = "mono"),
        axis.title = element_text(color = "#1A1A1A", size = 13, family = "mono"),
        axis.ticks = element_line(color = "#BBBBBB"),
        plot.margin = margin(16, 16, 16, 16)
      ) +
      scale_x_continuous(limits = lim, breaks = seq(-3,3,1), expand = c(0,0)) +
      scale_y_continuous(limits = lim, breaks = seq(-3,3,1), expand = c(0,0)) +
      coord_fixed() +
      xlab(expression(beta[1])) + ylab(expression(beta[2]))
    
    # Zero-axis highlights
    if (input$show_axes_zero) {
      p <- p +
        geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = .8) +
        geom_vline(xintercept = 0, color = "#AAAAAA", linewidth = .8)
    }
    
    # Shrinkage paths
    if (input$show_path) {
      if (input$show_ridge) {
        path_r <- shrinkage_path(XtX, b0, "Ridge")
        p <- p + geom_path(data=path_r, aes(x,y), color=COL$ridge, linewidth=.7,
                           linetype="dotted", alpha=.7)
      }
      if (input$show_lasso) {
        path_l <- shrinkage_path(XtX, b0, "Lasso")
        p <- p + geom_path(data=path_l, aes(x,y), color=COL$lasso, linewidth=.7,
                           linetype="dotted", alpha=.7)
      }
      if (input$show_enet) {
        path_e <- shrinkage_path(XtX, b0, "ElasticNet", alpha=alp())
        p <- p + geom_path(data=path_e, aes(x,y), color=COL$enet, linewidth=.7,
                           linetype="dotted", alpha=.7)
      }
    }
    
    # RSS ellipses
    if (input$show_ellipses) {
      c0 <- rss_at(XtX, b0, br)   # innermost touching ridge
      for (mult in c(1.8, 3.5, 6.5, 11)) {
        ell <- ellipse_points(XtX, b0, c0 * mult)
        p <- p + geom_path(data=ell, aes(x,y),
                           color=COL$ell, linewidth=.5, alpha=.25 + .05*(6/mult))
      }
      # Innermost ellipse (touching ridge constraint)
      ell_inner <- ellipse_points(XtX, b0, c0)
      p <- p + geom_path(data=ell_inner, aes(x,y),
                         color=COL$ell, linewidth=1, alpha=.9)
    }
    
    # Constraint shapes
    if (input$show_ridge) {
      circ <- constraint_circle(t_ridge_r())
      p <- p +
        geom_polygon(data=circ, aes(x,y), fill=COL$ridge, alpha=.08,
                     color=COL$ridge, linewidth=1)
    }
    if (input$show_lasso) {
      t_l <- sum(abs(bl))
      diam <- constraint_diamond(t_l)
      p <- p +
        geom_polygon(data=diam, aes(x,y), fill=COL$lasso, alpha=.08,
                     color=COL$lasso, linewidth=1, linetype="solid")
    }
    if (input$show_enet) {
      t_e <- alp() * sum(abs(be)) + (1-alp())/2 * sum(be^2)
      enet_shape <- constraint_elasticnet(t_e, alp())
      p <- p +
        geom_polygon(data=enet_shape, aes(x,y), fill=COL$enet, alpha=.08,
                     color=COL$enet, linewidth=1, linetype="dashed")
    }
    
    # Estimate points + leader lines from origin
    if (input$show_ridge) {
      p <- p +
        annotate("segment", x=0, y=0, xend=br[1], yend=br[2],
                 color=COL$ridge, linewidth=.5, linetype="dotted", alpha=.6) +
        geom_point(aes(x=br[1], y=br[2]), color=COL$ridge, size=4, shape=16) +
        geom_point(aes(x=br[1], y=br[2]), color="white",   size=1.5, shape=16) +
        annotate("text", x=br[1]+.13, y=br[2]+.13, label="Ridge",
                 color=COL$ridge, size=3.5, hjust=0, fontface="bold", family="mono")
    }
    if (input$show_lasso) {
      p <- p +
        annotate("segment", x=0, y=0, xend=bl[1], yend=bl[2],
                 color=COL$lasso, linewidth=.5, linetype="dotted", alpha=.6) +
        geom_point(aes(x=bl[1], y=bl[2]), color=COL$lasso, size=4, shape=16) +
        geom_point(aes(x=bl[1], y=bl[2]), color="white",   size=1.5, shape=16) +
        annotate("text", x=bl[1]+.13, y=bl[2]-.18, label="Lasso",
                 color=COL$lasso, size=3.5, hjust=0, fontface="bold", family="mono")
    }
    if (input$show_enet) {
      p <- p +
        annotate("segment", x=0, y=0, xend=be[1], yend=be[2],
                 color=COL$enet, linewidth=.5, linetype="dotted", alpha=.6) +
        geom_point(aes(x=be[1], y=be[2]), color=COL$enet, size=4, shape=16) +
        geom_point(aes(x=be[1], y=be[2]), color="white",  size=1.5, shape=16) +
        annotate("text", x=be[1]-.13, y=be[2]+.15, label="EN",
                 color=COL$enet, size=3.5, hjust=1, fontface="bold", family="mono")
    }
    
    # OLS estimate
    p <- p +
      geom_point(aes(x=b0[1], y=b0[2]), color=COL$ols, size=6, shape=16) +
      geom_point(aes(x=b0[1], y=b0[2]), color="white",  size=2,  shape=16) +
      annotate("text", x=b0[1]+.14, y=b0[2]+.14, label="OLS",
               color=COL$ols, size=3.8, hjust=0, fontface="bold", family="mono")
    
    print(p)
  }, bg = COL$bg)
  
  # ── Side panel outputs ─────────────────────────────────────────────────────
  fmt2 <- function(x) formatC(x, digits=3, format="f")
  
  output$estimates_ui <- renderUI({
    b0 <- b_ols_r(); br <- b_ridge_r(); bl <- b_lasso_r(); be <- b_enet_r()
    tagList(
      div(style="margin-bottom:10px;",
          div(class="metric-label", span(class="method-tag tag-ols","OLS")),
          div(class="metric-val", paste0("β₁=", fmt2(b0[1]), "  β₂=", fmt2(b0[2])))
      ),
      if (input$show_ridge) div(style="margin-bottom:10px;",
                                div(class="metric-label", span(class="method-tag tag-ridge","Ridge")),
                                div(class="metric-val", paste0("β₁=", fmt2(br[1]), "  β₂=", fmt2(br[2])))
      ),
      if (input$show_lasso) div(style="margin-bottom:10px;",
                                div(class="metric-label", span(class="method-tag tag-lasso","Lasso")),
                                div(class="metric-val", paste0("β₁=", fmt2(bl[1]), "  β₂=", fmt2(bl[2])))
      ),
      if (input$show_enet) div(
        div(class="metric-label", span(class="method-tag tag-enet",
                                       paste0("EN (α=", input$alpha, ")"))),
        div(class="metric-val", paste0("β₁=", fmt2(be[1]), "  β₂=", fmt2(be[2])))
      )
    )
  })
  
  output$shrinkage_ui <- renderUI({
    b0 <- b_ols_r(); br <- b_ridge_r(); bl <- b_lasso_r(); be <- b_enet_r()
    
    # Safe per-coefficient shrinkage string
    pct1 <- function(b_val, b0_val) {
      if (is.na(b0_val) || abs(b0_val) < 1e-10) {
        if (abs(b_val) < 1e-10) "0→0" else "—"
      } else {
        paste0(round(b_val / b0_val * 100, 1), "%")
      }
    }
    pct_str <- function(b) paste0("β₁: ", pct1(b[1], b0[1]),
                                  "  β₂: ", pct1(b[2], b0[2]))
    
    # Bar width based on mean absolute shrinkage (safe)
    bar_width <- function(b) {
      denom <- pmax(abs(b0), 1e-10)
      mean(pmin(abs(b) / denom, 1)) * 100
    }
    bar <- function(b, col) {
      w <- round(bar_width(b), 1)
      div(style = paste0(
        "height:4px;border-radius:2px;margin:2px 0 6px;",
        "width:", w, "%;background:", col, ";"))
    }
    
    tagList(
      if (input$show_ridge) div(
        div(class="metric-label", span(class="method-tag tag-ridge", "Ridge")),
        div(class="metric-val", pct_str(br)),
        bar(br, COL$ridge)
      ),
      if (input$show_lasso) div(
        div(class="metric-label", span(class="method-tag tag-lasso", "Lasso")),
        div(class="metric-val", pct_str(bl)),
        bar(bl, COL$lasso)
      ),
      if (input$show_enet) div(
        div(class="metric-label", span(class="method-tag tag-enet", "Elastic Net")),
        div(class="metric-val", pct_str(be)),
        bar(be, COL$enet)
      )
    )
  })
  
  output$radii_ui <- renderUI({
    tagList(
      if (input$show_ridge) div(style="margin-bottom:8px;",
                                div(class="metric-label", "Ridge  ||β||₂"),
                                div(class="metric-val", fmt2(t_ridge_r()))
      ),
      if (input$show_lasso) div(style="margin-bottom:8px;",
                                div(class="metric-label", "Lasso  ||β||₁"),
                                div(class="metric-val", fmt2(sum(abs(b_lasso_r()))))
      ),
      if (input$show_enet) div(
        div(class="metric-label", paste0("EN  α||β||₁ + (1-α)/2·||β||₂²")),
        div(class="metric-val", fmt2(t_enet_r()))
      )
    )
  })
  
  output$insight_ui <- renderUI({
    bl <- b_lasso_r(); be <- b_enet_r()
    lasso_sparse  <- any(abs(bl) < 0.01)
    enet_sparse   <- any(abs(be) < 0.01)
    msg <- if (lasso_sparse && input$show_lasso) {
      "Lasso has zeroed a coefficient — this is variable selection. The ellipse tangent hit a diamond corner."
    } else if (enet_sparse && input$show_enet) {
      "Elastic net zeroed a coefficient. Try α closer to 1 for more sparsity."
    } else if (input$show_lasso && !lasso_sparse) {
      "No zeros yet. Increase λ or move OLS closer to an axis to see Lasso select variables."
    } else {
      "Ridge never zeros coefficients — the circle has no corners. Lasso's diamond corners sit on the axes."
    }
    div(style="font-size:12px;color:#555555;line-height:1.6;", msg)
  })
}

shinyApp(ui, server)