xychoices <- function() {
    objects <- ls(envir = globalenv())
    if (length(objects) == 0)
        stop("No suitable object with spatial coordinates found", call. = FALSE)
    dataChoices <-
        objects[sapply(objects, function(x) {
            is.matrix(get(x)) | (substr(class(get(x))[1], 1, 7) == "Spatial")
        })]
    return(dataChoices)
}

ui <- shinyUI(fluidPage(
    theme = "style.css",
    # Application title
    titlePanel("Generate R code to create a spatial weighting matrix"),
    fluidRow(
        column(
            3,
            style = "background-color:#eaeaea;",
            # Sidebar with a slider input for number of bins
            h4("nb options"),
            selectInput("xy", "Sp object or coordinates:",
                        xychoices()),
            selectInput(
                "nb",
                "Graph type:",
                c(
                    "Delaunay" = 1,
                    "Gabriel" = 2,
                    "Relative" = 3,
                    "Minimum spanning tree" = 4,
                    "Distance" = 5,
                    "K-nearest" = 6
                )
            ),
            conditionalPanel("input.nb == 5",
                             fluidRow(
                                 column(6, numericInput("dmin", "dmin", 0)), column(6, numericInput("dmax", "dmax", NULL))
                             ), 
                             fluidRow(
                                 column(12, uiOutput("dthresh"))
                             ))
            
        ,
        conditionalPanel(
            "input.nb == 6",
            sliderInput(
                "knn",
                "K-nearest",
                min = 1,
                max = 5,
                value = 1
            )
        )    
        ),
        
        column(
            3,
            style = "background-color:#eaeaea;",
            h4("listw options"),
            selectInput(
                "style",
                "Standardization style:",
                c("W", "B", "C", "U", "minmax", "S")
            ),
            selectInput("glist",
                        "General weights:",
                        c("NULL", "1 - d/max(d)" = 1, "1 - (d/max(d))^a" = 2, "1 / d^a" = 3)),
            conditionalPanel(
                "input.glist > 1",
                sliderInput(
                    "a",
                    "Coefficient 'a':",
                    min = 0.1,
                    max = 3,
                    value = 1
                )
            )    
            
        ), 
        column(6, align = "center",
            h4("Spatial weights"),
            plotOutput("DistPlot", height = "200px")
            
        )
        
    ),
    fluidRow(
        column(6,
            h4("R code (copy & paste in the R  console):"), 
            tags$style(type = 'text/css', '#getcode {background-color: #999999; color: white;}'),
            verbatimTextOutput("getcode"), 
            br(),
            tags$style(type = 'text/css', '#summarybutton {white-space:nowrap;}'),
            checkboxInput("summarybutton", "Display summary", FALSE), 
            conditionalPanel("input.summarybutton == true", verbatimTextOutput("summary"))
            ),
        column(6, align = "center",
            h4("Connectivity"),
        plotOutput(
        "Plot",
        #height = "500px",
        dblclick = "plot_dblclick",
        brush = brushOpts(id = "plot_brush",
                          resetOnNew = TRUE)
    )
    
))))


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
    ranges <- reactiveValues(x = NULL, y = NULL)
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot_dblclick, {
        brush <- input$plot_brush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
            ranges$y <- c(brush$ymin, brush$ymax)
            
        } else {
            ranges$x <- NULL
            ranges$y <- NULL
        }
    })
    
    tmp <- reactiveValues(nb = NULL, lw = NULL, d1 = 0, d2 = NA)

    printRcode <- function(input){
        ## nb code
        mycode <- list()
        mycode[[1]] <- "library(adespatial);library(sp);library(spdep)"
        mycode[[2]] <- "nb <- chooseCN(coordinates("
        mycode[[2]] <- paste0(mycode[[2]], input$xy)
        mycode[[2]] <- paste0(mycode[[2]], "), type = ", input$nb)
        if (input$nb == 5)
            mycode[[2]] <- paste0(mycode[[2]], ", d1 = ", signif(tmp$d1, 5), ", d2 = ", signif(tmp$d2, 5))
        if (input$nb == 6)
            mycode[[2]] <- paste0(mycode[[2]], ", k = ", input$knn)
        mycode[[2]] <- paste0(mycode[[2]], ", plot.nb = FALSE)")

        ## listw code
        if (input$glist != "NULL") {
            mycode[[length(mycode) + 1]] <-  paste0("distnb <- nbdists(nb, ", input$xy, ")")

            if (input$glist == 1)
                mycode[[length(mycode) + 1]] <- paste0("fdist <- lapply(distnb, function(x) 1 - x/max(dist(", input$xy, ")))")
            if (input$glist == 2)
                mycode[[length(mycode) + 1]] <- paste0("fdist <- lapply(distnb, function(x) 1 - (x/max(dist(", input$xy, ")))^", input$a, ")")
            if (input$glist == 3)
                mycode[[length(mycode) + 1]] <- paste0("fdist <- lapply(distnb, function(x) 1 / x^", input$a, ")")
            }
        mycode[[length(mycode) + 1]] <-  paste0("lw <- nb2listw(nb, style = '", input$style, "'")
        if (input$glist != "NULL")
            mycode[[length(mycode)]] <- paste0(mycode[[length(mycode)]], ", glist = fdist")
        mycode[[length(mycode)]] <- paste0(mycode[[length(mycode)]], ", zero.policy = TRUE)")
        
        return(mycode)
    }  
    
    output$getcode <- renderPrint({
        mycode <- printRcode(input)
        cat(paste(mycode, collapse = "\n"))
        })
    
    output$summary <- renderPrint(summary(tmp$nb))
    
    output$DistPlot <- renderPlot({
    par(mar = c(5, 4, 0, 0))
        Spobj <- sp::coordinates(get(input$xy))
        dxy <- dist(Spobj)
        fdist <- NULL
        if (input$glist != "NULL") {
            distnb <- spdep::nbdists(tmp$nb, Spobj)
            
            if (input$glist == 1)
                fdist <- lapply(distnb, function(x) 1 - x/max(dxy))
            if (input$glist == 2)
                fdist <- lapply(distnb, function(x) 1 - (x/max(dxy))^input$a)
            if (input$glist == 3)
                fdist <- lapply(distnb, function(x) 1 / x^input$a)
        }
        tmp$lw <- spdep::nb2listw(tmp$nb, style = input$style, glist = fdist, zero.policy = TRUE)
        m1 <- as.matrix(dxy)
        m2 <- spdep::listw2mat(tmp$lw)

        plot(m1[!diag(ncol(m1))], m2[!diag(ncol(m2))], pch = 20, xlab = "distance", ylab = "spatial weights")
    })
    
    output$Plot <- renderPlot({
        Spobj <- sp::coordinates(get(input$xy))
        if (input$nb == 5) {
            dthresh <- adespatial::give.thresh(dist(Spobj))
            output$dthresh <- renderUI({tagList("give.thresh returns", tags$strong(signif(dthresh, 5)))})
            if (is.na(input$dmin))
                tmp$d1 <- 0
            else
                tmp$d1 <- input$dmin
            if (is.na(input$dmax)) {
                tmp$d2 <- dthresh
            }else
                tmp$d2 <- input$dmax
            tmp$nb <-
                adespatial::chooseCN(Spobj,
                                     type = input$nb,
                                     d1 = tmp$d1,
                                     d2 = tmp$d2, plot.nb = FALSE)
        } else if (input$nb == 6) {
            tmp$nb <-
                adespatial::chooseCN(Spobj, type = input$nb, k = input$knn, plot.nb = FALSE)
        } else {
            tmp$nb <- adespatial::chooseCN(Spobj, type = input$nb, plot.nb = FALSE)
        }
        par(mar = c(0, 4, 0, 0))
        spdep::plot.nb(
            tmp$nb,
            Spobj,
            xlim = ranges$x,
            ylim = ranges$y,
            pch = 20,
            cex = 2
        )
        
        text(Spobj, attr(tmp$nb, "region.id"), pos = 3)
        box()
    })
})

# Run the application
shinyApp(ui = ui, server = server)
