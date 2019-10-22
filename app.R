library( shiny )

ui <- fluidPage(
    theme = "owntheme.css",
    
    tags$head( 
        tags$meta( name = "description", content = paste0( "ARIMA-folyamatok vizsgálatát lehetővé tevő, tulajdonságait ",
                                                           "szemléltető alkalmazás. ",
                                                           "Írta: Ferenci Tamás." ) ),
        tags$meta( property = "og:title", content = "ARIMA folyamatok vizsgálata" ),
        tags$meta( property = "og:type", content = "website" ),
        tags$meta( property = "og:locale", content = "hu_HU" ),
        tags$meta( property = "og:url",
                   content = "http://research.physcon.uni-obuda.hu/ARIMAlab/" ),
        tags$meta( property = "og:image",
                   content = "http://research.physcon.uni-obuda.hu/ARIMAlab_Pelda.png" ),
        tags$meta( property = "og:description", content = paste0( "ARIMA-folyamatok vizsgálatát lehetővé tevő, tulajdonságait ",
                                                                  "szemléltető alkalmazás. ",
                                                                  "Írta: Ferenci Tamás." ) ),
        tags$meta( name = "DC.Title", content = "ARIMA folyamatok vizsgálata" ),
        tags$meta( name = "DC.Creator", content = "Ferenci Tamás" ),
        tags$meta( name = "DC.Subject", content = "idősorelemzés" ),
        tags$meta( name = "DC.Description", content = paste0( "ARIMA-folyamatok vizsgálatát lehetővé tevő, tulajdonságait ",
                                                              "szemléltető alkalmazás. " ) ),
        tags$meta( name = "DC.Publisher",
                   content = "http://research.physcon.uni-obuda.hu/ARIMAlab/" ),
        tags$meta( name = "DC.Contributor", content = "Ferenci Tamás" ),
        tags$meta( name = "DC.Language", content = "hu_HU" )
    ),
    
    tags$div( id="fb-root" ),
    tags$script( HTML( "(function(d, s, id) {
                     var js, fjs = d.getElementsByTagName(s)[0];
                     if (d.getElementById(id)) return;
                     js = d.createElement(s); js.id = id;
                     js.src = 'https://connect.facebook.net/hu_HU/sdk.js#xfbml=1&version=v2.12';
                     fjs.parentNode.insertBefore(js, fjs);
                     }(document, 'script', 'facebook-jssdk'));" ) ),
    
    tags$style( ".shiny-file-input-progress {display: none}" ),
    
    titlePanel( "ARIMA folyamatok vizsgálata" ),
    
    p( "A program használatát részletesen bemutató súgó, valamint a technikai részletek",
       a( "itt", href = "https://github.com/tamas-ferenci/ARIMAlab",
          target = "_blank" ), "olvashatóak el." ),
    div( class="fb-like",
         "data-href"="http://research.physcon.uni-obuda.hu/ARIMAlab/",
         "data-layout"="standard", "data-action"="like", "data-size"="small",
         "data-show-faces"="true", "data-share"="true"), p(),
    
    sidebarLayout(
        sidebarPanel(
            selectInput( "p", "p", 0:5 ),
            selectInput( "d", "d", 0:2 ),
            selectInput( "q", "q", 0:5 ),
            actionButton( "resimulate", "Újra-szimulálás" ),
            downloadButton( "AbraLetoltesPDF", "Az ábra letöltése (PDF)" ),
            uiOutput( "pselect" ),
            uiOutput( "qselect" )
        ),
        
        mainPanel(
            plotOutput( "plot1", height = "600px" )
        )
    ),
    
    h4( "Írta: Ferenci Tamás (Óbudai Egyetem, Élettani Szabályozások Kutatóközpont), v1.00" ),
    
    tags$script( HTML( "var sc_project=11601191; 
                     var sc_invisible=1; 
                     var sc_security=\"5a06c22d\";
                     var scJsHost = ((\"https:\" == document.location.protocol) ?
                     \"https://secure.\" : \"http://www.\");
                     document.write(\"<sc\"+\"ript type='text/javascript' src='\" +
                     scJsHost+
                     \"statcounter.com/counter/counter.js'></\"+\"script>\");" ),
                 type = "text/javascript" )
)

server <- function(input, output) {
    
    output$pselect <- renderUI({
        if( input$p>0 ) lapply( 1:input$p, function( i )
            sliderInput( paste0( "phi", i ), paste0( "\U03C6", stringi::stri_unescape_unicode( paste0( "\\u", 2080+i ) ) ),
                         -2, 2, 0, 0.01 ) )
    })
    output$qselect <- renderUI({
        if( input$q>0 ) lapply( 1:input$q, function( i )
            sliderInput( paste0( "theta", i ), paste0( "\U03B8", stringi::stri_unescape_unicode( paste0( "\\u", 2080+i ) ) ),
                         -2, 2, 0, 0.01 ) )
    })
    
    xsim <- reactive({
        input$resimulate
        rnorm( 1000 + 1000 )
    })
    
    plotInput <- function() {
        p <- as.integer( input$p )
        q <- as.integer( input$q )
        d <- as.integer( input$d )
        
        phi <- c( input$phi1, input$phi2, input$phi3, input$phi4, input$phi5 )[ 1:p ]
        theta <- c( input$theta1, input$theta2, input$theta3, input$theta4, input$theta5 )[ 1:q ]
        
        x <- xsim()
        
        if( q>0 ) {
            x <- as.numeric( filter( x, c( 1, theta ), sides = 1L ) )
            x[ seq_along( theta ) ] <- 0
        }
        if( p>0&!is.null(phi) ) x <- as.numeric( filter( x, phi, method = "recursive" ) )
        x <- x[ -( seq_len( 1000 ) ) ]
        
        if( d>0 ) x <- diffinv( x, d )[ -seq_len( d ) ]
        
        x2 <- x[ !is.infinite( x )&abs(x)<1e290 ]
        p1 <- lattice::xyplot( x2~(1:length( x2 ) ), type = "l", main = "Egy szimulált trajektória", xlab = "", ylab = "" )
        
        maroots <- NULL
        arroots <- NULL
        maoutside <- TRUE
        aroutside <- TRUE
        if( !is.null( theta ) ) {
            maroots <- polyroot( c( 1, theta ) )
            maoutside <- abs( maroots )>1+1e-10
        }
        if( !is.null( phi )|d!=0 ) {
            arroots <- polyroot( polynom::as.polynomial( c( 1, if( is.null( phi ) ) NULL else -phi ) )*
                                     polynom::as.polynomial( c( 1, -1 ) )^d )
            aroutside <- abs( arroots )>1+1e-10
        }
        
        p2 <- lattice::xyplot( NA~NA, xlim = c( -2, 2 ), ylim = c( -2, 2 ), xlab = "Valós rész", ylab = "Képzetes rész",
                               scales = list( at = seq( -2, 2, 0.5 ), labels = seq( -2, 2, 0.5 ) ), aspect = "iso",
                               key = list( space = "right", text = list( c( "AR", "MA" ) ), points = list( pch = c( 19, 4 ) ) ),
                               panel = function( marootsin = maroots, maoutsidein = maoutside, arrootsin = arroots,
                                                 aroutsidein = aroutside ) {
                                   lattice::panel.abline( h = 0, v = 0 )
                                   angle.inc <- 2 * pi/100
                                   angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
                                   lattice::panel.polygon( cos(angles), sin(angles) )
                                   lattice::panel.points( marootsin[maoutsidein], pch = 4, col = "blue" )
                                   lattice::panel.points( marootsin[!maoutsidein], pch = 4, col = "red" )
                                   lattice::panel.points( arrootsin[aroutsidein], pch = 19, col = "blue" )
                                   lattice::panel.points( arrootsin[!aroutsidein], pch = 19, col = "red" )
                               } )
        
        if( any( !aroutside )|any( !maoutside ) ) {
            specres <- with( spectrum( x, span = 50, plot = FALSE ), data.frame( spec, freq ) )
            emp <- TRUE
        } else {
            specres <- with( TSA::ARMAspec( list( ar = if( all( phi==0 ) ) NULL else phi, ma = theta ), plot = FALSE ),
                             data.frame( spec, freq ) )
            emp <- FALSE
        }
        maxspec <- max( specres$spec )
        p3 <- lattice::xyplot( spec ~ freq, data = specres, ylim = c( 0, if( !is.nan( maxspec ) ) maxspec else 1 ), type = "l",
                               xlab = "", ylab = "",
                               main = if( emp )
                                   "Nem stacioner folyamat - (értelmetlen) empirikus spektrum\na szimulált trajektóriára" else
                                       "Elméleti spektrum" )
        
        if( any( !aroutside )|any( !maoutside ) ) {
            emp <- TRUE
            acfres <- with( acf( x, plot = FALSE ), data.frame( acf, lag ) )
        } else {
            emp <- FALSE
            if( (p==0&q==0)|(is.null(phi)&is.null(theta)) )
                acfres <- data.frame( lag = 0:30, acf = c( 1, rep( 0, 30 ) ) )
            else 
                acfres <- data.frame( lag = 0:30, acf = ARMAacf( ar = phi, ma = theta, lag.max = 30 ) )
        }
        p4 <- lattice::xyplot( acf ~ lag, data = acfres, type = "h", ylim = c( -1, 1 ), xlab = "", ylab = "",
                               abline = list( h = 0 ), main = if( emp )
                                   "Nem stacioner folyamat - (értelmetlen) empirikus ACF\na szimulált trajektóriára" else
                                       "Elméleti ACF")
        
        
        x2 <- c( rep( 0, 5 ), 1, rep( 0, 30 ) )
        if( q>0 ) {
            x2 <- as.numeric( filter( x2, c( 1, theta ), sides = 1L ) )
            x2[ seq_along( theta ) ] <- 0
        }
        if( p>0&!is.null(phi) ) x2 <- as.numeric( filter( x2, phi, method = "recursive" ) )
        x2 <- x2[ -( seq_len( 5 ) ) ]
        
        if( d>0 ) x2 <- diffinv( x2, d )[ -seq_len( d ) ]
        p5 <- lattice::xyplot( x2~(0:30), type="h", abline = list( h = 0 ), xlab = "", ylab = "",
                               main = "Elméleti impulzusválasz-függvény" )
        
        if( any( !aroutside )|any( !maoutside ) ) {
            emp <- TRUE
            pacfres <- with( pacf( x[ !is.infinite( x ) ], plot = FALSE, na.action = na.omit ), data.frame( pacf = acf, lag ) )
        } else {
            emp <- FALSE
            if( (p==0&q==0)|(is.null(phi)&is.null(theta)) )
                pacfres <- data.frame( lag = 1:30, pacf = rep( 0, 30 ) )
            else 
                pacfres <- data.frame( lag = 1:30, pacf = ARMAacf( ar = phi, ma = theta, lag.max = 30, pacf = TRUE ) )
        }
        p6 <- lattice::xyplot( pacf ~ lag, data = pacfres, type = "h", ylim = c( -1, 1 ), xlab = "", ylab = "",
                               abline = list( h = 0 ),
                               main = if( emp )
                                   "Nem stacioner folyamat - (értelmetlen) empirikus PACF\na szimulált trajektóriára" else
                                       "Elméleti PACF")
        
        gridExtra::grid.arrange( p1, p2, p3, p4, p5, p6, layout_matrix = matrix( c( 1, 2, 2, 5, 1, 3, 4, 6 ), nc = 2 ) )
        
        grid::grid.text( "Ferenci Tamás, 2019", 0, 0.02, gp = grid::gpar( fontface = "bold" ), just = "left" )
        grid::grid.text( "http://research.physcon.uni-obuda.hu", 1, 0.02, gp = grid::gpar( fontface = "bold" ), just = "right" )
    }
    
    output$plot1 <- renderPlot({
        print( plotInput() )
    })
    
    output$AbraLetoltesPDF <- downloadHandler(
        filename = "ARIMAlab.pdf",
        content = function( file ) {
            cairo_pdf( file, width = 9, height = 8 )
            print( plotInput() )
            dev.off( )
        } )
}

shinyApp(ui = ui, server = server)
