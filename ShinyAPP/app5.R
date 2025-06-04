library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(MCMCpack)
library(progress)
library(coda)
library(ggplot2)
#library(tidyverse)
#library(rgdal)
n.iteration <- 6000
n.burnin <- 1000
n.thin <- 10
source("multinomial_MCMC.R")
load("thin_chain.RData")
load("World_ICdata.RData")
# Define UI ----
ui <- fluidPage(

  tags$head(
    tags$style(HTML("
      .input-title {
        background-color: #7EC0EE;
        padding: 5px;
        border-radius: 5px;
      }
      .down-title {
        background-color: #A4D2F3;
        padding: 5px;
        border-radius: 5px;
      }
      .blue-title {
        color: #2282b5;
      }
    "))
  ),

  fluidRow(
    column(8,
           tags$h1("Predicting population-level distribution of infant RSV hospitalisation at different months of age by birth month", class = "blue-title"),
           h5("Ling Guo, Nanjing Medical University (lingguo@stu.njmu.edu.cn) and You Li, Nanjing Medical University (you.li@njmu.edu.cn).")
    ),
    column(4,
           img(src="PROMISE.png",height= 50,width=160, style = "margin-top: 10px;"),
           img(src="NJMU.jpg",height= 80,width=80, style = "margin-top: 10px;"),
           img(src="TEAM.png",height= 80,width=80, style = "margin-top: 10px;")
    )
    
  ),

  div(h3("Introduction"),style="color:#437194"),
  p("This prediction tool is for understanding how infant RSV hospitalisation was distributed across different months of age 
    by month of birth. The prediction results can be used for identifying high-risk birth cohorts and high-risk chronological age windows 
    for prioritisation in RSV immunisation programmes to maximise their per-dose effectiveness. 
    Further details about the tool can be found in the publication (details to be updated). "),
  div(h3("Instruction for use"),style="color:#437194"),
  p("Please select the country for prediction and enter the number of RSV cases by each of the 12 calendar months in the region for prediction. 
    The entered number of RSV cases is ideally from a well-established surveillance or regionally representative 
    health-care providers (e.g., hospitals) with stable testing practice year round. 
    The application will return the predicted results in two forms. 
    The first shows the hospitalisation distribution (as proportion) by each month of age and each birth month. 
    The second shows the cumulative hospitalisation proportion across different months of age by birth month."),
  div("Pseudo data from Argentina are provided below for demonstration. 
      Users can reset and update the information below with their input data (may take a few seconds).",
      style = "color:green"),

  br(),

    h3("Input", class = "input-title"),
    selectInput("Country",h4("Select Country"),
                    choices = World_ICdata$region,
                    selected = "Argentina"),
   h4("Seasonality data"),
   div("Note: To improve the accuracy of the results, our model requires that the total number of RSV cases for the 12 months is not less than 100.",
            style = "color:red"),
   br(),
  flowLayout(
          numericInput("Jan", "January", value = 10, min = 0, step = 1),
          numericInput("Feb", "February", value = 8, min = 0, step = 1),
          numericInput("Mar", "March", value = 27, min = 0, step = 1),
          numericInput("Apr", "April", value = 110, min = 0, step = 1),
          numericInput("May", "May", value = 703, min = 0, step = 1),
          numericInput("Jun", "June", value = 1716, min = 0, step = 1),
          numericInput("Jul", "July", value = 1440, min = 0, step = 1),
          numericInput("Aug", "August", value = 620, min = 0, step = 1),
          numericInput("Sep", "September", value = 158, min = 0, step = 1),
          numericInput("Oct", "October", value = 61, min = 0, step = 1),
          numericInput("Nov", "November", value = 14, min = 0, step = 1),
          numericInput("Dec", "December", value = 8, min = 0, step = 1)
        ),
        actionButton("reset", "Reset"),
        actionButton("submit", "Update"),
  br(),
  br(),
        h3("Results", class = "input-title"),
        div(
          style = "display: flex; flex-direction: row;",
          column(width = 12,
                 radioButtons("selected_plots", "",
                              choices = list("Figure 1" = "Ribbonplot",
                                             "Figure 2" = "HeatMap",
                                             "Figure 3" = "HeatMap_cum"),
                              selected = "Ribbonplot",
                              inline = TRUE)
          )),
        br(),
        #style = "width: 1000px; height: 1200px;",
        plotOutput("BirthmonthPlot", width="900px", height = "600px"), #, width = "100%", height = "600px"), #输出图片自动填充
        textOutput("reminder"),
        br(),
        downloadButton("downloadData", "downloadData"),
  br(),
  br(),
  fluidRow(
    class = "down-title",
    column(4,
           br(),
           img(src="PROMISE_D.png",height= 60,width=300, style = "margin-left: 20px;")),
    column(8,
           p("Last update: May 20, 2025"),
           p("This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant agreement 101034339.
             This Joint Undertaking receives support from the European Union's Horizon 2020 research and innovation programme and EFPlA."))
  )
)

# Define server logic ----
server <- function(input, output,session) {
  BirthmonthSummary <- reactive({

    seasonality <- as.integer(c(input$Jan, input$Feb, input$Mar,
                                input$Apr, input$May, input$Jun,
                                input$Jul, input$Aug, input$Sep,
                                input$Oct, input$Nov, input$Dec))
    df_seasonality <- data.frame(region=input$Country,Jan=input$Jan,Feb =input$Feb,Mar =input$Mar,
                                 Apr=input$Apr, May=input$May, Jun=input$Jun,
                                 Jul =input$Jul, Aug=input$Aug, Sep=input$Sep,
                                 Oct =input$Oct, Nov=input$Nov,Dec =input$Dec)
    df_seasonality <- df_seasonality %>% left_join(World_ICdata,by="region")
    BirthmonthSummary <- genBirthRes(IC=df_seasonality$Income,Seasonality_each=seasonality)
    BirthmonthSummary$byBirthmonth.abb <-factor(month.abb[BirthmonthSummary$byBirthmonth.label], levels = month.abb)
    BirthmonthSummary <- BirthmonthSummary[,-1] %>% relocate(byBirthmonth.abb,.before = byAge.label) %>% rename(Birth_month=byBirthmonth.abb,
                                                                                                                Age=byAge.label,
                                                                                                                proportion.est=p_combined.est,
                                                                                                                proportion.lci=p_combined.lci,
                                                                                                                proportion.uci=p_combined.uci,
                                                                                                                cumulative_proportion.est=p_cumsum.est,
                                                                                                                cumulative_proportion.lci=p_cumsum.lci,
                                                                                                                cumulative_proportion.uci=p_cumsum.uci)
  })
  output$BirthmonthPlot <- renderPlot({
    selected_plots <- input$selected_plots
    if("Ribbonplot" %in% selected_plots){
      ggplot(data=BirthmonthSummary(), aes(x=Age)) + 
        geom_line(aes(y=proportion.est),colour="#428bca",linetype="dashed",linewidth=0.2)+
        geom_ribbon(alpha=0.3,aes(ymin=proportion.lci, ymax=proportion.uci),fill="#428bca")+
        geom_line(aes(y=cumulative_proportion.est),colour="#5cb85c",linetype="dashed",linewidth=0.2)+
        geom_ribbon(alpha=0.3,aes(ymin=cumulative_proportion.lci, ymax=cumulative_proportion.uci),fill="#5cb85c")+        
        labs(x="Age",y="Proportion",caption = "Blue indicated the proportion of RSV hospitalisations, and green indicated the cumulative proportion of RSV hospitalisations.\nLines indicated the medians, and shaded sections indicated 95% credible intervals.",
             title="Proportion and cumulative proportion of RSV hospitalisations by month of age in the first year of life \n for different birth months")+
        scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                           "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
        theme(legend.position = "bottom",
              text = element_text(size = 15),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11),
              axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
              axis.line.y = element_line(linetype=1,color="black",linewidth=0.25),
              plot.caption = element_text(hjust = 0,size=13))+
        facet_wrap(~ Birth_month) 
    } else{
      if("HeatMap" %in% selected_plots){
        ggplot(data = BirthmonthSummary(), aes(x=as.factor(Age), y=Birth_month, fill=proportion.est))+
          geom_tile()+
          scale_fill_gradient(low = "#f7fafc", high = "#005b96")+
          labs(x="Age",y="Birth Month",fill="Proportion",
               title="Proportion of RSV hospitalisations by birth month and by month of age")+
          scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                      "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
          theme(panel.background = element_blank(),
                text = element_text(size = 15),
                axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5, size = 12)) 
      }else{
        if("HeatMap_cum" %in% selected_plots){
          ggplot(data = BirthmonthSummary(), aes(x=as.factor(Age), y=Birth_month, fill=cumulative_proportion.est))+
            geom_tile()+
            scale_fill_gradient(low = "#fffdfd", high = "#c0281d")+
            labs(x="Age",y="Birth Month",fill="Cumulative Proportion",
                 title = "Cumulative proportion of RSV hospitalisations by birth month and by month of age")+
            scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                        "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
            theme(panel.background = element_blank(),
                  text = element_text(size = 15),
                  axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5, size = 12))
        } 
      }
    }
  })

  observeEvent(input$submit, {
    seasonality <- as.integer(c(input$Jan, input$Feb, input$Mar,
                                input$Apr, input$May, input$Jun,
                                input$Jul, input$Aug, input$Sep,
                                input$Oct, input$Nov, input$Dec))

    if (any(is.na(seasonality))) {
      output$reminder <- renderText({
        "Please enter full details"
      })
      output$BirthmonthPlot <- renderPlot(NULL)
      output$downloadData <- downloadHandler(NULL)
    } else {

      if(sum(seasonality) < 100) {
        showModal(modalDialog(
          title = "Warning",
          "Total number of RSV cases less than 100.",
          easyClose = TRUE
        ))
        output$BirthmonthPlot <- renderPlot(NULL)
        updateNumericInput(session, "Jan", value = NA)
        updateNumericInput(session, "Feb", value = NA)
        updateNumericInput(session, "Mar", value = NA)
        updateNumericInput(session, "Apr", value = NA)
        updateNumericInput(session, "May", value = NA)
        updateNumericInput(session, "Jun", value = NA)
        updateNumericInput(session, "Jul", value = NA)
        updateNumericInput(session, "Aug", value = NA)
        updateNumericInput(session, "Sep", value = NA)
        updateNumericInput(session, "Oct", value = NA)
        updateNumericInput(session, "Nov", value = NA)
        updateNumericInput(session, "Dec", value = NA)
        
      } else {
        output$reminder <- renderText(NULL)
        df_seasonality <- data.frame(region=input$Country,Jan=input$Jan,Feb =input$Feb,Mar =input$Mar,
                                     Apr=input$Apr, May=input$May, Jun=input$Jun,
                                     Jul =input$Jul, Aug=input$Aug, Sep=input$Sep,
                                     Oct =input$Oct, Nov=input$Nov,Dec =input$Dec)
        df_seasonality <- df_seasonality %>% left_join(World_ICdata,by="region")
        BirthmonthSummary <- genBirthRes(IC=df_seasonality$Income,Seasonality_each=seasonality)
        BirthmonthSummary$byBirthmonth.abb <-factor(month.abb[BirthmonthSummary$byBirthmonth.label], levels = month.abb)
        BirthmonthSummary <- BirthmonthSummary[,-1] %>% relocate(byBirthmonth.abb,.before = byAge.label) %>% rename(Birth_month=byBirthmonth.abb,
                                                                                                                    Age=byAge.label,
                                                                                                                    proportion.est=p_combined.est,
                                                                                                                    proportion.lci=p_combined.lci,
                                                                                                                    proportion.uci=p_combined.uci,
                                                                                                                    cumulative_proportion.est=p_cumsum.est,
                                                                                                                    cumulative_proportion.lci=p_cumsum.lci,
                                                                                                                    cumulative_proportion.uci=p_cumsum.uci)
        
        output$BirthmonthPlot <- renderPlot({
          selected_plots <- input$selected_plots
          if("Ribbonplot" %in% selected_plots){
            ggplot(data=BirthmonthSummary, aes(x=Age)) + 
              geom_line(aes(y=proportion.est),colour="#428bca",linetype="dashed",linewidth=0.2)+
              geom_ribbon(alpha=0.3,aes(ymin=proportion.lci, ymax=proportion.uci),fill="#428bca")+
              geom_line(aes(y=cumulative_proportion.est),colour="#5cb85c",linetype="dashed",linewidth=0.2)+
              geom_ribbon(alpha=0.3,aes(ymin=cumulative_proportion.lci, ymax=cumulative_proportion.uci),fill="#5cb85c")+        
              labs(x="Age",y="Proportion",caption = "Blue indicated the proportion of RSV hospitalisations, and green indicated the cumulative proportion of RSV hospitalisations.\nLines indicated the medians, and shaded sections indicated 95% credible intervals.",
                   title="Proportion and cumulative proportion of RSV hospitalisations by month of age in the first year of life \n for different birth months")+
              scale_x_continuous(breaks = seq(1,12,1),labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                                                 "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
              theme(legend.position = "bottom",
                    text = element_text(size = 15),
                    panel.background = element_blank(),
                    axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5,size = 11),
                    axis.line.x = element_line(linetype=1,color="black",linewidth=0.25),
                    axis.line.y = element_line(linetype=1,color="black",linewidth=0.25),
                    plot.caption = element_text(hjust = 0,size=13))+
              facet_wrap(~ Birth_month) 
          } else{
            if("HeatMap" %in% selected_plots){
              ggplot(data = BirthmonthSummary, aes(x=as.factor(Age), y=Birth_month, fill=proportion.est))+
                geom_tile()+
                scale_fill_gradient(low = "#f7fafc", high = "#005b96")+
                labs(x="Age",y="Birth Month",fill="Proportion",
                     title="Proportion of RSV hospitalisations by birth month and by month of age")+
                scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                            "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
                theme(panel.background = element_blank(),
                      text = element_text(size = 15),
                      axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5, size = 12)) 
            }else{
              if("HeatMap_cum" %in% selected_plots){
                ggplot(data = BirthmonthSummary, aes(x=as.factor(Age), y=Birth_month, fill=cumulative_proportion.est))+
                  geom_tile()+
                  scale_fill_gradient(low = "#fffdfd", high = "#c0281d")+
                  labs(x="Age",y="Birth Month",fill="Cumulative Proportion",
                       title = "Cumulative proportion of RSV hospitalisations by birth month and by month of age")+
                  scale_x_discrete(labels = c("0-<1m","1-<2m","2-<3m","3-<4m","4-<5m","5-<6m",
                                              "6-<7m","7-<8m","8-<9m","9-<10m","10-<11m","11-<12m"))+
                  theme(panel.background = element_blank(),
                        text = element_text(size = 15),
                        axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5, size = 12))
              } 
            }
          }
        })
        
        output$downloadData <- downloadHandler(
          filename = function() {
            paste("BirthmonthSummary", Sys.Date(), ".csv", sep = "")
          },
          content = function(file) {
            write.csv(BirthmonthSummary, file, row.names = FALSE)
          })
        #   output$BirthmonthTable <- renderTable({
        #   BirthmonthSummary()
        #   })
      }
    }  
    
  })
  

  observeEvent(input$reset, {
    
    output$reminder <- renderText({
      "Please enter full details"
    })
    output$BirthmonthPlot <- renderPlot(NULL)
    
    updateNumericInput(session, "Jan", value = NA)
    updateNumericInput(session, "Feb", value = NA)
    updateNumericInput(session, "Mar", value = NA)
    updateNumericInput(session, "Apr", value = NA)
    updateNumericInput(session, "May", value = NA)
    updateNumericInput(session, "Jun", value = NA)
    updateNumericInput(session, "Jul", value = NA)
    updateNumericInput(session, "Aug", value = NA)
    updateNumericInput(session, "Sep", value = NA)
    updateNumericInput(session, "Oct", value = NA)
    updateNumericInput(session, "Nov", value = NA)
    updateNumericInput(session, "Dec", value = NA)
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)
