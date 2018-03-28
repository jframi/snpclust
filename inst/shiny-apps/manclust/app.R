library(plotly)
library(shiny)
library(shinythemes)
#library(snpclust)
library(data.table)

options(warn =-1)
valid_file<-function(df,lc){
  if (lc){
    if (all(c("Position","Sample Name","Genotype","Dye")%in%colnames(df))){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else{
    return(TRUE)
  }
}
ui <- fluidPage(theme = shinytheme("sandstone"),title = "snpclust",
  navbarPage(title = "snpclust", id = "tabsetId",
             tabPanel("Load File",value = "load",
                      fluidRow(
                      h3("File format"),
                      column(width = 2,
                             checkboxInput('lc', 'LightCycler Format', TRUE)),
                      column(width = 2,
                             radioButtons('sep', 'Separator',
                                          c(Comma=',',
                                            Semicolon=';',
                                            Tab='\t'),
                                          '\t')),
                      column(2,
                      radioButtons('quote', 'Quote',
                                   c(None='',
                                     'Double Quote'='"',
                                     'Single Quote'="'"),
                                   '')),
                      column(2,
                      radioButtons('dec', 'Decimal separator',
                                   c(Comma=',',
                                     'Point'='.'),
                                   '.')),
                      column(2,
                      checkboxInput('header', 'Header', TRUE),
                      numericInput(inputId = 'skip',label = 'Number of lines to skip',value = 0)),
                      tags$hr(),
                      fileInput('file1', 'Choose file to upload',
                                accept = c(
                                  'text/csv',
                                  'text/comma-separated-values',
                                  'text/tab-separated-values',
                                  'text/plain',
                                  '.csv',
                                  '.tsv'
                                )
                      )),
                      tableOutput("df_data_out")

             ),
             tabPanel("Match Columns",value = "match",
                      #selectInput("kcol", label = "Identification Column", choices = NA),
                      selectizeInput("Xcol", label = "X Fluo Column", choices = NA),
                      selectizeInput("Ycol", label = "Y Fluo Column", choices = NA),
                      selectizeInput("Ccol", label = "Call Column", choices = NA),
                      selectizeInput("Pcol", label = "Plate Column", choices = NA),
                      selectizeInput("Scol", label = "SNP Column", choices = NA),
                      selectizeInput("Icol", label = "Sample Column", choices = NA),
                      actionButton(inputId = "ok_matchcol", label = "OK")
             ),
             tabPanel("Clustering",value="clust",
                      sidebarLayout(
                        sidebarPanel(
                          selectInput("SNP", label = "SNP", choices = NA),
                          selectInput("Plate", label = "Plate", choices = NA),
                          #selectInput("whichcall", label = "Show Call", choices = c("current","new"),selected = "new"),
                          radioButtons('whichcall', 'Show Call',
                                       c(Current='current',
                                         New='new'),
                                       'new'),
                          actionButton(inputId = "copycall", label = "Copy current to new"),
                          actionButton(inputId = "resetnewcall", label = "Reset new call"),
                          checkboxInput(inputId = "tetar",label = "Use Theta/R",value = 0),
                          tags$hr(),
                          actionButton(inputId = "updateY", label = "Score as Allele Y", style="color: #fff; background-color: #dc143c; border-color: #2e6da4"),#br(),br(),
                          actionButton(inputId = "updateH", label = "Score as Heterozygous", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),br(),br(),
                          actionButton(inputId = "updateX", label = "Score as Allele X", style="color: #fff; background-color: #3cb371; border-color: #2e6da4"),br(),br(),
                          actionButton(inputId = "updateU", label = "Score as Unknown", style="color: #fff; background-color: #ff7f50; border-color: #2e6da4"),#br(),br(),
                          actionButton(inputId = "updateN", label = "Score as Negative", style="color: #fff; background-color: #808080; border-color: #2e6da4"),br(),br(),
                          tags$hr(),
                          downloadButton('downloadData', 'Download new file')

                        ),
                        mainPanel(
                          plotlyOutput("plot", width = 800, height = 600)
                        )
                      )
             ))
)

server <- function(input, output, session) {

  values <- reactiveValues(df_data = NULL, newdf = NULL)
  observeEvent(input$lc,{
    if (input$lc){
      updateRadioButtons(session,inputId = "sep",selected = '\t')
      updateRadioButtons(session,inputId = "quote",selected = '')
      updateRadioButtons(session,inputId = "dec",selected = '.')
      updateCheckboxInput(session,inputId = "header", value=T)
      updateNumericInput(session, inputId = "skip", value = 0)
    }
  })
  observe( {
    inFile <- input$file1
    if (!is.null(inFile)){
      #df<-read.table(inFile$datapath, header = input$header,
      #               sep = input$sep, quote = input$quote, skip = input$skip, dec = input$dec, stringsAsFactors = F)
      df<-fread(inFile$datapath, header = input$header,
                     sep = input$sep, quote = input$quote, skip = input$skip, dec = input$dec, stringsAsFactors = F)

      if(!valid_file(df,input$lc)){
        showModal(modalDialog("File doesn't look like a LightCycler file"))
      } else{
        if (input$lc){
          df<-df[EPF!="-"]
          df$EPF<-as.numeric(df$EPF)
          DF<-df
          dyes<-unique(DF$Dye)
          #browser()
          DF<-DF[,.(Dyes=paste(Dye,collapse = "/"),X.Fluor=EPF[1],Y.Fluor=EPF[2],Call=Genotype[1],`Sample Name`=`Sample Name`[1],Notes=Notes[1],`Sample Prep Notes`=`Sample Prep Notes`[1]),Position]
          DF[Call==paste0("Homozygote: ",dyes[1]),Call:="Allele_X"]
          DF[Call==paste0("Homozygote: ",dyes[2]),Call:="Allele_Y"]
          DF[Call=="Heterozygote",Call:="Both_Alleles"]
          DF[Call=="-",Call:="Unknown"]
          DF[,Experiment_Name:=input$file1$name]
          values$df_data <- data.frame(DF)
        }else{
          values$df_data <- data.frame(df) #data.table(df)
        }
      }
    }
  })
  observe({
    updateSelectizeInput(session, "Xcol",choices = colnames(values$df_data), selected = "X.Fluor")
    updateSelectizeInput(session, "Ycol",choices = colnames(values$df_data), selected = "Y.Fluor")
    updateSelectizeInput(session, "Ccol",choices = c("",colnames(values$df_data)), selected = "Call")
    updateSelectizeInput(session, "Pcol",choices = c("",colnames(values$df_data)), selected = "Experiment_Name")
    updateSelectizeInput(session, "Scol",choices = c("",colnames(values$df_data)), selected = "Sonde")
    updateSelectizeInput(session, "Icol",choices = c("",colnames(values$df_data)), selected = "Name")
    #updateSelectInput(session, "kcol",choices = colnames(values$df_data), selected = "order")
  })

  observeEvent(input$ok_matchcol,{
    temp<-values$df_data
    #browser()
    if (input$Ccol==""){
      temp<-data.frame(temp,Call="Unknown", stringsAsFactors = F)
    }else{
      colnames(temp)[match(input$Ccol,colnames(temp))]<-"Call"
    }
    if (input$Pcol==""){
      temp<-data.frame(temp,Plate="Any Plate", stringsAsFactors = F)
    }else{
      colnames(temp)[match(input$Pcol,colnames(temp))]<-"Plate"
    }
    if (input$Scol==""){
      temp<-data.frame(temp,SNP="Any SNP", stringsAsFactors = F)
    }else{
      colnames(temp)[match(input$Scol,colnames(temp))]<-"SNP"
    }
    if (input$Icol==""){
      temp<-data.frame(temp,SampName="", stringsAsFactors = F)
    }else{
      colnames(temp)[match(input$Icol,colnames(temp))]<-"SampName"
    }
    colnames(temp)[match(c(input$Xcol,input$Ycol),colnames(temp))]<-c("X.Fluor","Y.Fluor")
    if (!any(colnames(temp)=="NewCall")){
      temp<-data.frame(temp,NewCall="Unknown",  stringsAsFactors = F)
    }
    temp<-data.frame(temp, Id = c(1:nrow(temp)), stringsAsFactors = F)
    #temp$X.Fluor<-temp$X.Fluor-min(temp$X.Fluor)
    #temp$Y.Fluor<-temp$Y.Fluor-min(temp$Y.Fluor)
    values$newdf<-temp
    updateSelectInput(session, "Plate",choices = sort(unique(temp$Plate)))
    updateSelectInput(session, "SNP",choices = sort(unique(temp$SNP)))
    updateNavbarPage(session, "tabsetId", selected = "clust")

  })
  #observeEvent(input$Plate,{
  #  temp<-values$newdf
  #  updateSelectInput(session, "SNP",choices = unique(temp[temp$Plate==input$Plate,"SNP"]))
  #})
  observeEvent(input$SNP,{
    temp<-values$newdf
    updateSelectInput(session, "Plate",choices = sort(unique(temp[temp$SNP==input$SNP,"Plate"])))
  })

  observeEvent(input$copycall,{
    temp<-values$newdf
    temp[temp$Plate==input$Plate & temp$SNP==input$SNP,"NewCall"]<-gsub(" ", "_", temp[temp$Plate==input$Plate & temp$SNP==input$SNP,"Call"])
    values$newdf<-temp

  })
  observeEvent(input$resetnewcall,{
    temp<-values$newdf
    temp[temp$Plate==input$Plate & temp$SNP==input$SNP,"NewCall"]<-"Unknown"
    values$newdf<-temp

  })

  observeEvent(input$updateX,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$Id%in%d$key]<- "Allele_X"
    values$newdf<-temp
  })
  observeEvent(input$updateY,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$Id%in%d$key]<- "Allele_Y"
    values$newdf<-temp
  })
  observeEvent(input$updateH,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$Id%in%d$key]<- "Both_Alleles"
    values$newdf<-temp
  })
  observeEvent(input$updateU,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$Id%in%d$key]<- "Unknown"
    values$newdf<-temp
  })
  observeEvent(input$updateN,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$Id%in%d$key]<- "Negative"
    values$newdf<-temp
  })

  output$df_data_out <- renderTable(head(values$df_data))

  observe({
    if (!is.null(values$newdf)){
      ptitle<-paste(ifelse(input$SNP=="Any SNP","",input$SNP),ifelse(input$Plate=="Any Plate","",paste("-",input$Plate)))
      if (input$tetar == TRUE){
      toplot<-values$newdf
      toplot<-cbind(toplot,xy2ThetaR(toplot[,c("X.Fluor","Y.Fluor")]))
      output$plot <- renderPlotly({
        if (input$whichcall=="new"){
          cols <- c("Allele_X" = "#3CB371FF", "Allele_Y" = "#DC143CFF", "Both_Alleles" = "#337AB7FF", "Unknown" = "#FF7F50FF", "Negative"="#808080FF")
          p <- ggplot(toplot[toplot$Plate==input$Plate & toplot$SNP==input$SNP,],aes(x=Theta, y=R, colour=NewCall, key= Id, text=paste("Sample:",SampName))) +  geom_point() #+facet_wrap(~Experiment_Name,ncol = 2)
          p <- p + scale_colour_manual(values = cols)
        }else{
          p <- ggplot(toplot[toplot$Plate==input$Plate & toplot$SNP==input$SNP,],aes(x=Theta, y=R, colour=Call, key= Id, text=paste("Sample:",SampName))) +  geom_point() #+facet_wrap(~Experiment_Name,ncol = 2)
        }
        ggplotly(p+ggtitle(ptitle)) %>% layout(dragmode = "lasso")
      })

    }else{
      toplot<-values$newdf
      toplot$X.Fluor<-toplot$X.Fluor-min(toplot$X.Fluor)
      toplot$Y.Fluor<-toplot$Y.Fluor-min(toplot$Y.Fluor)
      toplot<-toplot[toplot$Plate==input$Plate & toplot$SNP==input$SNP,]

      maxfluo<-max(c(toplot$X.Fluor,toplot$Y.Fluor))
      output$plot <- renderPlotly({
        if (input$whichcall=="new"){
          cols <- c("Allele_X" = "#3CB371FF", "Allele_Y" = "#DC143CFF", "Both_Alleles" = "#337AB7FF", "Unknown" = "#FF7F50FF", "Negative"="#808080FF")
          p <- ggplot(toplot,aes(x=X.Fluor, y=Y.Fluor, colour=NewCall, key = Id, text=paste("Sample:",SampName))) +  geom_point() #+facet_wrap(~Experiment_Name,ncol = 2)
          p <- p + coord_fixed(ratio = 1, xlim = c(0,maxfluo), ylim = c(0,maxfluo))+ scale_colour_manual(values = cols)
        }else{
          p <- ggplot(toplot,aes(x=X.Fluor, y=Y.Fluor, colour=Call, key = Id, text=paste("Sample:",SampName))) +  geom_point() + coord_fixed(ratio = 1,xlim = c(0,maxfluo), ylim = c(0,maxfluo)) #+facet_wrap(~Experiment_Name,ncol = 2)
        }
        ggplotly(p+ggtitle(ptitle)) %>% layout(dragmode = "lasso")
      })
    }
    }
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$file1$name, '.recoded.txt', sep='')
    },
    content = function(file) {
      write.table(values$newdf, file, sep="\t", col.names = T, row.names = F)
    }
  )
}

shinyApp(ui, server)
