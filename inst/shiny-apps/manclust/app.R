library(plotly)
library(shiny)
library(shinyBS)
library(shinythemes)
library(data.table)
library(DT)
library(brapirv2)
library(shinyWidgets)

options(warn =-1)
options(shiny.maxRequestSize=300*1024^2)
max_brapi_snp_number <- 10000

if (is.null(options()$brapi.cons)) {
  brapisupport <-FALSE
  brapi_connections <- NULL
} else {
  brapisupport <-TRUE
  brapi_connections <- names(options("brapi.cons")$brapi.cons)

}



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

ui <- fluidPage(theme = shinytheme("flatly"),
                title = "snpclust",
                shinysky::busyIndicator(wait = 1000, text = NULL),
                tags$link(rel = "stylesheet", type = "text/css", href = "custom-div.css"),
  navbarPage(title = "snpclust", id = "tabsetId",
             tabPanel("Load data",value = "load",
                      #navlistPanel(widths = c(1,11),"From",
                      h4("Load data from file or from BrAPI endpoint:"),
                      switchInput("brapiorfile",   label = "Click to choose",
                                  value = ifelse(brapisupport,TRUE,FALSE),
                                  onLabel = "BrAPI",
                                  offLabel = "File",labelWidth = 130, onStatus = "success", offStatus = "info"
                      ),
                      bsCollapse(id="loadfrom", open="From BrAPI",
                        bsCollapsePanel(title = "From file", style="info",
                        #tabPanel("File",
                                 h3("Load data from file"),
                                 #fluidRow(
                                   h4("File format"),
                                   checkboxInput('lc', 'LightCycler 96 Format', FALSE),
                                   checkboxInput('intertek_guess', 'Intertek format (will guess number of lines to skip)', FALSE, width = '100%'),
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
                                   fileInput('file1', 'Choose file to upload'
                                             #accept = c(
                                             #  'text/csv',
                                             #  'text/comma-separated-values',
                                             #  'text/tab-separated-values',
                                             #  'text/plain',
                                             #  '.csv',
                                             #  '.tsv'
                                             #)
                                   ),
                                 tableOutput("df_data_out")
                        ),
                        bsCollapsePanel(title = "From BrAPI", style="success",
                        #tabPanel("BrAPI",
                                 h3("Load data from BrAPI endpoint"),
                                 selectizeInput("mainbrapiendpoint","BrAPI end point", choices = brapi_connections),
                                 passwordInput("mainbrapitoken","Token"),
                                 actionButton("connect_brapi","Connect"),
                                 htmlOutput("mainbrapiendpoint_connect_res"),
                                 selectizeInput("brapi_program", "Program", choices = NULL, selected = NULL,
                                                options = list(placeholder = 'Select a database',
                                                            onInitialize = I('function() { this.setValue(""); }')
                                                            )
                                                ),
                                 selectizeInput("brapi_study", "Study", choices = NULL, selected = NULL,
                                                options = list(
                                                  placeholder = 'Select a project',
                                                  onInitialize = I('function() { this.setValue(""); }'))
                                                ),
                                htmlOutput("retrieve_variants_res")
                        )
                      )),
             tabPanel("Retrieve Samples information",value = "samples",
                        h3(
                        div(style="display:inline-block;",img(src="ibp.png", width="30px"), style="left;"),
                        div("BMS connection")),
                        column(width = 3,
                               selectizeInput("bmsendpoint","BMS end point", choices = brapi_connections),
                               passwordInput("bmstoken","BMS token"),
                               actionButton("connect_bms","Connect"),
                               htmlOutput("connect_res"),
                               tags$hr(),
                               #selectizeInput("crop","Crop", choices = NA),
                               selectizeInput("program","Programme", choices = NULL),
                               #tags$hr(),
                               #selectizeInput("sample_list","Sample List", choices = NULL),
                               #checkboxInput('loop_over_progs', 'Search in all programs', FALSE),
                               actionButton("fetch_samples","Fetch samples information")
                               ),
                        column(9,dataTableOutput("samples_info",height = "600px"),
                               div(style="display: inline-block;vertical-align:top;",actionButton("special_samples","Toggle selected as special samples")),
                               div(style="display: inline-block;vertical-align:top;",actionButton("update_samples","Update samples information")),
                               htmlOutput("update_samp_res")),
#                        column(2),
#                        column(2),
#
                               #checkboxInput('header', 'Header', TRUE),
                               #numericInput(inputId = 'skip',label = 'Number of lines to skip',value = 0)),
                        tags$hr()
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
                          selectizeInput("SNP", label = "SNP", choices = "",
                                         options = list(placeholder = 'Select a SNP',
                                                        onInitialize = I('function() { this.setValue(""); }')
                                         )),
                          selectizeInput("Plate", label = "Plate", choices = "", multiple=TRUE,options = list(plugins= list('remove_button'))),
                          #selectInput("whichcall", label = "Show Call", choices = c("current","new"),selected = "new"),
                          # radioButtons('whichcall', 'Display Call',
                          #              c(Current='current',
                          #                New='new'),
                          #              'current'),
                          switchInput("whichcall2",   label = "Switch current/new call",
                                      value = TRUE,
                                      onLabel = "Current",
                                      offLabel = "New", labelWidth = 180, onStatus = "info", offStatus = "warning"),
                          actionButton(inputId = "copycall", label = "Copy current to new"),
                          actionButton(inputId = "resetnewcall", label = "Reset new call"),
                          checkboxInput(inputId = "tetar",label = "Use Theta/R",value = 0),
                          checkboxInput(inputId = "fixed_ratio",label = "Fixed axes",value = 0),
                          tags$hr(),
                          bsCollapse(id="adv_geno_seetings", open=NULL,
                                     bsCollapsePanel(title = "Advanced Alleles/genotypes settings", style="primary",
                                                     h4("Alleles"),
                                                     selectizeInput("snp_x_allele","X Allele", choices = c("A","C","G","T","-","X"), selected="X"),
                                                     selectizeInput("snp_y_allele","Y Allele", choices = c("A","C","G","T","-","Y"), selected="Y"),
                                                     h4("Genotypes"),
                                                     numericInput("ploidy","Ploidy", value = 2, min = 1,max = 5, step = 1),
                                                     selectizeInput("allele_sep", "Allele separator",choices=c(":","/","|"), selected = ":"))
                          ),
                          tags$hr(),
                          uiOutput("score_buttons"),
                          tags$br(),
                          #actionButton(inputId = "updateY", label = "Score as Allele Y", style="color: #fff; background-color: #dc143c; border-color: #2e6da4"),br(),br(),
                          #actionButton(inputId = "updateH", label = "Score as Heterozygous", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),#br(),br(),
                          #actionButton(inputId = "updateX", label = "Score as Allele X", style="color: #fff; background-color: #3cb371; border-color: #2e6da4"),br(),br(),
                          actionButton(inputId = "updateU", label = "Score as Missing", style="color: #fff; background-color: #ff7f50; border-color: #ff7f50"),#br(),br(),
                          actionButton(inputId = "updateN", label = "Score as NTC", style="color: #fff; background-color: #E54FFF; border-color: #E54FFF"),br(),br(),
                          tags$hr(),
                          uiOutput("exportData")

                          #downloadButton('downloadData', 'Download new file')

                        ),
                        mainPanel(
                          plotlyOutput("plot", width = 800, height = 600)
                        )
                      )
             ))
)

server <- function(input, output, session) {

  values <- reactiveValues(df_data = NULL,
                           newdf = NULL,
                           samplesdfd = NULL,
                           toplot=NULL,
                           xcall=NULL,
                           ycall=NULL,
                           hcall=NULL,
                           cols=NULL,
                           recols=NULL,
                           intk_snpinfos=NULL,
                           snpinfos=NULL,
                           alls=NULL,
                           genots=NULL)
  scorebts <- reactiveValues()
  scorebts$ui <- list()
  o <- reactiveVal(list())


  observeEvent(input$lc,{
    if (input$lc){
      updateRadioButtons(session,inputId = "sep",selected = '\t')
      updateRadioButtons(session,inputId = "quote",selected = '')
      updateRadioButtons(session,inputId = "dec",selected = '.')
      updateCheckboxInput(session,inputId = "header", value=T)
      updateNumericInput(session, inputId = "skip", value = 0)
    }
  })

#  "mainbrapiendpoint"
#  "mainbrapitoken"
  observeEvent(input$brapiorfile,{
    if (input$brapiorfile){
        updateCollapse(session, id="loadfrom", open="From BrAPI", close = "From file")
        hideTab(inputId = "tabsetId", target = "samples")
        hideTab(inputId = "tabsetId", target = "match")
        output$exportData <- renderUI({
          actionButton(inputId = "pushtobrapi",label = "Save data to BrAPI endpoint", icon = icon(name = "cloud-upload-alt"))
        })
        if (!brapisupport){
        showNotification("To use BrAPI end-points, a list of brapi connections needs to be defined with options(brapi.cons= list(Connection1= brapirv2::brapi_connect(...))) before running the app", type="error",closeButton = TRUE, duration = NULL)
        #updateSwitchInput(session = session, inputId = "brapiorfile", value = FALSE)
      }
    }else{
      updateCollapse(session, id="loadfrom", open="From file", close = "From BrAPI")
      output$exportData <- renderUI({
        downloadButton('downloadData', 'Download recoded file')
      })

      showTab(inputId = "tabsetId", target = "samples")
      showTab(inputId = "tabsetId", target = "match")
    }
  })
  observeEvent(input$loadfrom,{
    if(input$loadfrom=="From file"){
      updateSwitchInput(session = session, inputId = "brapiorfile", value = FALSE)
    } else {
      updateSwitchInput(session = session, inputId = "brapiorfile", value = TRUE)
    }
  })
  observeEvent(input$connect_brapi,{

    if (input$mainbrapiendpoint!=""){
      brapicon <<- options()$brapi.cons[[input$mainbrapiendpoint]]
      brapicon$token <<- input$mainbrapitoken
      brapidbs <<- tryCatch(brapirv2::brapi_get_programs(brapicon),
                         error=function(e) e)
      # For offline testing
      #progs <<- data.table(name="toto")
      if ("error"%in%class(brapidbs)){
        output$mainbrapiendpoint_connect_res = renderText({paste("<span style=\"color:red\">Connection failed</span>")})
      }else {
        if (input$mainbrapitoken==""){
          output$mainbrapiendpoint_connect_res = renderText({paste("<span style=\"color:green\">Connection succeeded (No token provided, listing public datasets)</span>")})
        }else{
          output$mainbrapiendpoint_connect_res = renderText({paste("<span style=\"color:green\">Connection succeeded</span>")})
        }
        updateSelectizeInput(session, "brapi_program",choices = brapidbs$programDbId)
      }
    }
  })

  observeEvent(input$brapi_program,{
    if (input$mainbrapiendpoint!=""  & input$brapi_program!=""){
    brapicon <<- options()$brapi.cons[[input$mainbrapiendpoint]]
    brapicon$token <<- input$mainbrapitoken
    brapi_studies <<- tryCatch(brapirv2::brapi_get_studies(brapicon, trialDbId=input$brapi_program),
                                  error=function(e) e)
    updateSelectizeInput(session, "brapi_study",choices = brapi_studies$studyName, selected = NULL)
    }
  })

  observeEvent(input$brapi_study,{
    if (input$mainbrapiendpoint!="" & input$brapi_program!="" & input$brapi_study!=""){
      brapicon <<- options()$brapi.cons[[input$mainbrapiendpoint]]
      brapicon$token <<- input$mainbrapitoken
      selected_studyDbId <<- brapi_studies[brapi_studies$studyName==input$brapi_study, "studyDbId"]
      brapi_variantsets <<- tryCatch(brapirv2::brapi_get_variantsets(brapicon, studyDbId =  htmltools::urlEncodePath(selected_studyDbId)), error=function(e) e)
      brapi_variantsetsIds <- unique(brapi_variantsets$variantSetDbId)
      brapi_variantsetsIds <- brapi_variantsetsIds[!is.na(brapi_variantsetsIds)]
      brapi_variants <<- do.call(rbind,
                                 lapply(brapi_variantsetsIds,
                                        function(a) tryCatch({
                                          brapirv2::brapi_get_variants(brapicon, variantSetDbId = htmltools::urlEncodePath(a), pageSize = max_brapi_snp_number)
                                          },error=function(e) e)
                                 )
      )
      values$snpinfos <- data.table(brapi_variants)[,.(SNPID=variantNames, AlleleX=referenceBases, AlleleY=alternateBases)]
      #updateSelectizeInput(session, inputId = "SNP", choices = data.frame(label=brapi_variants$variantNames, value=brapi_variants$variantDbId), server = T)
      updateSelectizeInput(session, inputId = "SNP", choices = brapi_variants$variantNames, server = T, selected = "")
      if (nrow(brapi_variants)==max_brapi_snp_number){
        output$retrieve_variants_res = renderUI(HTML(paste("Found more than ",max_brapi_snp_number," variants:", paste(brapi_variants$variantNames[1:10],collapse = ", "), "...</br>", "Keeping only the first ",max_brapi_snp_number," variants")))
      }else{
        output$retrieve_variants_res = renderText({paste("Found", nrow(brapi_variants), "variants:", paste(brapi_variants$variantNames[1:10],collapse = ", "), "...")})
      }
    }
  })
  observe( {
    inFile <- input$file1
    if (!is.null(inFile)){
      #df<-read.table(inFile$datapath, header = input$header,
      #               sep = input$sep, quote = input$quote, skip = input$skip, dec = input$dec, stringsAsFactors = F)
      if (input$intertek_guess){
        rawfile<-scan(inFile$datapath, what = "character", sep = "\n",blank.lines.skip= F, quiet = T)
        updateNumericInput(session, "skip", value =  grep("^Data$",rawfile))
        updateRadioButtons(session, "sep",selected = ",")
        snpinforow <- grep("^SNPs$",rawfile)
        snpinforow_end <- grep("^Scaling$",rawfile)-2
        values$intk_snpinfos <- fread(inFile$datapath, header = input$header, sep = ",", quote = input$quote, skip = snpinforow, nrows = snpinforow_end-snpinforow)
      } else {
        updateNumericInput(session, "skip", value =  0)
      }
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
  observeEvent(input$connect_bms,{
      bmscon <<- options()$brapi.cons[[input$bmsendpoint]]
      if (!is.null(input$bmstoken)){
        bmscon$token <<- input$bmstoken
        progs <<- tryCatch(setDT(brapirv2::brapi_get_programs(bmscon, commonCropName = bmscon$commoncropname)),
                           error=function(e) e)
        if ("error"%in%class(progs)){
          output$connect_res = renderText({paste("<span style=\"color:red\">Connection failed</span>")})
        }else {
          output$connect_res = renderText({paste("<span style=\"color:green\">Connection succeeded</span>")})
          updateSelectizeInput(session, "program",choices = progs$programName, selected = progs$programName[1])
        }
      }
    })
  observe({
    if (input$program!=""){
      #selprogUUID <- progs[programName==input$program,programDbId]
      #samplelists <<- bmsapi_Get_sample_list_search(bmscon, crop = bmscon$commoncropname, programUUID = selprogUUID)
  #    updateSelectizeInput(session, "sample_list",choices = samplelists$listName)
    }
  })
  observeEvent(input$fetch_samples,{
    if (input$program!=""){
      # if (input$loop_over_progs == TRUE){
      #   progs_search <- progs$programDbId
      # } else {
      #   progs_search <- progs[programName==input$program,programDbId]
      # }
      # samples <- data.table(NULL)
      # for (pgDbId in progs_search){
      #   samplelists <- bmsapi_Get_sample_list_search(bmscon, crop = bmscon$commoncropname, programUUID = pgDbId)
      #   for (l in 1:nrow(samplelists)){
      #     samp <- bmsapi_Get_sample_list_download(con = bmscon,
      #                                             crop = bmscon$commoncropname,
      #                                             programUUID = pgDbId,
      #                                             listId = samplelists[l, id],
      #                                             listName = samplelists[l, listName])
      #     samples <- rbind(samples,samp)
      #   }
      # }
        sidslookup <- unique(values$df_data$SubjectID)
        sidslookup <- sidslookup[sidslookup!=""]
        sampsrchid <- brapi_post_search_samples(con = bmscon, sampleDbIds = sidslookup)
        samps <- brapi_get_search_samples_searchResultsDbId(con = bmscon, searchResultsDbId = sampsrchid$searchResultsDbId)
        nbpages <- attr(samps, which = "pagination")$totalPages
        if (nbpages > 1){
          samps <- rbind(samps, do.call(rbind,
                                        lapply(1:(nbpages-1),
                                               function(p) brapi_get_search_samples_searchResultsDbId(con = bmscon, searchResultsDbId = sampsrchid$searchResultsDbId, page = p)
                                               )
                                        )
          )
        }
        setDT(samps)
        samples <- unique(samps[,.(sampleDbId,sampleName,germplasmDbId)])
    }
      dfd<-unique(data.table(values$df_data)[,.(SubjectID, Found=FALSE,Special=FALSE)])
      samplesdfd<-samples[dfd, on=c(sampleDbId="SubjectID")]
      samplesdfd[!is.na(germplasmDbId), Found:=TRUE]
      output$update_samp_res = renderText({paste("<span style=\"color:green\">",nrow(samplesdfd[!is.na(germplasmDbId)]),"samples found out of ",nrow(samplesdfd)," samples. Use the 'Found' column to identify missing samples</span>")})
      values$samplesdfd <- samplesdfd
  })

  observeEvent(input$special_samples,{
    samplesdfd<-copy(values$samplesdfd)
    samplesdfd[input$samples_info_rows_selected, Special:=!Special]
    values$samplesdfd<-samplesdfd
  })

  observeEvent(input$update_samples,{
    dfd <- data.table(values$df_data)
    dfd <- values$samplesdfd[,.(sampleDbId,sampleName,germplasmDbId, Special)][dfd, on=c(sampleDbId="SubjectID")]
    setnames(dfd,old = "sampleDbId",new = "SubjectID")
    dfd [!is.na(germplasmDbId), Sample_Plot_Label:=paste0(sampleName," - GUID:",germplasmDbId)]
    dfd [is.na(Sample_Plot_Label), Sample_Plot_Label:=SubjectID]
    dfd [is.na(Special), Special:=FALSE]
    dfd [,Special:=c("Standard","Special")[as.numeric(Special)+1]]
    values$df_data <- dfd[,c(colnames(values$df_data),"Sample_Plot_Label", "Special"), with=F]
    output$update_samp_res = renderText({paste("<span style=\"color:green\">Samples information updated: use the Sample_Plot_Label column as Sample Column at next step</span>")})
    updateSelectizeInput(session, inputId = "Icol", selected = "Sample_Plot_Label")
    updateNavbarPage(session, "tabsetId", selected = "match")
  })

  output$samples_info<-renderDataTable({
    datatable(
      {values$samplesdfd},
      filter = list(position='top', clear=F),
      escape = F,
      rownames = F,
      extension = c("Scroller"),
      selection = 'multiple',
      option = list(
        scrollX = T, scrollY = 450, scrollCollapse = F, scroller = T,
        dom = 'Brti'
      ))
  })
  observe({
    updateSelectizeInput(session, "Xcol",choices = colnames(values$df_data), selected = "X")
    updateSelectizeInput(session, "Ycol",choices = colnames(values$df_data), selected = "Y")
    updateSelectizeInput(session, "Ccol",choices = c("",colnames(values$df_data)), selected = "Call")
    updateSelectizeInput(session, "Pcol",choices = c("",colnames(values$df_data)), selected = "MasterPlate")
    updateSelectizeInput(session, "Scol",choices = c("",colnames(values$df_data)), selected = "SNPID")
    updateSelectizeInput(session, "Icol",choices = c("",colnames(values$df_data)), selected = sort(colnames(values$df_data)[colnames(values$df_data)%in%c("Sample_Plot_Label","SubjectID")])[1])
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
    if (!any(colnames(temp)=="snpclustId")){
      temp<-data.frame(temp, snpclustId = c(1:nrow(temp)), stringsAsFactors = F)
    }
    #temp$X.Fluor<-temp$X.Fluor-min(temp$X.Fluor)
    #temp$Y.Fluor<-temp$Y.Fluor-min(temp$Y.Fluor)
    if (any(temp$Call%in%c("?","NA", "Uncallable", "Negative"))){
      temp[temp$Call%in%c("?","NA", "Uncallable", "Negative"),]$Call <- NA
    }
    values$newdf<-temp
    updateSelectizeInput(session, "Plate",choices = sort(unique(temp$Plate)))
    updateSelectizeInput(session, "SNP",choices = sort(unique(temp$SNP)))
    updateNavbarPage(session, "tabsetId", selected = "clust")
    if (input$intertek_guess){
      values$snpinfos <- values$intk_snpinfos[,.(SNPID,AlleleX, AlleleY)]
    } else {
      values$snpinfos <- data.table(SNPID=unique(temp$SNP), AlleleX="X", AlleleY="Y")
    }

  })
  observeEvent(input$Plate,{
    temp<-values$newdf
    selSNP<-input$SNP
    if (!is.null(input$Plate)){
      updateSelectizeInput(session, "SNP",selected=selSNP , choices = unique(temp[temp$Plate%in%input$Plate,"SNP"]),label = paste("SNP (Plates:",paste(input$Plate, collapse = ","),")",sep=""))
      if (!is.null(input$Plate)){
        values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
      }
      if (input$SNP!=""){
        values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
      }
    } else{
      updateSelectizeInput(session, "SNP",selected=selSNP , choices = unique(temp[,"SNP"]),label = "SNP")
      #updateSelectizeInput(session, "Plate" , choices = unique(temp[,"Plate"]))
    }
  })

  observeEvent(input$SNP,{
    temp<-values$newdf
    selplate<-input$Plate
    if (input$SNP!=""){
      updateSelectizeInput(session, "Plate", selected=selplate, choices = sort(unique(temp[temp$SNP==input$SNP,"Plate"])),label = paste("Plate (SNP:",input$SNP,")",sep=""))
    } else{
      updateSelectizeInput(session, "Plate", selected=selplate, choices = sort(unique(temp[,"Plate"])),label = "Plate")
    }
    if (input$SNP!=""){
      if (input$brapiorfile){
        brapi_variantsets <<- tryCatch(brapirv2::brapi_get_variantsets(brapicon, studyDbId =  htmltools::urlEncodePath(selected_studyDbId)), error=function(e) e)
        brapi_variantsetsIds <- unique(brapi_variantsets$variantSetDbId)
        brapi_variantsetsIds <- brapi_variantsetsIds[!is.na(brapi_variantsetsIds)]
        brapi_calls <<- do.call(rbind,
                                   lapply(brapi_variantsetsIds,
                                          function(a) tryCatch({
                                            brapi_get_calls(brapicon, variantDbId = htmltools::urlEncodePath(setDT(brapi_variants)[variantNames==input$SNP,variantDbId]), variantSetDbId = htmltools::urlEncodePath(a), expandHomozygotes = TRUE, sepPhased = "/", sepUnphased = "/", unknownString = "NA")
                                          },error=function(e) e)
                                   )
        )
        setDT(brapi_calls)
        #brapi_calls[nchar(genotype.values)==1 & genotype.values%in%c("A","C","G","T","N","-"), genotype.values:=paste0(genotype.values,"/",genotype.values)]
        if (any(colnames(brapi_calls)=="genotypeMetadata.fieldAbbreviation")){
          brapi_calls <- brapi_calls[genotypeMetadata.fieldAbbreviation=="FI" & !is.na(genotypeMetadata.fieldValue)]
          if (nrow(brapi_calls)>0){
            brapi_calls <- cbind(brapi_calls, brapi_calls[, tstrsplit(genotypeMetadata.fieldValue, split=",", names = c("X.Fluor","Y.Fluor"))])
            brapi_calls[, X.Fluor:=as.numeric(X.Fluor)]
            brapi_calls[, Y.Fluor:=as.numeric(Y.Fluor)]
            brapi_calls[, snpclustId:=1:.N]
            brapi_calls[, NewCall:="Unknown"]
            brapi_calls[, Plate:="Unknown"]
            brapi_calls[, variantName:=input$SNP]
            brapi_calls[, callSetName:=callSetDbId]
            #brapi_calls[,genotypeValue:=unlist(genotypeValue)]
            setnames(brapi_calls,
                     old=c("variantName",
                           "genotypeValue",
                           "callSetDbId",
                           "callSetName"),
                     new=c("SNP",
                           "Call",
                           "SubjectID",
                           "SampName"))
          values$newdf <- data.frame(brapi_calls[,.(SNP, Call,snpclustId, SubjectID, SampName, X.Fluor, Y.Fluor, NewCall, Plate)])
          }
        }else{
          values$newdf <- data.frame(SNP="", Call="",snpclustId="", SubjectID="", SampName="", X.Fluor=0, Y.Fluor=0, NewCall="", Plate="")
        }
      }
      if (!is.null(values$newdf)){
        values$toplot<-values$newdf
        if (!is.null(input$Plate)){
          values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
        }
        if (input$SNP!=""){
          values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
        }
          calls <- unique(values$toplot$Call)
          calls <- calls[!is.na(calls)]
          values$cols <- calls
          names(values$cols) <- calls
          values$cols[gsub("[:/|]","",calls)==paste(rep(values$snpinfos[SNPID==input$SNP, AlleleX],2),collapse = "")] <- "#DC143CFF"
          values$cols[gsub("[:/|]","",calls)==paste(rep(values$snpinfos[SNPID==input$SNP, AlleleY],2),collapse = "")] <- "#3CB371FF"
          values$cols[gsub("[:/|]","",calls)==paste(c(values$snpinfos[SNPID==input$SNP, AlleleY],values$snpinfos[SNPID==input$SNP, AlleleX]),collapse = "")] <- "#00CCC5FF"
          values$cols[gsub("[:/|]","",calls)==paste(c(values$snpinfos[SNPID==input$SNP, AlleleX],values$snpinfos[SNPID==input$SNP, AlleleY]),collapse = "")] <- "#00CCC5FF"
          values$cols[gsub("[:/|]","",calls)=="NTC"] <- "#E54FFF"
          values$cols[!grepl("^#([A-Fa-f0-9]{8}|[A-Fa-f0-9]{6})$",values$cols)] <- "#FF7F50FF"
          guess_seps <- gsub("^.*([:/|]).*$","\\1",na.omit(grep("[[:punct:]]",calls, value = T)))
          guess_sep <- unique(guess_seps)[which.max(table(guess_seps))]

#$$$
          updateSelectizeInput(session = session, inputId = "allele_sep", selected = guess_sep)
          updateSelectizeInput(session = session, inputId = "snp_x_allele", selected = values$snpinfos[SNPID==input$SNP, AlleleX])
          updateSelectizeInput(session = session, inputId = "snp_y_allele", selected = values$snpinfos[SNPID==input$SNP, AlleleY])
      }
    }
    })

  observeEvent(input$snp_x_allele,{
    if (!is.null(values$snpinfos)){
      values$snpinfos <- copy(values$snpinfos[SNPID==input$SNP, AlleleX:=input$snp_x_allele])
    }
  })
  observeEvent(input$snp_y_allele,{
    if (!is.null(values$snpinfos)){
      values$snpinfos <- copy(values$snpinfos[SNPID==input$SNP, AlleleY:=input$snp_y_allele])
    }
  })

  observe({
    if (!is.null(values$snpinfos) & input$SNP!=""){
      #values$alls <- c(values$snpinfos[SNPID==input$SNP, AlleleX], values$snpinfos[SNPID==input$SNP, AlleleY])
      values$alls <- c(input$snp_x_allele, input$snp_y_allele)
      n <- input$ploidy
      #values$genots <- unlist(lapply(0:n, function(a) paste(c(rep(values$alls[1],a), rep(values$alls[2],n-a)), collapse = input$allele_sep)))
      values$genots <- unlist(lapply(n:0, function(a) paste(c(rep(values$alls[2],a), rep(values$alls[1],n-a)), collapse = "")))
      isolate({
        values$recols <- c("#3CB371FF", colorRampPalette(c("#00CCC5", "#4E00D6"))(n-1), "#DC143CFF")
        #names(values$recols) <- values$genots
        names(values$recols) <- unlist(lapply(strsplit(values$genots, split = ""), function(a) paste(a,collapse=input$allele_sep)))
        values$recols <- c(values$recols, NTC="#E54FFF")
      })
      scorebts$ui <- lapply(values$genots, function(g) list(actionButton(inputId = paste0("scoreb",gsub(input$allele_sep,"",g)),
                                                                         label = paste0("Score as ",paste(unlist(strsplit(g, split = "")), collapse = input$allele_sep)),
                                                                         style=paste0("color: #fff; background-color: ",values$recols[[which(values$genots==g)]],"; border-color: ",values$recols[[which(values$genots==g)]]))))
      output$score_buttons <- renderUI({scorebts$ui})
    }
  })

  observe({
    req(!is.null(values$genots))
    isolate({
      lapply(o(), function(x){x$destroy()})
      #inputBtn <- paste0("scoreb", gsub(input$allele_sep,"",values$genots))
      inputBtn <- paste0("scoreb", values$genots)
      #o(lapply(values$genots, function(x){
      o(lapply(inputBtn, function(x){
          #observeEvent(input[[paste0("scoreb", gsub(input$allele_sep,"",x))]],{
          observeEvent(input[[x]],{
          d <- event_data("plotly_selected")
          temp<-values$newdf
          g <- gsub("scoreb","",x)
          temp$NewCall[temp$snpclustId%in%d$key]<- paste(unlist(strsplit(g, split = "")),collapse=input$allele_sep)
          values$newdf<-temp
          values$toplot<-values$newdf
          if (!is.null(input$Plate)){
            values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
          }
          if (input$SNP!=""){
            values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
          }
        })
      }))
    })
  })

  observeEvent(input$copycall,{
    temp<-values$newdf
    if (!is.null(input$Plate)){
      temp[temp$Plate%in%input$Plate & temp$SNP==input$SNP,"NewCall"]<-gsub(" ", "_", temp[temp$Plate%in%input$Plate & temp$SNP==input$SNP,"Call"])
    }else{
      temp[temp$SNP==input$SNP,"NewCall"]<-gsub(" ", "_", temp[temp$SNP==input$SNP,"Call"])
    }
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }
  })
  observeEvent(input$resetnewcall,{
    temp<-values$newdf
    if (!is.null(input$Plate)){
      temp[temp$Plate%in%input$Plate & temp$SNP==input$SNP,"NewCall"]<-"Unknown"
    }else{
      temp[temp$SNP==input$SNP,"NewCall"]<-"Unknown"
    }
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })

  observeEvent(input$updateX,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$snpclustId%in%d$key]<- values$xcall #"Allele_X"
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })
  observeEvent(input$updateY,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$snpclustId%in%d$key]<- values$ycall #"Allele_Y"
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })
  observeEvent(input$updateH,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$snpclustId%in%d$key]<- values$hcall #"Both_Alleles"
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })
  observeEvent(input$updateU,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$snpclustId%in%d$key]<- "NA"
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })
  observeEvent(input$updateN,{
    d <- event_data("plotly_selected")
    temp<-values$newdf
    temp$NewCall[temp$snpclustId%in%d$key]<- "NTC"
    values$newdf<-temp
    values$toplot<-values$newdf
    if (!is.null(input$Plate)){
      values$toplot<-values$toplot[values$toplot$Plate%in%input$Plate,]
    }
    if (input$SNP!=""){
      values$toplot<-values$toplot[values$toplot$SNP==input$SNP,]
    }

  })

  output$df_data_out <- renderTable(head(values$df_data))

  observe({
    if (!is.null(values$newdf)){
      if (input$SNP!=""){
        ptitle<-paste(ifelse(input$SNP%in%c("","Any SNP"),"",input$SNP))#,ifelse(input$Plate%in%c("","Any SNP"),"",paste("-",paste(input$Plate,collapse = ","))))
        if (input$tetar == TRUE){
          isolate({
            if (any(!c("R","Theta")%in%colnames(values$toplot))){
              values$toplot<-cbind(values$toplot,xy2ThetaR(values$toplot[,c("X.Fluor","Y.Fluor")]))
            }
          })
          output$plot <- renderPlotly({
            # if (input$whichcall2==FALSE){
            #   cols <- c("Allele_X" = "#3CB371FF", "Allele_Y" = "#DC143CFF", "Both_Alleles" = "#337AB7FF", "NA" = "#FF7F50FF", "Negative"="#808080FF")
            #   names(cols) <- c(values$xcall,values$ycall,values$hcall,"NA", "Negative" )
            #   p <- ggplot(values$toplot,aes(x=Theta, y=R, colour=NewCall, key= snpclustId, text=paste("Sample:",SampName))) +  geom_point() #+facet_wrap(~Experiment_Name,ncol = 2)
            #   p <- p + scale_colour_manual(values = cols)
            # }else{
            #   p <- ggplot(values$toplot[values$toplot$Plate%in%input$Plate & values$toplot$SNP==input$SNP,],aes(x=Theta, y=R, colour=Call, key= snpclustId, text=paste("Sample:",SampName))) +  geom_point()+ aes(shape=Special) + scale_shape_manual(values =c(Standard=16,Special=11), name="")  #+facet_wrap(~Experiment_Name,ncol = 2)
            # }
            if (input$whichcall2==FALSE){
              #cols <- c("Allele_X" = "#3CB371FF", "Allele_Y" = "#DC143CFF", "Both_Alleles" = "#337AB7FF", "NA" = "#FF7F50FF", "Negative"="#808080FF")
              #names(cols) <- c(values$xcall,values$ycall,values$hcall,"NA", "Negative" )
              if(any(colnames(values$toplot)=="Special")){
                p <- ggplot(values$toplot,aes(x=Theta, y=R, colour=NewCall, key = snpclustId, text=paste(paste("Sample:",SampName),paste("Call:",Call), sep="\n"))) + geom_point()+ aes(shape=Special) + scale_shape_manual(values = c(Standard=16,Special=11), name="") #+facet_wrap(~Experiment_Name,ncol = 2)
              }else{
                p <- ggplot(values$toplot,aes(x=Theta, y=R, colour=NewCall, key = snpclustId, text=paste(paste("Sample:",SampName),paste("Call:",Call), sep="\n")))+ geom_point()
              }
              p <- p +  scale_colour_manual(values = values$recols, na.value = "#FF7F50FF")
            }else{
              if(any(colnames(values$toplot)=="Special")){
                p <- ggplot(values$toplot,aes(x=Theta, y=R, colour=Call, key = snpclustId, text=paste(paste("Sample:",SampName),paste("NewCall:",NewCall), sep="\n"))) +  geom_point()+ aes(shape=Special) + scale_shape_manual(values = c(Standard=16,Special=11), name="")  + coord_fixed(ratio = 1,xlim = c(0,maxfluo), ylim = c(0,maxfluo)) #+facet_wrap(~Experiment_Name,ncol = 2)
                p <- p + scale_colour_manual(values = values$cols, na.value = "#FF7F50FF")
              }else{
                p <- ggplot(values$toplot,aes(x=Theta, y=R, colour=Call, key = snpclustId, text=paste(paste("Sample:",SampName),paste("NewCall:",NewCall), sep="\n"))) +  geom_point()   #+facet_wrap(~Experiment_Name,ncol = 2)
                p <- p + scale_colour_manual(values = values$cols, na.value = "#FF7F50FF")
              }
            }

            ggplotly(p+ggtitle(ptitle)) %>% layout(dragmode = "lasso")
          })

        }else{
          #toplot<-values$newdf
          #if (!is.null(input$Plate)){
          #  toplot<-toplot[toplot$Plate%in%input$Plate,]
          #}
          #if (input$SNP!=""){
          #  toplot<-toplot[toplot$SNP==input$SNP,]
          #}
          values$toplot$X.Fluor<-values$toplot$X.Fluor-min(values$toplot$X.Fluor)
          values$toplot$Y.Fluor<-values$toplot$Y.Fluor-min(values$toplot$Y.Fluor)

          maxfluo<-max(c(values$toplot$X.Fluor,values$toplot$Y.Fluor))
          minfluo<-min(c(values$toplot$X.Fluor,values$toplot$Y.Fluor))
          output$plot <- renderPlotly({
            if (input$whichcall2==FALSE){

              #names(cols) <- c(values$xcall,values$ycall,values$hcall,"NA", "Negative" )
              if(any(colnames(values$toplot)=="Special")){
                p <- ggplot(values$toplot,aes(x=X.Fluor, y=Y.Fluor, colour=NewCall, key = snpclustId, text=paste("Sample:",SampName)))+ geom_point()+ aes(shape=Special) + scale_shape_manual(values = c(Standard=16,Special=11), name="") #+facet_wrap(~Experiment_Name,ncol = 2)
              }else{
                p <- ggplot(values$toplot,aes(x=X.Fluor, y=Y.Fluor, colour=NewCall, key = snpclustId, text=paste("Sample:",SampName)))+ geom_point()
              }
              if (input$fixed_ratio){
                p <- p + coord_fixed(ratio = 1, xlim = c(0,maxfluo), ylim = c(0,maxfluo))
              }
              p <- p +  scale_colour_manual(values = values$recols, na.value = "#FF7F50FF")
            }else{
              if(any(colnames(values$toplot)=="Special")){
                p <- ggplot(values$toplot,aes(x=X.Fluor, y=Y.Fluor, colour=Call, key = snpclustId, text=paste("Sample:",SampName))) +
                     geom_point()+ aes(shape=Special) + scale_shape_manual(values = c(Standard=16,Special=11), name="")
                if (input$fixed_ratio){
                  p <- p + coord_fixed(ratio = 1, xlim = c(0,maxfluo), ylim = c(0,maxfluo))
                }
                p <- p + scale_colour_manual(values = values$cols, na.value = "#FF7F50FF")
              }else{
                p <- ggplot(values$toplot,aes(x=X.Fluor, y=Y.Fluor, colour=Call, key = snpclustId, text=paste("Sample:",SampName))) +  geom_point()
                if (input$fixed_ratio){
                  p <- p + coord_fixed(ratio = 1, xlim = c(0,maxfluo), ylim = c(0,maxfluo))
                }
                p <- p + scale_colour_manual(values = values$cols, na.value = "#FF7F50FF")

              }
            }
            ggplotly(p+ggtitle(ptitle)) %>% layout(dragmode = "lasso")
          })
        }
      }else{
        output$plot <- renderPlotly({NULL})
      }
      }
  })
  observeEvent(input$pushtobrapi,{
    showModal(modalDialog(
      "Push of data back to BrAPI endpoint is not yet implemented",
      easyClose = TRUE
    ))
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
