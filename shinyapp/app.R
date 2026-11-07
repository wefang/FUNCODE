library(shiny)
library(stringr)
library(GenomicRanges)
library(dplyr)
library(shinyjs)
library(shinythemes)
library(openxlsx)
library(shinyWidgets)
library(dipsaus)

# change the path of loaded data
load('./shinyapp/data/regions.rda')
align1 = readRDS('./shinyapp/data/align1.rds')
align2 = readRDS('./shinyapp/data/align2.rds')
align3 = readRDS('./shinyapp/data/align3.rds')
align = bind_rows(list(align1, align2, align3))
rm(align1, align2, align3)
unalign = readRDS('./shinyapp/data/unalign.rds')
# GH-CO-B: significant pairs, appended so each row carries a GH-CO-V or a GH-CO-B score
ghcob = readRDS('./shinyapp/data/unalign_ghcob.rds')
unalign = bind_rows(unalign, ghcob)
unalign_human = c(unalign_human, GRanges(ghcob$human_chromosome, IRanges(ghcob$human_start, ghcob$human_end), index_id = ghcob$index_id))
unalign_mouse = c(unalign_mouse, GRanges(ghcob$mouse_chromosome, IRanges(ghcob$mouse_start, ghcob$mouse_end), index_id = ghcob$index_id))
rm(ghcob)


find_align <- function(gr, mapinfo){
  hits  = findOverlaps(gr, mapinfo)
  if(length(hits) == 0){
    return(as.data.frame(matrix(NA,0,0)))
  }else{
    overlaps = pintersect(gr[queryHits(hits)], mapinfo[subjectHits(hits)])
    df1 = data.frame(gr[queryHits(hits)])
    df2 = data.frame(mapinfo[subjectHits(hits)])
    
    overlap_df = cbind(df1[,c(1,2,3,6)], df2[,c(6)])
    colnames(overlap_df) = c('seqnames','start', 'end', 'input_index', 'index_id')
    overlap_df$overlap_width = width(overlaps)
    overlap_df$input_region = paste0(overlap_df$seqnames,':', overlap_df$start, '-', overlap_df$end)
    overlap_df = overlap_df[,c('input_index','input_region', 'overlap_width','index_id')]
    result = merge(overlap_df, align, by = 'index_id')
    result[,c(1,2,3,4)] = result[,c('input_index','input_region', 'overlap_width','index_id')]
    colnames(result)[1:4] = c('input_index','input_region', 'overlap_width','index_id')
    result = result[,-4]
    return(result)
  }
}

find_un_align <- function(gr, mapinfo){
  hits  = findOverlaps(gr, mapinfo)
  if(length(hits) == 0){
    return(as.data.frame(matrix(NA,0,0)))
  }else{
    overlaps = pintersect(gr[queryHits(hits)], mapinfo[subjectHits(hits)])
    df1 = data.frame(gr[queryHits(hits)])
    df2 = data.frame(mapinfo[subjectHits(hits)])
    
    overlap_df = cbind(df1[,c(1,2,3,6)], df2[,c(6)])
    colnames(overlap_df) = c('seqnames','start', 'end', 'input_index', 'index_id')
    overlap_df$overlap_width = width(overlaps)
    overlap_df$input_region = paste0(overlap_df$seqnames,':', overlap_df$start, '-', overlap_df$end)
    overlap_df = distinct(overlap_df[,c('input_index','input_region', 'overlap_width','index_id')])
    sel.unalign = unalign[unalign$index_id %in% overlap_df$index_id,]
    result = left_join(sel.unalign, overlap_df, by = 'index_id')
    result$hg38_region = paste0(result$human_chromosome,':', result$human_start,'-',result$human_end)
    result$mm10_region = paste0(result$mouse_chromosome,':', result$mouse_start,'-',result$mouse_end)
    result = result[,c('input_index', 'input_region','overlap_width','human_dhs_id','hg38_region','human_symbol', 'mouse_dhs_id','mm10_region','mouse_symbol', 
                       'modality','cov','cov_fdr','shared_motif','cob','cob_fdr')]
    return(result)
  }
}

ui <- fluidPage(
  useShinyjs(),
  navbarPage("FUNCODE", theme = shinytheme("sandstone"),
      tabPanel("Overview",
               #change content of overview tab
               strong(h2("FUNCODE: Functional Conservation of DNA Elements")),
               img(src = 'Picture1.png'), height="50%", width="50%"),
      tabPanel("Search",
               sidebarPanel(
                 h4(strong('FUNCODE scores characterize functional conservation of DNA elements')),
                 #change paper link
                 HTML("<p>3.3 millions pairs of human-mouse functionally conserved DNA elements were identified based on chromatin accessibility and histone modifications, including H3K4me3, H3K27ac and H3K4me1. For more see our <a href='https://www.google.com/'>paper</a>.</p>"),
                 tags$hr(),
                 #change example file link
                 fileInput("file", HTML("<strong>Choose BED File </strong><em><a href='https://www.google.com/'>example file</a></em>"), multiple = FALSE, placeholder = "No file selected"),
                 textAreaInput('text', "Or Enter Regions(chrom:start-end)", rows = 5, placeholder = 'chr1:11111-22222'),
                 tags$hr(),
                 radioButtons("spec", "Choose Species",
                              choices = c(`Human: GRCh38(hg38)` = "hg38",
                                          `Mouse: GRCm38(mm10)` = "mm10"),
                              selected = "hg38"),
                 actionButton('submit',"Search"),
                 actionButtonStyled('reset',"Reset", type="info"),
                 tags$hr(),
                 h5(strong("Save Results")),
                 downloadButton("downloadData", "Download as xlsx", style='width:200px;')),
               
               mainPanel(
                 h3('Sequence Homology Results'),
                 fluidRow(
                   column(1,selectInput("index1","Index:",'')),
                   column(3,selectInput("moda1","Modality:",choices = c("All","chromatin accessibility", "H3K27ac" ,"H3K4me1", "H3K4me3"),
                                        selected = 'All',
                                        multiple = FALSE)),
                   column(3,textInput('range.v','Set the range of COV:')),
                   column(3,textInput('range.v.fdr','Set the range of COV.FDR:'))
                 ),
                 fluidRow(
                   column(4, em("Be sure to choose the modality to show first before setting your ranges")),
                   column(3,textInput('range.b','Set the range of COB:')),
                   column(3,textInput('range.b.fdr','Set the range of COB.FDR:'))
                 ),
                 uiOutput("content1"),
                 tags$hr(),
                 h3('Gene Homology Results'),
                 fluidRow(
                   column(1,selectInput("index2","Index:",'')),
                   column(3,selectInput("moda2","Modality:",'')),
                   column(2, textInput('range2','Set the range of COV:')),
                   column(2, textInput('range2.cob','Set the range of COB:')),
                   column(2, textInput('range2.fdr','Set the range of FDR:'))),
                 uiOutput('content2'))
      ),
      #change content of FAQ
      tabPanel("FAQ", "Temporarily blank"),
      selected="Search"
    )
)


server <- function(input, output, session){
  observeEvent(input$reset,{
    updateTextInput(session, "text", value = "")
    reset("file")
  })
  observeEvent(input$submit,{
    if(input$text != ""){
      mylist =  strsplit(input$text, split = "\n")[[1]]
      index = chrom = start = end = c()
      for(i in 1:length(mylist)){
        index = c(index, i)
        chrom = c(chrom, strsplit(mylist[i], split = ':', fixed = T)[[1]][1])
        start = c(start, strsplit(strsplit(mylist[i], split = ':')[[1]][2], split = '-', fixed = T)[[1]][1])
        end = c(end, strsplit(mylist[i], split = '-',fixed = T)[[1]][2])
      }
      df = data.frame(seqnames = chrom, start = start, end = end, index = index)
    }
    if(!is.null(input$file)){
      df = read.table(input$file$datapath, sep = '\t')
      if(ncol(df) == 3){
        colnames(df) = c('seqnames', 'start',' end')
        df$index = 1:nrow(df)
      }else{
        colnames(df) = c('seqnames', 'start', 'end', 'index')
      }
    }
    if(input$text == "" & is.null(input$file)){
            sendSweetAlert(
        session = session,
        title = "Error!!",
        text = "No Input",
        type = "error"
      )
      return(NULL)
    }
      mygr = makeGRangesFromDataFrame(df, keep.extra.columns = T)
      if(input$spec == 'hg38'){
        align.res = find_align(mygr, align_human)
        unalign.res = find_un_align(mygr, unalign_human)
      }else if(input$spec == 'mm10'){
        align.res = find_align(mygr, align_mouse)
        unalign.res = find_un_align(mygr, unalign_mouse)
      }
      
      if(nrow(align.res) != 0){
        updateSelectInput(session, "index1",choices = c('All',as.character(unique(align.res$input_index))))
      }
      if(nrow(unalign.res) != 0){
        updateSelectInput(session, "index2",choices = c('All',as.character(unique(unalign.res$input_index))))
        updateSelectInput(session, "moda2",choices = c('All',as.character(unique(unalign.res$modality))))
      }
      
      output$content1 <- renderUI({
        if(nrow(align.res)==0){
          HTML("<span style='font-size: 20px; color: Red;'>No Regulatory Element Found</span>")
        }else{
          output$show1 <- DT::renderDataTable({
            data1 <- align.res
            if (input$index1 != "All") {
              data1 <- data1[data1$input_index == input$index1,]}
            if (input$moda1 != 'All'){
              sel.columns = c(paste0('cov_',gsub(' ', '_', input$moda1)),
                              paste0('cov_',gsub(' ', '_', input$moda1),'_fdr'),
                              paste0('cob_',gsub(' ', '_', input$moda1)),
                              paste0('cob_',gsub(' ', '_', input$moda1),'_fdr'))
              data1 <- align.res[,c('input_index', 'input_region','overlap_width','human_dhs_id','hg38_region','mm10_region', sel.columns)]
            }
            default.v.min = round(as.numeric(min(na.omit(data1[,7]))),3)
            default.v.max = round(as.numeric(max(na.omit(data1[,7]))),3)
            default.v.fdr.min = round(as.numeric(min(na.omit(data1[,8]))),3)
            default.v.fdr.max = round(as.numeric(max(na.omit(data1[,8]))),3)
            default.b.min = round(as.numeric(min(na.omit(data1[,9]))),3)
            default.b.max = round(as.numeric(max(na.omit(data1[,9]))),3)
            default.b.fdr.min = round(as.numeric(min(na.omit(data1[,10]))),3)
            default.b.fdr.max = round(as.numeric(max(na.omit(data1[,10]))),3)
            updateTextInput(session, 'range.v', placeholder = paste0(default.v.min,'-',default.v.max))
            updateTextInput(session, 'range.v.fdr', placeholder = paste0(default.v.fdr.min,'-',default.v.fdr.max))
            updateTextInput(session, 'range.b', placeholder = paste0(default.b.min,'-',default.b.max))
            updateTextInput(session, 'range.b.fdr', placeholder = paste0(default.b.fdr.min,'-',default.b.fdr.max))
            
            if (input$range.v != ''){
              min.v = as.numeric(strsplit(input$range.v, '-', fixed = TRUE)[[1]][1])
              max.v = as.numeric(strsplit(input$range.v, '-', fixed = TRUE)[[1]][2])
              data1 = data1[(data1[,7] >= min.v) &(data1[,7] <= max.v),]
            }
            if (input$range.v.fdr != ''){
              min.v.fdr = as.numeric(strsplit(input$range.v.fdr, '-', fixed = TRUE)[[1]][1])
              max.v.fdr = as.numeric(strsplit(input$range.v.fdr, '-', fixed = TRUE)[[1]][2])
              data1 = data1[(data1[,8] >= min.v.fdr) &(data1[,8] <= max.v.fdr),]
            }
            if (input$range.b != ''){
              min.b = as.numeric(strsplit(input$range.b, '-', fixed = TRUE)[[1]][1])
              max.b = as.numeric(strsplit(input$range.b, '-', fixed = TRUE)[[1]][2])
              data1 = data1[(data1[,9] >= min.b) &(data1[,9] <= max.b),]
            }
            if (input$range.b.fdr != ''){
              min.b.fdr = as.numeric(strsplit(input$range.b.fdr, '-', fixed = TRUE)[[1]][1])
              max.b.fdr = as.numeric(strsplit(input$range.b.fdr, '-', fixed = TRUE)[[1]][2])
              data1 = data1[(data1[,10] >= min.b.fdr) &(data1[,10] <= max.b.fdr),]
            }
            DT::datatable(data1, options = list(scrollX = TRUE))
          })
          DT::dataTableOutput("show1")
        }
      })
      
      output$content2 <- renderUI({
        if(nrow(unalign.res)==0){
          HTML("<span style='font-size: 20px; color: Red;'>No Regulatory Element Found</span>")
        }else{
          output$show2 <- DT::renderDataTable({
            data2 <- unalign.res
            if (input$index2 != "All") {
              data2 <- data2[data2$input_index == input$index2,]}
            if (input$moda2 != "All") {
              data2 <- data2[data2$modality == input$moda2,]}
            
            # a row has a GH-CO-V or a GH-CO-B score, so ranges skip the NA of the other one
            row.fdr = ifelse(is.na(data2$cov_fdr), data2$cob_fdr, data2$cov_fdr)
            span = function(x){
              x = na.omit(as.numeric(x))
              if(length(x) == 0) return('')
              paste0(round(min(x),3),'-',round(max(x),3))
            }
            updateTextInput(session, 'range2', placeholder = span(data2$cov))
            updateTextInput(session, 'range2.cob', placeholder = span(data2$cob))
            updateTextInput(session, 'range2.fdr', placeholder = span(row.fdr))

            in.range = function(x, range){
              lim = as.numeric(strsplit(range, '-', fixed = TRUE)[[1]])
              !is.na(x) & x >= lim[1] & x <= lim[2]
            }
            keep = rep(TRUE, nrow(data2))
            if (input$range2 != '') keep = keep & in.range(data2$cov, input$range2)
            if (input$range2.cob != '') keep = keep & in.range(data2$cob, input$range2.cob)
            if (input$range2.fdr != '') keep = keep & in.range(row.fdr, input$range2.fdr)
            data2 = data2[keep,]
            DT::datatable(data2, options = list(scrollX = TRUE))
          })
          DT::dataTableOutput("show2")
        }
      })
      
      wb = createWorkbook()
      sheet1 = addWorksheet(wb, "Sequence Homology")
      sheet2 = addWorksheet(wb, "Gene Homology")
      writeData(wb, align.res, sheet="Sequence Homology")
      writeData(wb, unalign.res, sheet="Gene Homology")
      
      output$downloadData <- downloadHandler(
        filename = function() {"result.xlsx"},
        content = function(file) {
          saveWorkbook(wb, file, overwrite = TRUE)
        }
      )

  })
}

shinyApp(ui, server)
