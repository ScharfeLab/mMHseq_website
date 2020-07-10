library(shiny)
library(DT)
library(markdown)

all_amp_name=read.table("all_90_amplicon_name.txt", header=FALSE)
all_sample_name=read.table("all_sample_name.txt", header=TRUE)
all_pop_name = sort(unique(all_sample_name$population))


fluidPage(
  #create UI for the application
  #input variable will be defined here and calculated in server() function
  #output plot returned from server() function should be displayed in main panel defined here.
  #Use multi-tab set up
  titlePanel("mMHseq Full", title = "Multiplex Microhaplotype Sequencing (mMHseq)"),
  tabsetPanel(
    tabPanel("MH individual figure",
             sidebarLayout(
               sidebarPanel(
                 fluidRow(column(12, selectInput(inputId = "mh_region",
                                                 label = "MH region",
                                                 choices = as.character(all_amp_name$V1)))),
                            fluidRow(column(12, selectInput(inputId="population",
                                                 label="Population",
                                                 choices= as.character(all_pop_name)))),#select population
                            fluidRow(column(12, htmlOutput("selectSample"))),#select sample ID
                 
                 # left side bar width = 2           
                 fluidRow(tags$img(style="height:80 px; width:90%", src='image/color_legend.png')),
                 width = 2
                 # left side bar width = 3
                 # fluidRow(tags$img(style="height:100 px; width:55%", src='image/color_legend.png')),
                 # width = 3
               ),
               mainPanel(#textOutput("txtOutput")
                 fluidRow(plotOutput("mh_figure_single",click = "plot_click")),
                 br(),
                 fluidRow(column(width= 6, verbatimTextOutput("info")),
                          column(width= 6, verbatimTextOutput("info2"))),
                 br(),
                 #br(),
                 fluidRow(#column(6, align = "right", uiOutput('ui.download.figure.simple')),
                          column(12, align = "center", uiOutput('ui.download.figure.complete'))),
                 #        column(4, align = "center", uiOutput('ui.download.click.info')))
                 br(),
                 fluidRow(column(12, align = "center", uiOutput('ui.download.click.info')))
               )
             )
    ),
    tabPanel("All MH figures by population",
             sidebarLayout(
               sidebarPanel(tags$style(".well {background-color:#FFF;}"),
                            #includeMarkdown("content/Population_tab.md"),
                            #hr(),
                            fluidRow(column(12, selectInput(inputId = "mh_region2",
                                                 label = "MH region",
                                                 #width = "95%", # will align to the left, so removed
                                                 
                                                 choices = as.character(all_amp_name$V1)))),
                            fluidRow(column(12, selectInput(inputId="population2",
                                                 label="Population",
                                                 #width = "95%",
                                                 
                                                 choices=as.character(all_pop_name)))),#select population
                            includeMarkdown("content/Population_tab.md"),
                            #br(),
                            
                            # left side bar width = 2           
                            fluidRow(tags$img(style="height:80 px; width:90%", src='image/color_legend.png')),
                            width = 2
                            # left side bar width = 3
                            # fluidRow(tags$img(style="height:100 px; width:55%", src='image/color_legend.png')),
                            # width = 3
               ),
               mainPanel(
                 fluidRow(uiOutput("mh_figure_all")) #depend on mh_figure_all variable
             ))
    ),
    tabPanel("MH SNP table",
             sidebarLayout(
               sidebarPanel(selectInput(inputId = "mh_region3",
                                        label = "MH region",
                                        choices = as.character(all_amp_name$V1)),
                            width = 2
               ),
               mainPanel(
                 fluidRow(DT::dataTableOutput("mh_table")),
                 fluidRow((column(12, align = "right", uiOutput('ui.download.mh.table'))))
                 )
             )),
    tabPanel("MH haplotype frequency summary table",
             sidebarLayout(
               sidebarPanel(selectInput(inputId="mh_region4",
                                        label="MH region",
                                        choices=as.character(all_amp_name$V1)),
                            selectInput(inputId="population4",
                                        label="Population",
                                        choices=as.character(all_pop_name)),
                            width = 2
               ),
               mainPanel(DT::dataTableOutput("mh_frequency_table")))
    ),
    tabPanel("QC figure",
             sidebarLayout(#How the plot should be organized based on data
               sidebarPanel(
                            #selectInput(inputId = "dataset",
                            #            label = "Dataset",
                            #            choices = c("mMHOldPool","mMHPoolV2-1","mMHPoolV2-2")),
                            #            choices = c("run1","run2","run3","run4")),
                            radioButtons(inputId = "qc_pdf",
                                         label = "Quality Control Type",
                                         choices = c(
                                                     "Sample coverage",
                                                     # "Amplicon coverage by amplicon",
                                                     "Amplicon coverage by sample",
                                                     "Base pair coverage"
                                                     )),
                            #radioButtons(inputId = "qc_num",
                            #             label = "QC_NUM",
                            #             choices = c("20X", "50X", "100X")),
                            width = 3
               ),
               mainPanel(uiOutput("qc_measure"))
             )
    )
  )
)
# Unused code
# radioButtons(inputId = "dataset_2",
#              label = "dataset",
#              choices = c("mMHOldPool",
#                          "mMHPoolV2-1",
#                          "mMHPoolV2-2"))
