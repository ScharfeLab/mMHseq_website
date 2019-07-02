library(shiny)
library(DT)
all_amp_name=read.table("all_90_amplicon_name.txt",header=FALSE)
all_sample_name=read.table("all_sample_name.txt",header=TRUE)
fluidPage(
  #create UI for the application
  #input variable will be defined here and calculated in server() function
  #output plot returned from server() function should be displayed in main panel defined here.
  #Use multi-tab set up
  titlePanel("mMHseq Project-prototype"),
  tabsetPanel(
    tabPanel("MH individual figure",
             sidebarLayout(
               sidebarPanel(fluidRow(selectInput(inputId = "mh_region",
                                                 label = "MH region",
                                                 choices = as.character(all_amp_name$V1))),
                            fluidRow(selectInput(inputId="population",
                                                 label="Population",
                                                 choices=c("Sandawe",
                                                           "Chagga",
                                                           "EuroAmer",
                                                           "Biaka",
                                                           "TWChinese",
                                                           "Zaramo"))),#select population
                            fluidRow(htmlOutput("selectSample")),#select sample ID
                            fluidRow(tags$img(src='image/color_legend.png')),
                            width = 3
               ),
               mainPanel(#textOutput("txtOutput")
                 fluidRow(plotOutput("mh_figure_single",click = "plot_click")),
                 fluidRow(column(width=6,verbatimTextOutput("info")),
                          column(width=6,verbatimTextOutput("info2")))
               )
             )
    ),
    tabPanel("All MH figures by population",
             sidebarLayout(
               sidebarPanel(tags$style(".well {background-color:#FFF;}"),
                            fluidRow(selectInput(inputId = "mh_region3",
                                                 label = "MH region",
                                                 choices = as.character(all_amp_name$V1))),
                            fluidRow(selectInput(inputId="population2",
                                                 label="Population",
                                                 choices=c("Sandawe",
                                                           "Chagga",
                                                           "EuroAmer",
                                                           "Biaka",
                                                           "TWChinese",
                                                           "Zaramo"))),#select population
                            fluidRow(tags$img(src="image/color_legend.png")),
                            width = 3
               ),
               mainPanel(
                 fluidRow(uiOutput("mh_figure_all"))#depend on mh_figure_all variable
               )
             )
    ),
    tabPanel("MH SNP table",
             sidebarLayout(
               sidebarPanel(selectInput(inputId = "mh_region2",
                                        label = "MH region",
                                        choices = as.character(all_amp_name$V1)),
                            width = 3
               ),
               mainPanel(DT::dataTableOutput("mh_table"))
             )),
    tabPanel("MH haplotype frequency summary table",
             sidebarLayout(
               sidebarPanel(selectInput(inputId="mh_region4",
                                        label="MH region",
                                        choices=as.character(all_amp_name$V1)),
                            selectInput(inputId="population3",
                                        label="Population",
                                        choices=c("Sandawe",
                                                  "Chagga",
                                                  "EuroAmer",
                                                  "Biaka",
                                                  "TWChinese",
                                                  "Zaramo")),
                            width = 3
               ),
               mainPanel(DT::dataTableOutput("mh_frequency_table")))
    ),
    tabPanel("QC figure",
             sidebarLayout(#How the plot should be organized based on data
               sidebarPanel(selectInput(inputId = "dataset",
                                        label = "Dataset",
                                        choices = c("mMHOldPool","mMHPoolV2-1","mMHPoolV2-2")),
                            radioButtons(inputId = "qc_pdf",
                                         label = "QC",
                                         choices = c("SampleCoverage",
                                                     "AmpliconCoverageBySample",
                                                     "AmpliconCoverageByAmplicon",
                                                     "BasePairCoverage")),
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