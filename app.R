# R包加载----
library(shiny)
library(shinythemes)
library(survival)
library(survminer)
library(DT)
library(reshape2)
library(ggsci)
library(ggplot2)


# 选项设置
options(shiny.maxRequestSize = 500*1024^2)  # 500MB

# 自定义CSS样式----
custom_css <- "
#home_text {
  font-size: 18px;
  font-family: Arial, sans-serif;
  text-align: center;
  margin-top: 20px;
}
#home_logo {
  display: flex;
  justify-content: center;
  margin-top: 30px;
}
#plot_container {
  border: 2px solid #ccc; /* 边框颜色 */
  padding: 15px;         /* 内边距 */
  margin-bottom: 20px;   /* 底部外边距 */
  border-radius: 5px;    /* 圆角边框 */
  background-color: #f9f9f9; /* 背景颜色 */
}
#divider {
  border-top: 1px solid #ccc; /* 浅灰色分割线 */
  margin-top: 20px;
  margin-bottom: 20px;
}
"

# UI 部分----
ui <- navbarPage(
  title = "LUAD-Kla.Sig",
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML(custom_css))),  # 添加自定义CSS样式
  
  # Model页
  tabPanel(
    "Model",  
    
    sidebarPanel(
      width = 4,  
      tabsetPanel(
        id = "sidebar_tabs",
        
        h4("Data Selection"),
        selectInput("data_source", "Data Source:", choices = c("Use Provided Datasets", "Upload Your File")),
        
        conditionalPanel(
          condition = "input.data_source == 'Use Provided Datasets'",
          selectInput("dataset", "Choose Dataset:", choices = c("TCGA-LUAD", "GSE14814", "GSE29016", "GSE30219", "GSE31210",
                                                                "GSE37745", "GSE42127", "GSE50081", "GSE68465", "GSE72094"))
        ),
        
        conditionalPanel(
          condition = "input.data_source == 'Upload Your File'",
          fileInput("upload_file", "Upload CSV File:", accept = c(".csv")),
          helpText("Upload a CSV file containing 'Sample_ID', 'status', and 'time' in top 3 columns.")
        ),
        
        actionButton("run_fit", "Run Kla.Sig Model", class = "btn-primary"),
        
        tags$hr(style = "border-top: 2px solid #ccc;"),
        
        # h4("Table Settings"),
        # radioButtons("table_format", "Table Format:",
        #              choices = c(".csv", ".txt", ".xlsx"),
        #              selected = ".csv"),
        # downloadButton("download_table", "Download Table", class = "btn-primary")
      )
    ),
    
    mainPanel(
      width = 8,  
      
      # 标签页部分
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel(
          "K-M Plot", 
          fluidRow(
            column(12, 
                   plotOutput("km_plot", height = "500px", width = "500px")  # K-M plot
            )
          ),
          # fluidRow(
          #   column(12, 
          #          selectInput("km_plot_format", "Choose Image Format:", choices = c("png", "jpg", "jpeg", "pdf")),
          #          downloadButton("download_km_plot", "Download K-M Plot", class = "btn-primary")  # 下载K-M plot按钮
          #   )
          # )
        ),
        
        tabPanel(
          "ROC Curve", 
          fluidRow(
            column(12, 
                   plotOutput("roc_plot", height = "500px", width = "500px")  # ROC plot
            )
          ),
          # fluidRow(
          #   column(12, 
          #          selectInput("roc_plot_format", "Choose Image Format:", choices = c("png", "jpg", "jpeg", "pdf")),
          #          downloadButton("download_roc_plot", "Download ROC Curve", class = "btn-primary")  # 下载ROC plot按钮
          #   )
          # )
        )
      ),
      
      # 数据表格部分
      fluidRow(
        column(12, 
               h4("Dataset Table"), 
               DTOutput("dataset_table")  # 显示数据表格
        )
      )
    )
  )
)


# Server 部分----
server <- function(input, output, session) {
  
  # 数据处理：根据用户上传或选择的数据集进行加载
  dataset_input <- reactive({
    req(input$data_source)
    
    if (input$data_source == "Use Provided Datasets") {
      # 在这里加载选中的数据集
      if (input$dataset == "TCGA-LUAD") {
        # 示例：加载TCGA-LUAD数据集
        data <- readRDS("data/TCGA_LUAD.rds")  # 更改为您的实际路径
      } else if (input$dataset == "GSE14814") {
        data <- readRDS("data/GSE14814.rds")
      } else if (input$dataset == "GSE29016") {
        data <- readRDS("data/GSE29016.rds")
      } else if (input$dataset == "GSE30219") {
        data <- readRDS("data/GSE30219.rds")
      } else if (input$dataset == "GSE31210") {
        data <- readRDS("data/GSE31210.rds")
      } else if (input$dataset == "GSE37745") {
        data <- readRDS("data/GSE37745.rds")
      } else if (input$dataset == "GSE42127") {
        data <- readRDS("data/GSE42127.rds")
      } else if (input$dataset == "GSE50081") {
        data <- readRDS("data/GSE50081.rds")
      } else if (input$dataset == "GSE68465") {
        data <- readRDS("data/GSE68465.rds")
      } else if (input$dataset == "GSE72094") {
        data <- readRDS("data/GSE72094.rds")
      } 
      
      data2 <- scale(data[,-c(1:3)])
      data <- cbind(data[,c(1:3)], data2)
      return(data)
    } else if (input$data_source == "Upload Your File") {
      # 读取用户上传的文件
      req(input$upload_file)
      data <- read.csv(input$upload_file$datapath, check.names = F)
      data2 <- scale(data[,-c(1:3)])
      data <- cbind(data[,c(1:3)], data2)
      return(data)
    }
  })
  
  # 计算risk score
  rs_data <- reactive({
    req(dataset_input())
    data <- dataset_input()
    
    # 读取 Kla.Sig 预训练模型
    fit <- readRDS("data/KlaSigFit.rds")  # 确保文件路径正确
    pre_var <- fit[["xvar.names"]]  # 获取预测变量
    
    # 确保数据集中包含所有必要的列
    req(all(c("Sample_ID", "status", "time") %in% colnames(data)))
    
    # 选择相关变量
    if (all(pre_var %in% colnames(data))) {
      data3 <- data[, c("Sample_ID", "status", "time", pre_var)]
    } else {
      stop("Error: Missing predictor variables in dataset.")
    }
    
    # 计算风险评分
    risk_score <- predict(fit, newdata = data3)$predicted
    risk_group <- ifelse(risk_score > median(risk_score), "High-Risk", "Low-Risk")
    
    # 生成 `rs` 数据框
    rs <- data.frame(
      Sample_ID = data3$Sample_ID,
      time = data3$time,
      status = data3$status,
      RiskScore = risk_score,
      Risk_group = risk_group
    )
    
    return(rs)
  })
  
  
  # 生成K-M Plot
  output$km_plot <- renderPlot({
    # 计算kmfit
    kmfit <- reactive({
      req(rs_data())
      rs <- rs_data()  # 使用rs_data()获取数据
      km <- survival::survfit(Surv(time = time, event = status) ~ Risk_group, data = rs)
      return(km)
    })
    
    req(rs_data())  # 确保 rs 数据已计算
    rs <- rs_data()  # 获取最新的 rs 数据
    req(kmfit())
    kmfit2 <- kmfit()
    
    # 生成 K-M 曲线
    library(survival)
    library(survminer)
    
    #kmfit <- survival::survfit(Surv(time = rs$time, event = rs$status) ~ Risk_group, data = rs)
    
    # 颜色
    high_color <- ggsci::pal_npg("nrc")(10)[1]
    low_color <- ggsci::pal_npg("nrc")(10)[2]

    # 绘制 K-M 曲线
    plot_surv <- ggsurvplot(fit = kmfit2, data = rs,
                               pval = TRUE,
                               pval.method = TRUE,
                               conf.int = TRUE,
                               risk.table = TRUE,
                               legend.labs = c("High", "Low"),
                               legend.title = "Risk-group",
                               risk.table.col = "strata",
                               linetype = 1,
                               surv.median.line = "hv",
                               ggtheme = theme_bw(base_size = 12) + theme(panel.grid = element_blank()),
                               palette = c(high_color, low_color)
    )
    
    return(plot_surv)
  })
  
  # 下载K-M Plot
    # output$download_km_plot <- downloadHandler(
    #   filename = function() {
    #     paste("km_plot.", input$km_plot_format, sep = "")
    #   },
    #   content = function(file) {
    #     # 将K-M图形保存为文件
    #     ggsave(file, plot = last_plot(), device = input$km_plot_format)
    #   }
    # )
    
    
  # 生成ROC曲线
  output$roc_plot <- renderPlot({
    req(rs_data())
    rs <- rs_data()
    
    library(timeROC)
    
    # 计算1、3、5年生存的ROC曲线
    roc_1year <- timeROC(T = rs$time, delta = rs$status, marker = rs$RiskScore, cause = 1, times = 1, iid = TRUE)
    roc_3year <- timeROC(T = rs$time, delta = rs$status, marker = rs$RiskScore, cause = 1, times = 3, iid = TRUE)
    roc_5year <- timeROC(T = rs$time, delta = rs$status, marker = rs$RiskScore, cause = 1, times = 5, iid = TRUE)
    
    # 绘制 ROC 曲线
    plot(roc_1year, time = 1, col = pal_npg("nrc")(10)[1], title = "Time-dependent ROC", lwd = 2)
    plot(roc_3year, time = 3, add = TRUE, col = pal_npg("nrc")(10)[2], lwd = 2)
    plot(roc_5year, time = 5, add = TRUE, col = pal_npg("nrc")(10)[3], lwd = 2)
    
    # 添加图例
    legend("bottomright", 
           legend = c(
             paste0("1 Year AUC: ", sprintf("%.3f", roc_1year$AUC[2])),
             paste0("3 Year AUC: ", sprintf("%.3f", roc_3year$AUC[2])),
             paste0("5 Year AUC: ", sprintf("%.3f", roc_5year$AUC[2]))
           ),
           col = c(pal_npg("nrc")(10)[1], pal_npg("nrc")(10)[2], pal_npg("nrc")(10)[3]), 
           lty = 1, lwd = 2)
  })
  
  # 下载ROC曲线
  # output$download_roc_plot <- downloadHandler(
  #   filename = function() {
  #     paste("roc_curve.", input$roc_plot_format, sep = "")
  #   },
  #   content = function(file) {
  #     ggsave(file, plot = output$roc_plot , device = input$roc_plot_format)  # 保存为指定格式的图像
  #   }
  # )
  
  # 显示数据表格
  output$dataset_table <- renderDT({
    req(rs_data())
    datatable(rs_data())  # 使用DT包显示数据表格
  })
  
  # 数据表格下载
  # output$download_table <- downloadHandler(
  #   filename = function() {
  #     paste("dataset_table", input$table_format, sep = "")
  #   },
  #   content = function(file) {
  #     req(rs_data())
  #     data <- rs_data()  # 获取最新的 rs 数据
  #     
  #     # 根据选择的格式处理文件
  #     if (input$table_format == ".csv") {
  #       write.csv(data, file, row.names = FALSE, quote = FALSE)
  #     } else if (input$table_format == ".txt") {
  #       write.table(data, file, row.names = FALSE, quote = FALSE, sep = "\t")
  #     } else if (input$table_format == ".xlsx") {
  #       openxlsx::write.xlsx(data, file)  # 使用 file 来保存 XLSX 文件
  #     }
  #   }
  # )
}


# 启动应用----
shinyApp(ui = ui, server = server)
