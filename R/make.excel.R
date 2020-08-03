MAKE.EXCEL <- function(pathinput,fileinput,contrast,filename, pvalue,padj, FC){


  require(openxlsx)
  require(grDevices)
  ## Read Data
  muestra<-readRDS(file=file.path(pathinput, fileinput))
  # Change GO columns to eliminate the space at the beginning
  if ("GO.BP" %in% colnames(muestra)){
    muestra$GO.BP <- gsub(" ","", muestra$GO.BP)
  }
  if ("GO.CC" %in% colnames(muestra)){
    muestra$GO.CC <- gsub(" ","", muestra$GO.CC)
  }
  if("GO.MF" %in% colnames(muestra)){
    muestra$GO.MF <- gsub(" ","", muestra$GO.MF)
  }
  #Subset columns which shows the PValue and filter pvalues

  #Select heatmap colors columns
  color.values <- muestra[,grep(colnames(muestra),pattern=".scaled",fixed = TRUE)]

  ## Create a new workbook
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "AllData")
  writeData(wb, 1  , muestra)

  if (length(contrast) <= 1){  # We have to separate one contrast from two or more
    if (is.null(padj)){
      muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
      logFC <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
      muestra.final <- muestra[which(muestra.filter.pval <= pvalue & abs(logFC) > FC ),]
      muestra <- muestra.final
    }else{
      muestra.filter.adj<-muestra[,grep(colnames(muestra),pattern="adj.P.Val",fixed = TRUE)]
      logFC <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
      muestra.final <- muestra[which(muestra.filter.adj <= padj & abs(logFC) > FC),]
      muestra <- muestra.final
    }
    #order by FC columns
    muestra.f<-muestra[order(muestra[paste("FC",contrast[[1]][1],"vs",contrast[[1]][2],sep=".")]),]
    # Add new sheet
    addWorksheet(wb, sheetName = paste(contrast[[1]][1],"vs",contrast[[1]][2],sep="."))
    writeData(wb, paste(contrast[[1]][1],"vs",contrast[[1]][2],sep=".") , muestra.f)

  }else{
    for (i in 1:length(contrast)){
      if (is.null(padj)){
        muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
        logFC <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
        muestra.final <- muestra[which(muestra.filter.pval[i] <= pvalue & abs(logFC[i]) > FC ),]
        muestra <- muestra.final
      }else{
        muestra.filter.adj<-muestra[,grep(colnames(muestra),pattern="adj.P.Val",fixed = TRUE)]
        logFC <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
        muestra.final <- muestra[which(muestra.filter.adj[i] <= padj & abs(logFC[i]) > FC ),]
        muestra <- muestra.final
      }
      #order by FC columns
      muestra.f <- muestra[order(muestra[paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")]),]
      # add new sheet
      addWorksheet(wb, sheetName = paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
      writeData(wb, paste(contrast[[i]][1],"vs",contrast[[i]][2],sep=".") , muestra.f)
    }
  }
  sheet.num <- length(contrast)+1
  #per tots els sheets
  for (i in 1:sheet.num){
    # Create several styles for columns and rows
    headerStyle1 <- createStyle(fontSize = 10,halign = "center",textDecoration = "Bold",
                                wrapText = TRUE, textRotation = 90)

    addStyle(wb, sheet = i, headerStyle1, rows = 1, cols = 1:length(colnames(muestra)),
             gridExpand = TRUE)

    bodyStyle1 <- createStyle(fontSize = 10,
                              wrapText = TRUE, valign = "top", halign = "left")
    addStyle(wb, sheet = i, bodyStyle1, rows = 2:length(rownames(muestra)),
             cols = 1:length(colnames(muestra)), gridExpand = TRUE)

    headerStyle2 <- createStyle(fontSize = 8, valign = "center",textDecoration = "Bold",
                                wrapText = TRUE, textRotation = 90)
    addStyle(wb, sheet = i, headerStyle2, rows = 1,
             cols = 1:dim(color.values)[2], gridExpand = TRUE)

    NumberStyle <- createStyle( numFmt = "0.00", fontSize = 10)

    number.col <- ncol(muestra) - dim(color.values)[2]
    FCcols <- grep("FC", colnames(muestra))
    meanCols <- grep("mean", colnames(muestra))
    addStyle(wb, sheet = i, NumberStyle, rows = 2:nrow(muestra),
             cols = number.col : ncol(muestra), gridExpand = TRUE)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:nrow(muestra),
             cols = FCcols, gridExpand = TRUE)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:nrow(muestra),
             cols = meanCols, gridExpand = TRUE)
    # Set Heights and Widths
    setRowHeights(wb, sheet = i, rows = 1, heights = 70)

    setColWidths(wb, sheet = i, cols = ncol(color.values):ncol(muestra), widths = 8)
    setColWidths(wb, sheet = i, cols = 1:ncol(color.values), widths = 3)
    setColWidths(wb, sheet = i, cols =  number.col : ncol(muestra), widths = 5)
    setColWidths(wb, sheet = i, cols =  FCcols, widths = 5)
    setColWidths(wb, sheet = i, cols =  meanCols, widths = 5)

    if ("Description" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Description")
      setColWidths(wb, sheet = i, cols = number.col, widths = 25)
    }
    # Heatmap :
    conditionalFormatting( wb,
                           sheet = i,
                           cols = 1:as.numeric(dim(color.values)[2]),
                           rows = 1:as.numeric(dim(color.values)[1] + 1),
                           rule = as.numeric(c(min(color.values),0,max(color.values))),
                           style = c("blue","white", "red"),
                           type = "colorScale"
    )
  }

  stats<-list()
  #Vector of contrast columns colors
  colors4stats <- c("#FFEA00", "#FFC000", "#00B0F0", "#92D050", "#FF6600", "#CCFF99","#CC99FF", "#FF5252", "#5C45FF", "#45FFC7","#fc79f4","#00B0F0", "#9458d1","#c2a03a", "#d1589b","#b3a7cc","#ccf1ff","#1fad66", "#ffeacc", "#f0a1a1" )
  #només per sheets dels contrasts
  for (i in 1:length(contrast)){
    # Select contrast columns
    stats[[i]] <- grep(colnames(muestra),pattern= paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."),fixed = TRUE)
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1:(nrow(muestra)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i]),
                           type = "contains"
    )# Colouring cols
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = "vs",
                           style = createStyle(bgFill = colors4stats[i]),
                           type = "contains"
    )

  }
  # només pel sheet de ALL DATA:
  for (i in 1:length(contrast)){
    # Select contrast columns
    stats[[i]] <- grep(colnames(muestra),pattern= paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."),fixed = TRUE)
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1:(nrow(muestra)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i]),
                           type = "contains"
    )# Colouring col
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = "vs",
                           style = createStyle(bgFill = colors4stats[i]),
                           type = "contains"
    )

  }

  # Legend:

  addWorksheet(wb, sheetName = "Legend")
  a<-min(color.values)/5
  b<-max(color.values)/5
  vector.min <- c(min(color.values),a*4,a*3,a*2,a)
  vector.max <- c(b,b*2,b*3,b*4, max(color.values))
  legend.df<-as.numeric(c(vector.min, 0 , vector.max))
  options("openxlsx.borderColour" = "black")
  options("openxlsx.borderStyle" = "thin")
  writeData(wb,  "Legend" , legend.df, borderStyle = getOption("openxlsx.borderStyle", "thin")
            ,startCol = 1,startRow = 1, borders = "columns")

  conditionalFormatting( wb,
                         sheet = "Legend",
                         cols = 1,
                         rows = 1:12,
                         rule = as.numeric(c(min(color.values),0,max(color.values))),
                         style = c("blue","white", "red"),
                         type = "colorScale"
  )

  saveWorkbook(wb, file = filename,overwrite = TRUE)

}
