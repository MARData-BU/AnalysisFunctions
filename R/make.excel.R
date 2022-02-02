make.excel <- function(pathinput,fileinput,contrast,pathoutput,filename, pvalue = NULL,padj = 0.05, logFC = 1,add.colors=NULL){
  
  require(wrapr)
  require(openxlsx)
  require(grDevices)
  ## Read Data
  message("Read RDS file ...")
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
  #Select heatmap colors columns
  color.values <- muestra[,grep(colnames(muestra),pattern=".scaled",fixed = TRUE)]
  
  ## Create a new workbook
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "AllData")
  
  writeData(wb, 1  , muestra)
  
  
  ########### FILTER ############
  
  message("Filter ...")
  if (length(contrast) <= 1){  # We have to separate one contrast from two or more
    vector.adj = c()
    vector.noadj = c()
    if (is.null(pvalue)){ # Pvalue = NULL
      muestra.filter.adj<-muestra[,grep(colnames(muestra),pattern="adj.P.Val",fixed = TRUE)]
      logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
      if(!any(muestra.filter.adj < padj)){  
        warning("PValue = 0.05 was used to filter data.")
        muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
        logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
        muestra.final <- muestra[which(muestra.filter.pval <= 0.05 & abs(logFC.col) > logFC ),]
        #order by FC columns
        muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[1]][1],"vs",contrast[[1]][2],sep=".")]),]
        # add new sheet
        addWorksheet(wb, sheetName = paste(contrast[[1]][1],"vs",contrast[[1]][2],sep="."))
        writeData(wb, paste(contrast[[1]][1],"vs",contrast[[1]][2],sep=".") , muestra.f)
      }else{ # if muestra.filter.adj < padj : 
        warning("Padj = 0.05 was used to filter data.")
        muestra.final <- muestra[which(muestra.filter.adj <= padj & abs(logFC.col) > logFC),]
        #order by FC columns
        muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[1]][1],"vs",contrast[[1]][2],sep=".")]),]
        # add new sheet
        addWorksheet(wb, sheetName = paste(contrast[[1]][1],"vs",contrast[[1]][2],sep="."))
        writeData(wb, paste(contrast[[1]][1],"vs",contrast[[1]][2],sep=".") , muestra.f)
      }
    }else{ # if pvalue != NULL
      warning("Pvalue = 0.05 was used to filter data.")
      muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
      logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
      muestra.final <- muestra[which(muestra.filter.pval <= pvalue & abs(logFC.col) > logFC ),]
      #order by FC columns
      muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[1]][1],"vs",contrast[[1]][2],sep=".")]),]
      # add new sheet
      addWorksheet(wb, sheetName = paste(contrast[[1]][1],"vs",contrast[[1]][2],sep="."))
      writeData(wb, paste(contrast[[1]][1],"vs",contrast[[1]][2],sep=".") , muestra.f)
    }
  }else{ #if length contrast > 1
    vector.adj = c()
    vector.noadj = c()
    for (i in 1:length(contrast)){
      if (is.null(pvalue)){ # Pvalue = NULL
        muestra.filter.adj<-muestra[,grep(colnames(muestra),pattern="adj.P.Val",fixed = TRUE)]
        logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
        if(!any(muestra.filter.adj < padj)){  
          warning("PValue = 0.05 was used to filter data.")
          muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
          logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
          muestra.final <- muestra[which(muestra.filter.pval[i] <= 0.05 & abs(logFC.col[i]) > logFC ),]
          #order by FC columns
          muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")]),]
          # add new sheet
          addWorksheet(wb, sheetName = paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
          writeData(wb, paste(contrast[[i]][1],"vs",contrast[[i]][2],sep=".") , muestra.f)
        }else{ # if muestra.filter.adj < padj : 
          warning("Padj = 0.05 was used to filter data.")
          muestra.final <- muestra[which(muestra.filter.adj[i] <= padj & abs(logFC.col[i]) > logFC),]
          #order by FC columns
          muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")]),]
          # add new sheet
          addWorksheet(wb, sheetName = paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
          writeData(wb, paste(contrast[[i]][1],"vs",contrast[[i]][2],sep=".") , muestra.f)
        }
      }else{ # if pvalue != NULL
        warning("Pvalue = 0.05 was used to filter data.")
        muestra.filter.pval<-muestra[,grep(colnames(muestra),pattern="P.Value",fixed = TRUE)]
        logFC.col <- muestra[,grep(colnames(muestra),pattern="logFC",fixed = TRUE)]
        muestra.final <- muestra[which(muestra.filter.pval[i] <= pvalue & abs(logFC.col[i]) > logFC ),]
        #order by FC columns
        muestra.f <- muestra.final[orderv(muestra.final[paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")]),]
        # add new sheet
        addWorksheet(wb, sheetName = paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
        writeData(wb, paste(contrast[[i]][1],"vs",contrast[[i]][2],sep=".") , muestra.f)
      }
    }
  }
  
  message("Style ...")
  
  sheet.num <- length(contrast)+1
  # For all sheets 
  for (i in 1:sheet.num){
    # Create several styles for columns and rows
    headerStyle1 <- createStyle(fontSize = 10,textDecoration = "Bold",
                                wrapText = TRUE,textRotation = 90,halign = "center")
    
    addStyle(wb, sheet = i, headerStyle1, rows = 1, cols = 1:length(colnames(muestra)),
             gridExpand = TRUE)
    
    bodyStyle1 <- createStyle(fontSize = 10,
                              wrapText = TRUE, valign = "top", halign = "left")
    addStyle(wb, sheet = i, bodyStyle1, rows = 2:(length(rownames(muestra))+1),
             cols = 1:length(colnames(muestra)), gridExpand = TRUE)
    
    headerStyle2 <- createStyle(fontSize = 10, halign = "center",textDecoration = "Bold",
                                wrapText = TRUE, textRotation = 90)
    addStyle(wb, sheet = i, headerStyle2, rows = 1,
             cols = 1:dim(color.values)[2], gridExpand = TRUE)
    
    NumberStyle <- createStyle( fontSize = 10, numFmt = "0.00")
    
    NumberStyle_adj.pval <- createStyle (fontSize = 10, numFmt = "SCIENTIFIC")
    
    HeatmapStyle <- createStyle(fontSize = 10, numFmt = "0")
    addStyle(wb,sheet=i, HeatmapStyle, rows = 2:(nrow(muestra)+1) , cols = 1:(dim(color.values)[2]),gridExpand=T)
    number.col <- ncol(muestra) - dim(color.values)[2]
    FCcols <- grep("FC", colnames(muestra))
    meanCols <- grep("mean", colnames(muestra))
    adjPvalCols <- grep("adj.P.Val", colnames(muestra))
    scaleCols <- grep(".scaled", colnames(muestra))
    
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(muestra)+1),
             cols = number.col : ncol(muestra),gridExpand=T)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(muestra)+1),
             cols = FCcols,gridExpand=T)
    addStyle(wb, sheet = i, NumberStyle, rows = 2:(nrow(muestra)+1),
             cols = meanCols,gridExpand=T)
    addStyle(wb, sheet = i , NumberStyle_adj.pval, rows = 2:(nrow(muestra)+1), cols = adjPvalCols,gridExpand=T)
    # Set Heights and Widths
    setRowHeights(wb, sheet = i, rows = 1, heights = 150)
    setRowHeights(wb, sheet = i, rows = 2:(nrow(muestra)+1), heights = 14)
    
    setColWidths(wb, sheet = i, cols = (ncol(color.values)-1):ncol(muestra), widths = 8)
    setColWidths(wb, sheet = i, cols = 1:ncol(color.values), widths = 2)
    setColWidths(wb, sheet = i, cols =  number.col : ncol(muestra), widths = 5)
    setColWidths(wb, sheet = i, cols =  FCcols, widths = 5)
    setColWidths(wb, sheet = i, cols =  meanCols[1]:ncol(muestra), widths = 6)
    # Fix first row
    freezePane(wb, sheet = i , firstRow = TRUE)
    
    # Change some widths according to specific columns
    if ("AffyID" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "AffyID")
      setColWidths(wb, sheet = i, cols = number.col, widths = 18)
    }
    if ("Symbol" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Symbol")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("Geneid" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Geneid")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("mrna" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "mrna")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("UCSC_symbols" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "UCSC_symbols")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GOBP" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GOBP")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO.BP" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO.BP")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO:BP" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO:BP")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GOCC" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GOCC")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO.CC" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO.CC")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO:CC" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO:CC")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GOMF" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GOMF")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO.MF" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO.MF")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("GO:MF" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "GO:MF")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("^path" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "^path")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("Description" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Description")
      setColWidths(wb, sheet = i, cols = number.col, widths = 40)
    }
    if ("Length" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Length")
      setColWidths(wb, sheet = i, cols = number.col, widths = 6)
    }
    if ("Strand" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Strand")
      setColWidths(wb, sheet = i, cols = number.col, widths = 3)
    }
    if ("Chr" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Chr")
      setColWidths(wb, sheet = i, cols = number.col, widths = 5)
    }
    if ("Chrom" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Chrom")
      setColWidths(wb, sheet = i, cols = number.col, widths = 4)
    }
    if ("Start" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Start")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
    }
    if ("Stop" %in% colnames(muestra)){
      number.col<-which(colnames(muestra) == "Stop")
      setColWidths(wb, sheet = i, cols = number.col, widths = 8)
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
  
  # COLOURING CONTRASTS: 
  stats<-list()
  #Vector of contrast columns colors
  colors4stats <- c(c("#FFEA00", "#FFC000", "#00B0F0", "#92D050", "#FF6600", "#CCFF99","#CC99FF", "#FF5252", "#5C45FF", "#45FFC7","#fc79f4","#00B0F0", "#9458d1","#c2a03a", "#d1589b","#b3a7cc","#ccf1ff","#1fad66", "#ffeacc", "#f0a1a1" ),add.colors)
  #Only for sheets in contrasts
  for (i in 1:length(contrast)){
    adjp=paste("adj.P.Val",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    fc=paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    # Select contrast columns
    stats[[i]] <- which(colnames(muestra) %in% c(adjp,p,fc,logfc))
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1:(nrow(muestra)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )# Colouring cols
    conditionalFormatting( wb,
                           sheet = i+1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )
    
  }
  # Only for sheet ALL DATA:
  for (i in 1:length(contrast)){
    adjp=paste("adj.P.Val",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    fc=paste("FC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    # Select contrast columns
    stats[[i]] <- which(colnames(muestra) %in% c(adjp,p,fc,logfc))
    # Colouring rows
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1:(nrow(muestra)+1),
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )# Colouring col
    conditionalFormatting( wb,
                           sheet = 1,
                           cols = stats[[i]],
                           rows = 1,
                           rule = ".",
                           style = createStyle(bgFill = colors4stats[i], fontSize = 10),
                           type = "contains"
    )
    
  }
  
  # Legend:
  message("Legend ...")
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
  saveWorkbook(wb, file.path(pathoutput,paste(filename, "xlsx", sep=".")),overwrite = TRUE)
  message("The function was performed successfully")
}



# 09/02/2021: change pvalue of contrast = 1
# change letter size
# 10/02/2021: change widths of columns 
# 21/10/2021: change widths of columns 
# 02/02/2022: with the updated version of R, appeared a warning message. We could solve it replacing order() by orderv() : expalantion -> https://www.r-bloggers.com/2021/02/it-has-always-been-wrong-to-call-order-on-a-data-frame/
