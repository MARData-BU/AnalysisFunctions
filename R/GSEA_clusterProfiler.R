##' Run GSEA for a list of contrasts with cluterProfiler wrapper for fgsea 
##'
##' Function that runs gene-set enrichment analysis (GSEA) using the function GSEA
##' from clusterProfiler package. This function calls fgsea package, a fast implementation 
##' of GSEA.
##' @param data4Tyers data4tyers dataframe. Must contain a Geneid column and 
##' P.Value and logFC columns for each contrast (eg. "P.Value.LES_ACT.vs.LES_INACT")
##' @param contrast List with one vector of length 2 for each contrast.
##' @param gmt Data frame with gene sets. Can be obtained with function 
##' clusterProfiler::read.gmt. Must have a column "term" and a column "gene".
##' @param resultsDir Character vector with output results directory. Default is working directory.
##' @param specie Character vector with the name of the specie. If not "human", 
##' must provide an homology file (HOM_MouseHuman). 
##' @param collection_name Name of the collection to append to output files.
##' @param minGSSize minimal size of each geneSet for analyzing (default 15)
##' @param maxGSSize maximal size of genes annotated for testing (Default 500)
##' @param pvalueCutoff adjusted pvalue cutoff (default 1 to return all results)
##' @param HOM_MouseHuman dataframe with HOM_MouseHumanSequence.txt. Must have 
##' the following columns: HomoloGene.ID (common ID for homologues), Common.Organism.Name, Symbol
##' @return Returns a list with gsea results for each contrast. The list is saved An RData object 
##' in the resultsDir directory, along with a folder for each contrast containing an 
##' excel file with enrichment results.
##' @author Julia Perera Bel <jperera@imim.es>
##' @export
##' @import clusterProfiler
##' @import openxlsx 
GSEA.run <- function(data4Tyers,contrast,gmt,resultsDir=getwd(),specie="human",
                     collection_name="",minGSSize = 15,maxGSSize = 500,pvalueCutoff=1,HOM_MouseHuman=""){
  
  require(clusterProfiler)
  require(openxlsx)
  headerStyle1 <- createStyle(halign = "center",valign = "center",textDecoration = "Bold",
                              wrapText = TRUE) #textRotation=90
  # Transform species homologues to human Symbols
  if (specie!="human"){
    #HOM_MouseHuman <- read.table(file="/bicoh/MARGenomics/annotationData/HOM_MouseHumanSequence.txt",sep="\t",stringsAsFactors = F, header=T)
    # Divide human and mouse data
    Human.HOM <- HOM_MouseHuman[HOM_MouseHuman$Common.Organism.Name == "human",] #19124    13
    Mouse.HOM <- HOM_MouseHuman[HOM_MouseHuman$Common.Organism.Name == "mouse, laboratory",] #20943    13
    
    sum(data4Tyers$Geneid %in% Mouse.HOM$Symbol)#14080 of 20943 of data4tyers Geneids are in mouse hom. DB
    
    # Subset mouse DB to data4tyers Geneids Seleccionem els ID de mouse que es troben al nostre data4Tyers
    Mouse.ID <- Mouse.HOM[Mouse.HOM$Symbol %in% data4Tyers$Geneid,] #14080    13
    
    # Change colnames to make them unique for further merging (except first column "HomoloGene.ID")
    colnames(Mouse.ID)[2:ncol(Mouse.ID)] <- paste(colnames(Mouse.ID)[2:ncol(Mouse.ID)], "Mouse", sep=".")
    colnames(Human.HOM)[2:ncol(Human.HOM)] <- paste(colnames(Human.HOM)[2:ncol(Human.HOM)], "Human", sep=".")
    
    Hum.nd.Mouse <- merge(Mouse.ID, Human.HOM, all.x = T, all.y=F, sort=FALSE)
    dim(Hum.nd.Mouse)# 14126    25
    
    #Hi ha simbols repetits en Mouse i Human
    # Hum.nd.Mouse$Symbol.Mouse[duplicated(Hum.nd.Mouse$Symbol.Mouse)]
    # Hum.nd.Mouse$Symbol.Human[duplicated(Hum.nd.Mouse$Symbol.Human)]
    # Remove NA's from human
    Hum.nd.Mouse <- Hum.nd.Mouse[!is.na(Hum.nd.Mouse$Symbol.Human),]
    dim(Hum.nd.Mouse)#13538    25
    
    Hum.nd.Mouse$Geneid <- Hum.nd.Mouse$Symbol.Mouse
    dim(data4Tyers) #17668    99
    data4Tyers.HOM <- merge(data4Tyers[,c("Geneid",
                                          grep("logFC", colnames(data4Tyers), value=T),
                                          grep("P.Value", colnames(data4Tyers), value=T),
                                          grep("adj.P.Val", colnames(data4Tyers), value=T))],
                            Hum.nd.Mouse[,c("Geneid","Symbol.Human")])
    dim(data4Tyers.HOM)#13538    26
    #Hi ha simbols repetits en Mouse i Human
    # data4Tyers.HOM$Geneid[duplicated(data4Tyers.HOM$Geneid)]
    # data4Tyers.HOM$Symbol.Human[duplicated(data4Tyers.HOM$Symbol.Human)]
    
    save(data4Tyers.HOM, file=file.path(resultsDir,"data4Tyers.HOM.RData"))
    #load(file.path(dir,"data4Tyers.HOM.RData"))
    data4Tyers=data4Tyers.HOM
    colnames(data4Tyers)[colnames(data4Tyers)=="Geneid"]="Geneid.mouse" 
    colnames(data4Tyers)[colnames(data4Tyers)=="Symbol.Human"]="Geneid" 
    
  }
  
  # Create object to return gsea results
  gsea=list()
  dir=resultsDir
  for (i in 1:length(contrast)){
    # Create a folder for each contrast
    resultsDir=file.path(dir,paste("GSEA",collection_name,contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
    dir.create(resultsDir,showWarnings = F)
    ### Prepare geneList from data4tyers
    p=paste("P.Value",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    logfc=paste("logFC",contrast[[i]][1],"vs",contrast[[i]][2],sep=".")
    
    sign=sign(data4Tyers[,logfc])
    logP=-log10(data4Tyers[,p])
    metric=logP*sign
    geneList = metric
    names(geneList) = data4Tyers$Geneid
    geneList = sort(geneList, decreasing = TRUE)
    
    ### GSEA
    set.seed(123)
    gsea[[i]] <- clusterProfiler::GSEA(geneList, TERM2GENE=gmt, verbose=FALSE,seed = T,
                      minGSSize = minGSSize,maxGSSize = maxGSSize,pvalueCutoff=pvalueCutoff) # Volem tots els resultats
    
    ## Save excel with results
    wb <- createWorkbook()
    for (j in 1:length(contrast[[i]])){ # Loop through positive and negative
      # Positive NES
      addWorksheet(wb, sheetName = paste("Enriched in",contrast[[i]][j]))
      if(j==1){writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES>0,-c(2,8)])} # POSITIVE NES
      if(j==2){writeData(wb, j, gsea[[i]]@result[gsea[[i]]@result$NES<0,-c(2,8)])} # NEGATIVE NES
      setColWidths(wb, sheet = j, cols = 1:9, widths =c(35,6,10,10,10,10,5,15,25))
      addStyle(wb, sheet = j, headerStyle1, rows = 1, cols = 1:9, gridExpand = TRUE)
      setRowHeights(wb, sheet = j, rows = 1, heights =30)
    }
    
    saveWorkbook(wb, file.path(resultsDir, paste("GSEA",collection_name,contrast[[i]][1],"vs",contrast[[i]][2],"xlsx",sep=".")), overwrite = TRUE)
  }
  resultsDir=dir
  save(gsea,file=file.path(resultsDir, paste0("GSEA",collection_name,".RData")))
  return(gsea)
}


##' Create plots summarizing GSEA results
##'
##' Function that creates barplots, dotplots, GSEA plots, gene-concept networks 
##' and enrichmentMAP from clusterProfiler and enrichplot packages. 
##' @param gsea GSEA results obtained with GSEA.run function
##' @param contrast List with one vector of length 2 for each contrast.
##' @param collection_name Name of the collection to append to output files.
##' @param resultsDir Character vector with output results directory. Default is working directory.
##' @param plot_top maximum number of GSEA plots (default 50 to return all results)
##' @param p.adjust adjusted pvalue cutoff (default 0.05). Gene sets with p.adjust 
##' below this threshold will be included in all plots, unless maximum for visualization is reached.
##' @return This function creates a folder for each contrast and generates barplots, 
##' dotplots, GSEA plots, gene-concept networks and enrichmentMAP for both phenotypes.
##' @author Julia Perera Bel <jperera@imim.es>
##' @export
##' @import clusterProfiler
##' @import openxlsx 
##' @import gridExtra
##' @import png 
GSEA.plots <- function(gsea,contrast,collection_name="",resultsDir=getwd(),
                       plot_top=10,p.adjust=0.05,make_cnet=F){
  require(clusterProfiler)
  require(enrichplot)
  require(gridExtra)
  require(png)
  require(Cairo)
  options(bitmapType = "cairo")
  dir=resultsDir
  for (i in 1:length(contrast)){
    cat("Plotting", paste(contrast[[i]][1],"vs",contrast[[i]][2],sep="."),"\n")
    resultsDir=file.path(dir,paste("GSEA",collection_name,contrast[[i]][1],"vs",contrast[[i]][2],sep="."))
    dir.create(resultsDir,showWarnings = F)
    # Create a list with significant [[1]] positive and [[2]] negative results
    gsea.L=list()
    gsea.L[[contrast[[i]][1]]]=gsea[[i]]
    gsea.L[[contrast[[i]][1]]]@result=gsea.L[[contrast[[i]][1]]]@result[gsea.L[[contrast[[i]][1]]]@result$NES>0&gsea.L[[contrast[[i]][1]]]@result$p.adjust<p.adjust,]
    gsea.L[[contrast[[i]][2]]]=gsea[[i]]
    gsea.L[[contrast[[i]][2]]]@result=gsea.L[[contrast[[i]][2]]]@result[gsea.L[[contrast[[i]][2]]]@result$NES<0&gsea.L[[contrast[[i]][2]]]@result$p.adjust<p.adjust,]
    
    # Generate a general Barplot with NES and pval (maximum to top 30 up and down)
    up=gsea.L[[contrast[[i]][1]]]@result[,c("ID","NES","p.adjust")]
    up=up[1:ifelse(nrow(up)>30,30,nrow(up)),]
    down=gsea.L[[contrast[[i]][2]]]@result[,c("ID","NES","p.adjust")]
    down=down[1:ifelse(nrow(down)>30,30,nrow(down)),]
    data=na.omit(rbind(up,down))
    data$ID=strtrim(data$ID, 70) # maximum label length
    data$ID=factor(data$ID,data$ID[order(data$NES)],ordered=T)
    if(nrow(data)!=0){
      p=ggplot(data, aes(x=ID, y=NES,fill=p.adjust)) + 
        geom_bar(stat = "identity") +
        scale_fill_continuous(low='red', high='blue')+
        ylim(ifelse(min(data$NES)< (-2),min(data$NES),-2),
             ifelse(max(data$NES)> (2),max(data$NES),2)) +
        coord_flip() +
        theme_bw() +
        ggtitle(paste0("GSEA of ", collection_name)) +
        theme(plot.title = element_text(hjust = 0.5,face="bold"))
      
      ggsave(file.path(resultsDir, 
                       paste0("GSEA.",collection_name,".BarplotNES.",
                              paste0(contrast[[i]][1],"vs",contrast[[i]][2]),".png")),
             plot=p,width = 10,height =8*(ceiling(nrow(data)/80)))
      
    }else{
      png::writePNG(array(0, dim = c(1,1,4)), 
                    file.path(resultsDir,
                              paste0("GSEA.",collection_name,".BarplotNES.",
                                     paste0(contrast[[i]][1],"vs",contrast[[i]][2]),".png")))
    }
    for (j in 1:length(contrast[[i]])){ # Loop through positive and negative
      if(nrow(gsea.L[[j]]@result)!=0){ # If there are significant results; do plots
        # GSEA Dotplot
        gsea.L[[j]]@result$Description=strtrim(gsea.L[[j]]@result$Description, 70) # maximum label length
        
        p=clusterProfiler::dotplot(gsea.L[[j]], showCategory=40,font.size=8,title= paste0("Enriched in ",names(gsea.L)[j],"\n p.adjust<0.05"))
        
        ggsave(file.path(resultsDir, paste0("GSEA.",collection_name,".Dotplot.", names(gsea.L)[j],".png")),
               plot=p,width =9,height =7)
        
        # GSEA Running score
        p=list()
        # Create GSEA dir for Running Score results
        GSEA_score_dir=paste0(resultsDir,"/","RunningScore");dir.create(GSEA_score_dir)
        ## Set max num. of plots to all DE or to plot_top
        n_plot=ifelse(nrow(gsea.L[[j]]@result)>plot_top,plot_top,nrow(gsea.L[[j]]@result))
        for (u in 1:n_plot){
          p=enrichplot::gseaplot2(gsea.L[[j]],color="green", base_size = 4.5, 
                           geneSetID = u,title=gsea.L[[j]]$Description[u])
          print(p)
          ggsave(filename = file.path(GSEA_score_dir, paste0(gsea.L[[j]]$Description[u],".RunningScore.",names(gsea.L)[j],".png")),
                 width = 3,height =2.5)
        }

        # GSEA Gene-Concept networks (top 5 gene sets by default)
        if(make_cnet==TRUE){
          p=clusterProfiler::cnetplot(gsea.L[[j]], foldChange=gsea.L[[j]]@geneList,cex_label_gene = 0.5,
                     cex_label_category = 0.7,cex_category = 0.7,layout = "kk",showCategory = 5)
          p=p+ scale_color_gradient2(name = "-log(p.val)*signFC", low = "blue", mid = "white", high = "red")
          ggsave(file.path(resultsDir, paste0("GSEA.",collection_name,".GeneConceptNetworks.", 
                                              names(gsea.L)[j], ".png")), plot=p)
        }
        # GSEA EnrichmentMAP (Jaccard index. Plot top 30 by default)
        if(nrow(gsea.L[[j]]@result)>1){
          pt=enrichplot::pairwise_termsim(gsea.L[[j]], method = "JC", semData = NULL, showCategory = 200)
          p <- clusterProfiler::emapplot(pt,cex_label_category = 0.5,showCategory = 30)
          ggsave(file.path(resultsDir, paste0("GSEA.",collection_name,".EnrichmentMAP.", 
                                              names(gsea.L)[j], ".png")),plot=p)
        }else{png::writePNG(array(0, dim = c(1,1,4)), file.path(resultsDir, paste0("GSEA.",collection_name,".EnrichmentMAP.", 
                                                                                   names(gsea.L)[j], ".png")))}
      }else{
        png::writePNG(array(0, dim = c(1,1,4)), file.path(resultsDir, paste0("GSEA.",collection_name,".Dotplot.", names(gsea.L)[j],".png")))
        png::writePNG(array(0, dim = c(1,1,4)), file.path(resultsDir, paste0("GSEA.",collection_name,".GeneConceptNetworks.", 
                                                                             names(gsea.L)[j], ".png")))
        png::writePNG(array(0, dim = c(1,1,4)), file.path(resultsDir, paste0("GSEA.",collection_name,".EnrichmentMAP.", 
                                                                             names(gsea.L)[j], ".png")))
        
      }
    }
  }
}

