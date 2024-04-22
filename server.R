library(shiny)
library(shinymanager)
library(Seurat)
library(tidyverse)


#### login credentials ####
source('credentials.R')

#### load data ####
plasma_l2 <- readRDS('./data/plasma_l2.rds') %>%
  mutate(Age = factor(Age, levels = c('young', 'old')))

plasma_wilcox <- readRDS('./data/plasma_wilcox.rds')

km_l2 <- readRDS('./data/km_l2.rds') %>%
  mutate(Age = factor(Age, levels = c('young', 'old')))

km_wilcox <- readRDS('./data/km_wilcox.rds')

scdat <- readRDS('./data/allintegrated.rds')

metadata <- readRDS('./data/sc_metadata.rds') %>%
  mutate(label2 = gsub('_cl','_cluster',label2))

allgeneids <- readRDS('./data/geneids.rds') %>%
  mutate(uppersym = toupper(external_gene_name))

geneids <- allgeneids

#### extra functions & setup ####

plot_error <- function(msg = 'gene not found'){
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_label(ggplot2::aes(label = msg, x = 0, y = 0), size = 6, 
                        fill = 'gray75', alpha = 0.25,
                        label.padding = ggplot2::unit(10,'pt'), label.size = 0, 
                        fontface = 'italic') 
}

theme_set(theme_bw())

#### server function ####

shinyServer(function(input, output, session){
  
  #### check credentials ####
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  #### read gene, prep data ####
  
  observeEvent(input$filter, {
    geneids <- allgeneids
    if(nrow(geneids)!=0 & 'km_up' %in% input$filtervar){
      geneids <- filter(km_wilcox, l2FC>0) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    if(nrow(geneids)!=0 & 'km_down' %in% input$filtervar){
      geneids <- filter(km_wilcox, l2FC<0) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    if(nrow(geneids)!=0 & 'km_sig' %in% input$filtervar){
      geneids <- filter(km_wilcox, padj<=0.05) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    if(nrow(geneids)!=0 & 'pl_up' %in% input$filtervar){
      geneids <- filter(plasma_wilcox, l2FC>0) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    if(nrow(geneids)!=0 & 'pl_down' %in% input$filtervar){
      geneids <- filter(plasma_wilcox, l2FC<0) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    if(nrow(geneids)!=0 & 'pl_sig' %in% input$filtervar){
      geneids <- filter(plasma_wilcox, p<=0.05) %>%
        inner_join(geneids) %>%
        select(colnames(geneids)) 
    }
    updateSelectInput(session, "gene",
                      label = 'Gene Symbol',
                      choices = sort(geneids$displaygenes))
  })
  
  gene <- eventReactive(input$submit, input$gene)
  
  plasma_data <- eventReactive(input$submit, {
    filter(geneids, displaygenes %in% gene()) %>%
      select(uppersym) %>% unique() %>%
      left_join(plasma_l2) %>%
      left_join(plasma_wilcox) %>%
      na.omit()
  })
  
  km_data <- eventReactive(input$submit, {
    filter(geneids, displaygenes %in%  gene()) %>%
      select(uppersym) %>% unique() %>%
      left_join(km_l2) %>%
      left_join(km_wilcox) %>%
      na.omit()
  })
  
  sc_data <- eventReactive(input$submit, {
    scgene <- unique(filter(geneids, displaygenes %in%  gene())$datagenes)
    if(!is.na(scgene)){
      scgenes <- grep(paste('^',scgene,sep=''),rownames(scdat@assays$RNA),v=T)
      } else {scgenes = c()}
    if(length(scgenes)>1){
      xx = as.matrix(scdat@assays$RNA[scgenes,])
      reshape2::melt(apply(xx,2,median)) %>% 
        set_names('Normalized Expression') %>%
        mutate(cell = colnames(xx)) %>%
        left_join(metadata) %>%
        mutate(SampleID = setNames(c('Young1', 'Young2', 'Old1', 'Old2', 'Old3'),
                                   c('A', 'C', 'D', 'E', 'F'))[as.character(SampleID)]) %>%
        mutate(SampleID = factor(SampleID, 
                                 levels = c('Young1', 'Young2', 'Old1', 'Old2', 'Old3')))
    } else if(length(scgenes)==1){
      reshape2::melt(as.matrix(scdat[['RNA']]$data[scgenes,])) %>%
        set_names('cell', 'Gene', 'Normalized Expression') %>%
        left_join(metadata) %>%
        mutate(SampleID = setNames(c('Young1', 'Young2', 'Old1', 'Old2', 'Old3'),
                                   c('A', 'C', 'D', 'E', 'F'))[as.character(SampleID)]) %>%
        mutate(SampleID = factor(SampleID, 
                                 levels = c('Young1', 'Young2', 'Old1', 'Old2', 'Old3')))
      
    } else { data.frame()}
  })
  
  tsnedata <- eventReactive(input$plot_click, {
    data.frame(scdat@reductions$tsne@cell.embeddings) %>%
      mutate(cell = rownames(scdat@reductions$tsne@cell.embeddings))  %>%
      left_join(sc_data())
    })
  
  #### plasma boxplot ####
  output$plasma_box <- renderPlot({
    if(nrow(plasma_data())>0){
      plasma_data() %>%
        ggplot(aes(x = Age, y = log2Expression)) +
        geom_boxplot(outlier.shape = NA, fill = 'gray70', width = 0.3) +
        geom_jitter(width = 0.1, size = 2) +
        ggpubr::stat_compare_means() +
        ylab('Log 2 Protein Expression') + xlab(NULL) +
        ggtitle('Plasma Proteomics',
                paste('Log2 Fold Change = ',
                      round(unique(plasma_data()$l2FC), 2), sep=''))
      
    } else {
      return(plot_error(paste(gene(),'\nis not found in\nplasma proteomics.')))
      }
  })
  
  
  #### km boxplot ####
  output$km_box <- renderPlot({
    if(nrow(km_data())>0){
      km_data() %>%
        ggplot(aes(x = Age, y = log2Expression)) +
        geom_boxplot(outlier.shape = NA, fill = 'gray70', width = 0.3) +
        geom_jitter(width = 0.1, size = 2) +
        ggpubr::stat_compare_means() +
        ylab('Log 2 Protein Expression') + xlab(NULL) +
        ggtitle('Kidney Marrow Proteomics',
                paste('Log2 Fold Change = ',
                      round(unique(km_data()$l2FC), 2), sep=''))
      
    } else{return(plot_error(paste(gene(),'\nis not found in\nkidney marrow proteomics.')))}
  })
  
  
  #### tsne ####
  output$tsne <- renderPlot({
    if(nrow(sc_data())>0){
      
      tsnedat <- data.frame(scdat@reductions$tsne@cell.embeddings) %>%
        mutate(cell = rownames(scdat@reductions$tsne@cell.embeddings))
     
       tsnelabels = left_join(metadata,tsnedat) %>%
        group_by(seurat_clusters) %>%
        summarise(tSNE_1 = median(tSNE_1),
                  tSNE_2 = median(tSNE_2))
     
       tsnedat %>%
        left_join(sc_data()) %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2, color = `Normalized Expression`)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_color_gradient(high = '#800026', low = "#FFFFCC") +
        geom_text(data = tsnelabels, aes(label = seurat_clusters),
                  color = 'black',
                  size = 4) +
        guides(color = guide_colorbar(title= 'Normalized\nExpression**')) +
        theme_classic() +
        ggtitle('Kidney Marrow scRNAseq*')
    } else { 
      return(plot_error(paste(gene(),'\nis not found in\nkidney marrow scRNAseq.'))) 
      }
  })
  
  output$info <- renderText({
    names(sort(table(nearPoints(tsnedata(), 
                                input$plot_click, 
                                xvar = "tSNE_1", 
                                yvar = "tSNE_2")$label2),dec=T))[1]
  })
  
  
  #### dotplot ####
  output$dotplot <- renderPlot({
    if(nrow(sc_data())>0){
      sc_data() %>%
        group_by(Gene, SampleID, label2) %>%
        summarise(meanExp = mean(`Normalized Expression`),
                  percExp = mean(`Normalized Expression`>0)*100) %>%
        ggplot(aes(x = SampleID, y = label2, fill = meanExp, size = percExp)) +
        geom_point(shape=21) +
        scale_fill_gradient(high = '#800026', low = "#FFFFCC") +
        ylab(NULL) + xlab(NULL) +
        geom_vline( xintercept = 2.5) +
        guides(fill = guide_colorbar( title = 'Mean\nNormalized\nExpression'),
               size = guide_legend( title = '% Expressed')) +
        coord_flip() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
        ggtitle('Kidney Marrow scRNAseq')
      
    } 
  }
  )
  
  
})