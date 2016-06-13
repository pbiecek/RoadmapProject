library(shiny)
library(ggplot2)
library(survminer)
library(survival)

load("clinical_expression_mut.rda")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    cohort <- "TCGA Breast Cancer"
    hspgene = "DNAJB2"
    hspgene = input$hspgene
    
    if (input$cohort == "PANCAN12") {
      df <- na.omit(clinical_expression_mut[, c(hspgene), drop=FALSE])
    } else {
      df <- na.omit(clinical_expression_mut[clinical_expression_mut$X_cohort == input$cohort,
                                            c(hspgene), drop=FALSE])
    }

    pl <- ggplot(df, aes_string(hspgene)) +
      geom_histogram()+
      ggtitle(paste(hspgene, "in\n", cohort,"\n", nrow(df), "cases"))
    if (input$median) {
      pl <- pl + geom_vline(xintercept = 0, color="red", size=2) 
    } else { 
      pl <- pl + geom_vline(xintercept = median(df[,hspgene], na.rm = TRUE), color="red", size=2) 
    }
    pl
  })
  
  output$distPlot2 <- renderPlot({
    cohort <- "TCGA Breast Cancer"
    hspgene = "DNAJB2"
    
    if (input$cohort == "PANCAN12") {
      df <- (clinical_expression_mut[, , drop=FALSE])
    } else {
      df <- (clinical_expression_mut[clinical_expression_mut$X_cohort == input$cohort, , drop=FALSE])
    }
    df <- df[df$X18475 == "1", ]
    
    hspgene = input$hspgene
    if (input$median) {
      df$bin <- cut(df[,hspgene], breaks = c(-100,0,100), labels = paste(hspgene,c("low", "high"), "/ mut53"))
    } else { 
      df$bin <- cut(df[,hspgene], breaks = c(-100,median(df[,hspgene], na.rm = TRUE),100), labels = paste(hspgene,c("low", "high"), "/ mut53"))
    }

    df <- na.omit(df[,c("X_TIME_TO_EVENT","X_EVENT","bin")])
    dfClean <- df
    colnames(dfClean)[3] <- "HSP" 
    dfClean$cohort <- input$cohort
    write.table(dfClean, file="www/data/tmp2.csv", sep=";", row.names = F)
    
    model <- survfit(Surv(X_TIME_TO_EVENT,X_EVENT) ~ bin, data=df)
    pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ bin, data=df)$chisq, 1, lower.tail = F)

    ggsurvplot(model, xlim=c(0,4000))$plot + 
      ggtitle(paste0("High/Low ",hspgene, "\nOnly mut tp53\np-value: ", signif(pp,2))) +
      theme(legend.position=c(0.2,0.1)) + xlim(0,4000)+
      scale_color_brewer(type = "qual", palette = 4)
  })
  
  # tylko HIGH
  output$distPlot3 <- renderPlot({
    if (input$cohort == "PANCAN12") {
      df <- (clinical_expression_mut[, , drop=FALSE])
    } else {
      df <- (clinical_expression_mut[clinical_expression_mut$X_cohort == input$cohort, , drop=FALSE])
    }
    
    hspgene = input$hspgene
    if (input$median) {
      df <- df[which(df[,hspgene] > 0),]
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,0,100), labels = paste("MDM2",c("low", "high")))
    } else { 
      df <- df[which(df[,hspgene] > median(df[,hspgene], na.rm = TRUE)),]
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,median(df[,"MDM2"], na.rm = TRUE),100), labels = paste("MDM2",c("low", "high")))
    }
    
    df$TP53 = ifelse(df$X18475 == "1", "TP53 mut", "TP53 wild")
    
    df <- na.omit(df[,c("X_TIME_TO_EVENT", "X_EVENT", "TP53", "MDM2b")])
    df$g <- factor(paste(df$MDM2b, df$TP53))
    
    if (input$groups42) {
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 3, lower.tail = F)
      
      dfClean <- df
      colnames(dfClean)[5] ="Both"
      dfClean$cohort <- input$cohort
      write.table(dfClean[,c(1:4, 6)], file="www/data/tmp3.csv", sep=";", row.names = F)
    } else {
      df$g <- factor(ifelse(grepl(as.character(df$g), pattern = "high.*mut"), paste0("TP53 mut, MDM2 high"), "other"))
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 1, lower.tail = F)
      
      dfClean <- df
      colnames(dfClean)[5] ="Both"
      dfClean$cohort <- input$cohort
      write.table(dfClean[,c(1,2,5,6)], file="www/data/tmp3.csv", sep=";", row.names = F)
    }
    model <- survfit(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)

    ggsurvplot(model, xlim=c(0,4000))$plot + ggtitle(paste0("p-value: ", signif(pp,2))) +
      ggtitle(paste0("Only HIGH ",hspgene, "\np-value: ", signif(pp,2))) +
      theme(legend.position=c(0.3,0.15)) + xlim(0,4000)
  })
  
  # tylko LOW
  output$distPlot4 <- renderPlot({
    if (input$cohort == "PANCAN12") {
      df <- (clinical_expression_mut[, , drop=FALSE])
    } else {
      df <- (clinical_expression_mut[clinical_expression_mut$X_cohort == input$cohort, , drop=FALSE])
    }
    
    hspgene = input$hspgene
    if (input$median) {
      df <- df[which(df[,hspgene] < 0),]
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,0,100), labels = paste("MDM2",c("low", "high")))
    } else { 
      df <- df[which(df[,hspgene] < median(df[,hspgene], na.rm = TRUE)),]
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,median(df[,"MDM2"], na.rm = TRUE),100), labels = paste("MDM2",c("low", "high")))
    }
    
    df$TP53 = ifelse(df$X18475 == "1", "TP53 mut", "TP53 wild")
    
    df <- na.omit(df[,c("X_TIME_TO_EVENT", "X_EVENT", "TP53", "MDM2b")])
    df$g <- factor(paste(df$MDM2b, df$TP53))

    if (input$groups42) {
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 3, lower.tail = F)
      
      dfClean <- df
      dfClean$cohort <- input$cohort
      write.table(dfClean[,c(1,2,3,4,6)], file="www/data/tmp4.csv", sep=";", row.names = F)
    } else {
      df$g <- factor(ifelse(grepl(as.character(df$g), pattern = "high.*mut"), paste0("TP53 mut, MDM2 high"), "other"))
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 1, lower.tail = F)
      
      dfClean <- df
      dfClean$cohort <- input$cohort
      write.table(dfClean[,c(1,2,5,6)], file="www/data/tmp4.csv", sep=";", row.names = F)
    }
    model <- survfit(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)
    pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 3, lower.tail = F)
    
    ggsurvplot(model, xlim=c(0,4000))$plot + ggtitle(paste0("p-value: ", signif(pp,2))) +
      ggtitle(paste0("Only LOW ",hspgene, "\np-value: ", signif(pp,2))) +
      theme(legend.position=c(0.3,0.15)) + xlim(0,4000)
  })
  
  # tylko 
  output$distPlot5 <- renderPlot({
    if (input$cohort == "PANCAN12") {
      df <- (clinical_expression_mut[, , drop=FALSE])
    } else {
      df <- (clinical_expression_mut[clinical_expression_mut$X_cohort == input$cohort, , drop=FALSE])
    }
    
    hspgene = input$hspgene
    if (input$median) {
      df$gene <- cut(df[,hspgene], breaks = c(-100,0,100), labels = paste(hspgene,c("low", "high")))
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,0,100), labels = paste("MDM2",c("low", "high")))
    } else { 
      df$gene <- cut(df[,hspgene], breaks = c(-100,median(df[,hspgene], na.rm = TRUE),100), labels = paste(hspgene,c("low", "high")))
      df$MDM2b <- cut(df[,"MDM2"], breaks = c(-100,median(df[,"MDM2"], na.rm = TRUE),100), labels = paste("MDM2",c("low", "high")))
    }
    
    df$TP53 = ifelse(df$X18475 == "1", "TP53 mut", "TP53 wild")
    
    df <- na.omit(df[,c("X_TIME_TO_EVENT", "X_EVENT", "TP53", "MDM2b", "gene")])
    df$g <- factor(paste(df$MDM2b, df$TP53))

    if (input$groups42) {
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 3, lower.tail = F)
    } else {
      df$g <- factor(ifelse(grepl(as.character(df$g), pattern = "high.*mut"), paste0("TP53 mut, MDM2 high"), "other"))
      pp <- pchisq(survdiff(Surv(X_TIME_TO_EVENT,X_EVENT) ~ g, data=df)$chisq, 1, lower.tail = F)
    }
    
    tmp <- data.frame(table(df$g, df$gene))
    dfClean <- table(df$g, df$gene)
    write.table(dfClean, file="www/data/tmp5.csv", sep=";")
    
    ggplot(tmp, aes(Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity") +
      coord_flip() + theme(legend.position="bottom") + xlab("") + ylab("") +
      scale_fill_brewer(type = "qual")
  })
  
})




# 
# clinical <- clinical.cb[,c("X_PATIENT", "X_TIME_TO_EVENT", "X_EVENT", "X_cohort")]
# clinical$X_PATIENT <- substr((clinical$X_PATIENT), 1, 12)
# clinical$X_PATIENT <- gsub(clinical$X_PATIENT, pattern="-", replacement=".")
# 
# # expression
# expression <- rbind(expression.cb1, expression.cb2)
# rownames(expression) <- expression[,1]
# expression <- expression[,-1]
# expression <- t(expression)
# expression <- as.data.frame(expression)
# expression$X_PATIENT <- substr(rownames(expression), 1, 12)
# expression$X_PATIENT <- gsub(expression$X_PATIENT, pattern="-", replacement=".")
# 
# # mutation
# clinical_expression <- merge(clinical, expression[,c("X_PATIENT", "MDM2", "TP53", "TP63", "TP73", "DNAJB1", "DNAJB2", "DNAJB4", "DNAJB5", "DNAJB6", "DNAJB9", "DNAJB11", "DNAJB12", "DNAJB13", "DNAJB14", "MDM2", "TP73", "TP63",
#                                                      "HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA5", "HSPA6", "HSPA7", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14", "HSP90AA1", "HSP90AB1", "HSP90B1", "TRAP1")],  by="X_PATIENT", all.x=TRUE)
# 
# TP53 <- mutation.cb[grep(mutation.cb[,1], pattern="TP53$", value = FALSE),-1]
# TP53v <- data.frame(X_PATIENT = substr(names(TP53), 1, 12),TP53=t(TP53))
# 
# clinical_expression_mut <- merge(clinical_expression, TP53v, by="X_PATIENT", all.x=TRUE)

# save(clinical_expression_mut, file="clinical_expression_mut.rda")
# 
