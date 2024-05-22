#social transitions withcontrols. 

# libraries 
library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(organism, character.only = TRUE)
library(DOSE)
library(tidyverse)
grcm38 # mouse genes

my_logFC_threshold = 0.2

#descenders
limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CDOM.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom
cdes <-limma_list$cdes
cdom <- limma_list$cdom

#Ascenders
limma_list<- readRDS("manuscript/brain/results_RDS/limma_MEA_CSUB.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

as <- limma_list$asub
ca <- limma_list$casc
cs <- limma_list$csub


#upreg
dd_up <- dd %>% filter(logFC >= 0.2)
as_up <- as %>% filter(logFC >= 0.2)
cdes_up <- cdes %>% filter(logFC >= 0.2)
ca_up <- ca %>% filter(logFC >= 0.2)
cdom_up <- cdom %>% filter(logFC >= 0.2)
cs_up <- cs %>% filter(logFC >= 0.2)

#downreg
dd_down <- dd %>% filter(logFC <= -0.2)
as_down <- as %>% filter(logFC <= -0.2)
cdes_down <- cdes %>% filter(logFC <= -0.2)
ca_down <- ca %>% filter(logFC <= -0.2)
cdom_down <- cdom %>% filter(logFC <= -0.2)
cs_down <- cs %>% filter(logFC <= -0.2)


# start with des dom controls 
# upsetter Plot 
library(UpSetR)
library(workflowr)
library(ComplexUpset)
listInput <- list(cdes_up= cdes_up$symbol, dd_up=dd_up$symbol, cdes_down = cdes_down$symbol, 
                  dd_down = dd_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

# then asc sub controls 
listInput <- list(ca_up= ca_up$symbol, as_up=as_up$symbol, ca_down = ca_down$symbol, 
                  as_down = as_down$symbol)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)

#getting overlap genes with controls for desenders 
d_up <- cdes_up$symbol[cdes_up$symbol %in% dd_up$symbol] %>% as.data.frame()
d_down <- cdes_down$symbol[cdes_down$symbol %in% dd_down$symbol]%>% as.data.frame()

#getting overlap genes with controls for ascenders 
a_up <- ca_up$symbol[ca_up$symbol %in% as_up$symbol] %>% as.data.frame()
a_down <- ca_down$symbol[ca_down$symbol %in% as_down$symbol]%>% as.data.frame()


#Now go transition genes 
listInput <- list(a_up= a_up$., d_up=d_up$., a_down = a_down$., 
                  d_down = d_down$.)


UpSetR::upset(fromList(listInput), nsets = 4, order.by = "freq", keep.order = F)



#significant test
mat <- matrix(c(5,0,0,3),ncol=2)

# Perform chi-squared test
chi_square <- chisq.test(mat)

# Calculate phi-coefficient manually
n <- sum(mat)
phi <- sqrt(chi_square$statistic / n)

# phi-coefficient and p-value
phi #0.7333333
chi_square$p.value  #0.03806261

#get trans genes with controls 
ctrans_up <- a_up$.[a_up$. %in% d_up$.] %>% as.data.frame()
td_up<- dd[dd$symbol %in% ctrans_up$.,] %>% dplyr::select(symbol, DES_logFC = logFC, DES_pvalue = P.Value)
td_upx<- cdes[cdes$symbol %in% ctrans_up$.,] %>% dplyr::select(symbol, CDES_logFC = logFC, CDES_pvalue = P.Value)
ta_up <- as[as$symbol %in% ctrans_up$., ] %>% dplyr::select(symbol, ASC_logFC = logFC, ASC_pvalue = P.Value)
ta_upx <- ca[ca$symbol %in% ctrans_up$., ] %>% dplyr::select(symbol, CASC_logFC = logFC, CASC_pvalue = P.Value)


ctrans_down <- a_down$.[a_down$. %in% d_down$.] %>% as.data.frame()
td_down<- dd[dd$symbol %in% ctrans_down$.,] %>% dplyr::select(symbol, DES_logFC = logFC, DES_pvalue = P.Value)
td_downx<- cdes[cdes$symbol %in% ctrans_down$.,] %>% dplyr::select(symbol, CDES_logFC = logFC, CDES_pvalue = P.Value)
ta_down <- as[as$symbol %in% ctrans_down$., ] %>% dplyr::select(symbol, ASC_logFC = logFC, ASC_pvalue = P.Value)
ta_downx <- ca[ca$symbol %in% ctrans_down$., ] %>% dplyr::select(symbol, CASC_logFC = logFC, CASC_pvalue = P.Value)


ct_up <- td_up %>% full_join(td_upx) %>%  full_join(ta_up) %>% full_join(ta_upx) 
ct_down <- td_down %>% full_join(td_downx) %>%  full_join(ta_down) %>% full_join(ta_downx) 

ct <- ct_up %>% rbind(ct_down)
# write.csv(ct,"manuscript/brain/results_tables/TransitionGenes_withControls.csv", row.names = F)


library(Hmisc)
##go terms for just dom/sub
#getting overlap genes with controls for desenders 
d_up <- cdes_up$symbol[cdes_up$symbol %in% dd_up$symbol] %>% as.data.frame() %>% unique(.)
dx_up <- dd_up[dd_up$symbol %in% d_up$.,] 
#remove trans genes
dxx_up <- dx_up[dx_up$symbol %nin% ctrans_up$.,]

d_down <- cdes_down$symbol[cdes_down$symbol %in% dd_down$symbol]%>% as.data.frame()
dx_down<- dd_down[dd_down$symbol %in% d_down$.,] 

#remove trans genes
dxx_down <- dx_down[dx_down$symbol %nin% ctrans_down$.,]

go_dx <- dxx_up %>% rbind(dxx_down)



#getting overlap genes with controls for ascenders 
a_up <- ca_up$symbol[ca_up$symbol %in% as_up$symbol] %>% as.data.frame()
ax_up<- as_up[as_up$symbol %in% a_up$.,] 
#remove trans genes
axx_up <- ax_up[ax_up$symbol %nin% ctrans_up$.,]

a_down <- ca_down$symbol[ca_down$symbol %in% as_down$symbol]%>% as.data.frame()
ax_down<- as_down[as_down$symbol %in% a_down$.,] 
#remove trans genes
axx_down <- ax_down[ax_down$symbol %nin% ctrans_down$.,]
go_ax <- axx_up %>% rbind(axx_down)

# #goterms 
# source("functions/gettop10GO.R")
# 
# gettop10GO(go_dx, my_showCategory) %>% 
#   mutate(comparison = "DES_DOM+CDOM") -> top10go1
# 
# gettop10GO(go_ax, my_showCategory) %>% 
#   mutate(comparison = "ASC-SUB+CSUB") -> top10go2
# 
# rbind(top10go1,top10go2) -> top10_GOterms

# write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_Trans_withControls.csv", row.names = F)

