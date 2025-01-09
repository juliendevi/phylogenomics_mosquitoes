############################################################################

############################################################################
##################        Figures phylogenomics          ###################
############################################################################

############################################################################





############################################################################
################        Trees and ancestral trait          #################
############################################################################


# Libraries ----
library(RRphylo)   # CRAN v2.8.0
library(phytools)  # CRAN v2.0-3
library(geiger)    # CRAN v2.0.11
library(MetBrewer) # CRAN v0.2.0
library(ggtree)
library(readxl)

# Import and setup the tree ----
tree = read.newick('~/Desktop/data/CAFE/astral/tree.tre')

for (i in 1:length(tree$tip.label)){
  tree$tip.label[i]=unlist(strsplit(tree$tip.label[i], "[<>]"))[1]
}
for (i in 1:length(tree$node.label)){
  tree$node.label[i]=unlist(strsplit(tree$node.label[i], "[<>]"))[2]
}

informations=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='R_table')

for (i in 1:length(tree$tip.label)){
  tree$tip.label[i]=informations$new_name[which(tree$tip.label[i]==informations$tree_name)]
}


# Run ancestral trait reconstruction ----
biting_astral=cbind(informations$biting, informations$nonbiting)
rownames(biting_astral)=informations$new_name
name.check(tree,biting_astral)
trait_astral=RRphylo(tree, biting_astral)

# Visualisation ancestral states ----

# Prepare tree and table for plot
p=ggtree(trait_astral$tree)+geom_tiplab(aes(x=x+3))+xlim(c(0,450))
# Extract ancestral states
ances_astral=rbind(biting_astral, trait_astral$aces)
lab=c(rownames(ances_astral)[1:64],p$data$label[65:126])
ances_astral=as.data.frame(cbind(lab, ances_astral))
# Extract the final likely state
ances_astral[,4]=apply(ances_astral[,2:3], 1, which.max)
# reorder as tree in p
end=ances_astral[match(ances_astral[,1], p$data$label),]
# join tree and table
p2= p %<+% end

offsetlab=125

p2 + scale_color_manual(values=met.brewer("Egypt")[1:2], labels=c('Blood feeding', 'Non blood feeding'), name='Behaviour')+
  geom_point(aes(col=as.factor(V4)), size=2.5)+
  geom_highlight(node=121, fill='red', alpha=0.2)+
  geom_highlight(node=70, fill='red', alpha=0.2)+
  geom_cladelabel(node=70, label='Culicidae (mosquitoes)', align=T, color='orange', offset = offsetlab, angle = -90, fontsize = 4, offset.text = 3, hjust='center')+
  geom_cladelabel(node=109, label='Chironomidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=117, label='Ceratopogonidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 2)+
  geom_cladelabel(node=120, label='Psychodidae', align=T, offset = offsetlab, color='orange', fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=106, label='Chaoboridae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_text2(aes(subset= node == 126, label='Outgroups'), hjust=-10.4, size=4, color='grey', vjust = 1.33)+
  geom_cladelabel(node=60, label='Blephariceridae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=39, label='Corethrellidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=40, label='Simuliidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=41, label='Dixidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_nodepoint(aes(subset = node ==70), size=3, color='red')+
  geom_nodepoint(aes(subset = node ==121), size=3, color='red')+
  theme(legend.position = 'bottom')












############################################################################
##################        Boxplot phylogenomics          ###################
############################################################################


# Libraries ----
library(ggplot2)   # CRAN v3.5.0
library(readr)     # CRAN v2.1.5
library(readxl)    # CRAN v1.4.3
library(stringr)   # CRAN v1.5.1
library(MetBrewer) # CRAN v0.2.0
library(cowplot)   # CRAN v1.1.3
library(plyr)      # CRAN v1.8.9

# Import gene counts ----
count=read_delim("~/Desktop/data/CAFE/astral/Gamma_count.tab", delim = "\t", escape_double = FALSE)
OG_interest=c("OG0000128", "OG0001113", "OG0001405", "OG0001557", "OG0002006")
count=count[which(count$FamilyID %in% OG_interest),]
to_keep=str_detect(colnames(count), '^[A-Z]')
count=count[,to_keep]
newnames=sapply(colnames(count[-1]), function(x) unlist(str_split(x, "[<]"))[1])
colnames(count)=c('FamilyID', newnames)

get_table <- function(i){
  df=(t(count[i,]))
  colnames(df)='nb genes'
  sp=rownames(df)
  OG=rep(df[1,1], nrow(df))
  df2=as.data.frame(cbind(df, sp, OG))
  df2=df2[-1,]
  return(df2)
}

OG1=get_table(1)
OG2=get_table(2)
OG3=get_table(3)
OG4=get_table(4)
OG5=get_table(5)
df=rbind(OG1, OG2, OG3, OG4, OG5)


# Import sample info ----
informations=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='R_table')
for (i in 1:length(rownames(df))){
  df[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  df[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(df$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
df$V5=new


# Get mean ----

# OG0001113
for (i in 1:length(rownames(OG1))){
  OG1[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  OG1[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(OG1$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
OG1$V5=new
means_OG0001113_b_n <- aggregate(as.numeric(`nb genes`) ~  V5, OG1, mean)
means_OG0001113_group <- aggregate(as.numeric(`nb genes`) ~  V4, OG1, mean)

# OG0001405
for (i in 1:length(rownames(OG2))){
  OG2[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  OG2[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(OG2$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
OG2$V5=new
means_OG0001405_b_n <- aggregate(as.numeric(`nb genes`) ~  V5, OG2, mean)
means_OG0001405_group <- aggregate(as.numeric(`nb genes`) ~  V4, OG2, mean)

# OG0001557
for (i in 1:length(rownames(OG3))){
  OG3[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  OG3[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(OG3$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
OG3$V5=new
means_OG0001557_b_n <- aggregate(as.numeric(`nb genes`) ~  V5, OG3, mean)
means_OG0001557_group <- aggregate(as.numeric(`nb genes`) ~  V4, OG3, mean)

# OG0002006
for (i in 1:length(rownames(OG4))){
  OG4[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  OG4[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(OG4$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
OG4$V5=new
means_OG0002006_b_n <- aggregate(as.numeric(`nb genes`) ~  V5, OG4, mean)
means_OG0002006_group <- aggregate(as.numeric(`nb genes`) ~  V4, OG4, mean)

# OG0000128
for (i in 1:length(rownames(OG5))){
  OG5[i,4]=informations$boxplot[which(df$sp[i]==informations$tree_name)]
  OG5[i,5]=informations$biting[which(df$sp[i]==informations$tree_name)]
}
new=revalue(as.factor(OG5$V5), c("1"="Blood \n feeding", "0"="Non-blood \n feeding"))
OG5$V5=new
means_OG0000128_b_n <- aggregate(as.numeric(`nb genes`) ~  V5, OG5, mean)
means_OG0000128_group <- aggregate(as.numeric(`nb genes`) ~  V4, OG5, mean)


# Boxplot ----

p1=ggplot(df, aes(x=V4, y=as.numeric(`nb genes`), fill=V4))+
  geom_boxplot()+
  facet_grid(~factor(OG), scales='free_x')+
  scale_fill_manual(values=met.brewer("Lakota", 3))+
  xlab('')+
  ylab('Number of genes')+
  theme_linedraw()+
  theme(legend.position = 'none', axis.text.x = element_text(size=12), axis.title.y = element_text(size=12))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize =  0.5, position=position_dodge(1))
p2=ggplot(df, aes(x=as.factor(V5), y=as.numeric(`nb genes`), fill=as.factor(V5)))+
  geom_boxplot()+
  facet_grid(~factor(OG), scales='free_x')+
  scale_fill_manual(values=met.brewer("Lakota", 3)[c(1,3)])+
  xlab('')+
  ylab('Number of genes')+
  theme_linedraw()+
  theme(legend.position = 'none', axis.text.x = element_text(size=12), axis.title.y = element_text(size=12))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize =  0.5, position=position_dodge(1))
plot_grid(p1,p2, nrow=2, labels = c('A','B'))























############################################################################
########################        Statistics          ########################
############################################################################


# Libraries ----
library(geiger)    # CRAN v2.0.11
library(picante)   # CRAN v1.8.2
library(phytools)  # CRAN v2.1-1
library(readr)     # CRAN v2.1.5
library(readxl)    # CRAN v1.4.3
library(stringr)   # CRAN v1.5.1
library(ggtree)    # Bioconductor v3.8.2
library(plyr)      # CRAN v1.8.9


# Get data ----

## Gene counts
count=read_delim("~/Desktop/data/CAFE/astral/Gamma_count.tab", delim = "\t", escape_double = FALSE)
OG_interest=c("OG0000128", "OG0001113", "OG0001405", "OG0001557", "OG0002006")
count=count[which(count$FamilyID %in% OG_interest),]
to_keep=str_detect(colnames(count), '^[A-Z]')
count=count[,to_keep]
newnames=sapply(colnames(count[-1]), function(x) unlist(str_split(x, "[<]"))[1])
colnames(count)=c('FamilyID', newnames)
count_df=t(count)
colnames(count_df)=count_df[1,]
count_df=count_df[-1,]

## Tree
tree = read.newick('~/Desktop/data/Phylogenomics/r8s/astral_rooted.tre')

## Information 
information=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='R_table')

# Preprocessing ----

## Match tree and data
nb_genes=match.phylo.data(tree, as.data.frame(count_df))$data
tree2=match.phylo.data(tree, as.data.frame(count_df))$phy
name.check(tree2, nb_genes)
behaviour=vector()
for (i in 1:length(tree2$tip.label)){
  behaviour[i]=information$biting[which(tree2$tip.label[i]==information$tree_name)]
}
behaviour_f=revalue(as.factor(behaviour), c("1"="Blood feeding", "0"="Non-blood feeding"))
names(behaviour_f)=tree2$tip.label

group_id=vector()
behaviour=vector()
for (i in 1:length(rownames(count_df))){
  behaviour[i]=information$biting[which(rownames(count_df)[i]==information$tree_name)]
  group_id[i]=information$boxplot[which(rownames(count_df)[i]==information$tree_name)]
}
behaviour_f=revalue(as.factor(behaviour), c("1"="Blood feeding", "0"="Non-blood feeding"))
names(behaviour_f)=rownames(count_df)
names(group_id)=rownames(count_df)



# Run phylogenetic ANOVA ----
nb1=as.vector(nb_genes$OG0001113)
names(nb1)=rownames(nb_genes)
nb2=as.vector(nb_genes$OG0001405)
names(nb2)=rownames(nb_genes)
nb3=as.vector(nb_genes$OG0001557)
names(nb3)=rownames(nb_genes)
nb4=as.vector(nb_genes$OG0002006)
names(nb4)=rownames(nb_genes)
nb5=as.vector(nb_genes$OG0000128)
names(nb5)=rownames(nb_genes)

t1=aov.phylo(nb1~behaviour_f, tree2, nsim = 100) 
t2=aov.phylo(nb2~behaviour_f, tree2, nsim = 100)
t3=aov.phylo(nb3~behaviour_f, tree2, nsim = 100)
t4=aov.phylo(nb4~behaviour_f, tree2, nsim = 100)
t5=aov.phylo(nb5~behaviour_f, tree2, nsim = 100)
#print(attributes(t1)$summary)
#print(attributes(t2)$summary)
#print(attributes(t3)$summary)
#print(attributes(t4)$summary)
#print(attributes(t5)$summary)

# Run ANOVA ----
nb1=as.vector(count_df[,1])
names(nb1)=rownames(count_df)
res_aov <- aov(nb1~behaviour_f + group_id)
summary(res_aov)


nb2=as.vector(count_df[,2])
names(nb2)=rownames(count_df)
res_aov <- aov(nb2~behaviour_f + group_id)
summary(res_aov)

nb3=as.vector(count_df[,3])
names(nb3)=rownames(count_df)
res_aov <- aov(nb3~behaviour_f + group_id)
summary(res_aov)

nb4=as.vector(count_df[,4])
names(nb4)=rownames(count_df)
res_aov <- aov(nb4~behaviour_f + group_id)
summary(res_aov)


nb5=as.vector(count_df[,5])
names(nb5)=rownames(count_df)
res_aov <- aov(nb5~behaviour_f + group_id)
summary(res_aov)


TukeyHSD(res_aov, 'group_id')









############################################################################
###################        Trees phylogenomics          ####################
############################################################################


# Libraries ----
library(readr)     # CRAN v2.1.5
library(ggtree)    # Bioconductor v3.8.2
library(treeio)    # Bioconductor v1.24.3
library(tidyr)     # CRAN v1.3.1
library(stringr)   # CRAN v1.5.1
library(MetBrewer) # CRAN v0.2.0
library(ggplot2)   # CRAN v3.5.0
library(dplyr)     # CRAN v1.1.4
library(readxl)    # CRAN v1.4.3
library(cowplot)   # CRAN v1.1.3


# Working directory ----
setwd('~/Desktop/data/CAFE/astral/')

# Function transpose table ----
table_to_table <- function (tableau, tr){
  tableau=t(tableau)
  colnames(tableau)=tableau[1,]
  tableau=tableau[-1,]
  for (i in 1:length(rownames(tableau))){
    if (unlist(strsplit(rownames(tableau)[i], "[<>]"))[1]!=""){
      rownames(tableau)[i]=unlist(strsplit(rownames(tableau)[i], "[<>]"))[1]
    }
    else{
      rownames(tableau)[i]=unlist(strsplit(rownames(tableau)[i], "[<>]"))[2]
    }
  }
  data.frame(tableau[match(tr$data$label,row.names(tableau)),])
}


# Import and setup the tree ----
tree = read.newick('tree.tre')

for (i in 1:length(tree$tip.label)){
  tree$tip.label[i]=unlist(strsplit(tree$tip.label[i], "[<>]"))[1]
}
for (i in 1:length(tree$node.label)){
  tree$node.label[i]=unlist(strsplit(tree$node.label[i], "[<>]"))[2]
}

informations=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='R_table')

# Preparation visualisation
tr=ggtree(tree, size=0.1)+geom_tiplab(aes(x=x+10), size=1.55)+xlim(c(0,340))#+geom_text(aes(label=node), size=2)#+geom_nodelab()
ggtree(tree)+geom_tiplab(aes(x=x+10), size=1.55)+xlim(c(0,340))+geom_nodelab(size=3, hjust=-0.2)
# Import and setup the changes ----
families_changes=read_delim("~/Desktop/data/CAFE/astral/Gamma_change.tab", delim = "\t", escape_double = FALSE)
families_changes=table_to_table(families_changes, tr)

# Import and setup p-values ----
pval=read_delim("~/Desktop/data/CAFE/astral/Gamma_branch_probabilities.tab", delim = "\t", escape_double = FALSE)
pval=table_to_table(pval, tr)

for (i in 1:64){
  tr$data$label[i]=informations$new_name[which(tr$data$label[i]==informations$tree_name)]
}

tr2=collapse(tr, 70, 'max') 
tr3=collapse(tr2, 121, 'max')

# Import and setup counts of copies ----

#nb_gene=read_delim("~/Desktop/data/CAFE/astral/Gamma_count.tab", delim = "\t", escape_double = FALSE)
#of_interest=c("OG0000128", "OG0001113", "OG0001405", "OG0001557", "OG0002006")
#to_keep=which(nb_gene$FamilyID %in% of_interest)
#nb_gene=nb_gene[to_keep,]
#for (i in 1:length(colnames(nb_gene))){
#  if (unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[1]!=""){
#    colnames(nb_gene)[i]=unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[1]
#  }
#  else{
#    colnames(nb_gene)[i]=unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[2]
#  }
#}
#colnames(nb_gene)=str_replace(colnames(nb_gene), '_', ' ')
#write_csv(nb_gene, '~/Desktop/nb_genes.csv', col_names = T)
#
nb_gene=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='OG')
nb_gene=table_to_table(nb_gene, tr)


# Function visualization duplications and looses history ----
visu_OG_history <- function(OG, p, nb=F){
  position=which(colnames(families_changes)==as.character(OG))
  nbgenes=which(colnames(nb_gene)==as.character(OG))
  pv=which(colnames(pval)==as.character(OG))
  tr3+ggtitle(OG)+
    {if(nb)geom_text(aes(label=nb_gene[,nbgenes]), hjust=-.3, size=1.2)}+
    geom_text(aes(label=if_else(as.numeric(pval[,pv])<p, as.character('*'), '')), vjust=-0.25, hjust=-0.35, color='black', size=3)+
    geom_text(aes(label=if_else(as.numeric(families_changes[,position])!=0, families_changes[,position],''), color=ifelse(as.numeric(families_changes[,position])>0, 'positif', 'negatif')),vjust=-0.5, size=1.8)+
    scale_color_manual(values=c(met.brewer('Austria', n=3, type = 'discrete')[1], met.brewer('Austria', n=3, type = 'discrete')[3]))+ 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, face='bold', size=12))
}



p1=visu_OG_history(OG='OG0000128', p=0.05, nb=T)
p2=visu_OG_history(OG='OG0001113', p=0.05, nb=T)
p3=visu_OG_history(OG='OG0001405', p=0.05, nb=T)
p4=visu_OG_history(OG='OG0001557', p=0.05, nb=T)
p5=visu_OG_history(OG='OG0002006', p=0.05, nb=T)
plot_grid(p1,p2,p3,p4,p5, labels = c('A', 'B', 'C', 'D', 'E'), nrow = 3, ncol = 2)











# Working directory ----
setwd('~/Desktop/data/CAFE/busco/')

# Function transpose table ----
table_to_table <- function (tableau, tr){
  tableau=t(tableau)
  colnames(tableau)=tableau[1,]
  tableau=tableau[-1,]
  for (i in 1:length(rownames(tableau))){
    if (unlist(strsplit(rownames(tableau)[i], "[<>]"))[1]!=""){
      rownames(tableau)[i]=unlist(strsplit(rownames(tableau)[i], "[<>]"))[1]
    }
    else{
      rownames(tableau)[i]=unlist(strsplit(rownames(tableau)[i], "[<>]"))[2]
    }
  }
  data.frame(tableau[match(tr$data$label,row.names(tableau)),])
}


# Import and setup the tree ----
tree = read.newick('tree.tre')

for (i in 1:length(tree$tip.label)){
  tree$tip.label[i]=unlist(strsplit(tree$tip.label[i], "[<>]"))[1]
}
for (i in 1:length(tree$node.label)){
  tree$node.label[i]=unlist(strsplit(tree$node.label[i], "[<>]"))[2]
}

informations=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='R_table')

# Preparation visualisation
tr=ggtree(tree, size=0.1)+geom_tiplab(aes(x=x+10), size=1.55)+xlim(c(0,340))#+geom_text(aes(label=node), size=2)#+geom_nodelab()
ggtree(tree)+geom_tiplab(aes(x=x+10), size=1.55)+xlim(c(0,340))+geom_nodelab(size=3, hjust=-0.2)
# Import and setup the changes ----
families_changes=read_delim("~/Desktop/data/CAFE/busco/Gamma_change.tab", delim = "\t", escape_double = FALSE)
families_changes=table_to_table(families_changes, tr)

# Import and setup p-values ----
pval=read_delim("~/Desktop/data/CAFE/busco/Gamma_branch_probabilities.tab", delim = "\t", escape_double = FALSE)
pval=table_to_table(pval, tr)

for (i in 1:64){
  tr$data$label[i]=informations$new_name[which(tr$data$label[i]==informations$tree_name)]
}

tr2=collapse(tr, 72, 'max') 
tr3=collapse(tr2, 121, 'max')

# Import and setup counts of copies ----

#nb_gene=read_delim("~/Desktop/data/CAFE/busco/Gamma_count.tab", delim = "\t", escape_double = FALSE)
#of_interest=c("OG0000128", "OG0001113", "OG0001405", "OG0001557", "OG0002006")
#to_keep=which(nb_gene$FamilyID %in% of_interest)
#nb_gene=nb_gene[to_keep,]
#for (i in 1:length(colnames(nb_gene))){
#  if (unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[1]!=""){
#    colnames(nb_gene)[i]=unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[1]
#  }
#  else{
#    colnames(nb_gene)[i]=unlist(strsplit(colnames(nb_gene)[i], "[<>]"))[2]
#  }
#}
#colnames(nb_gene)=str_replace(colnames(nb_gene), '_', ' ')
#write_csv(nb_gene, '~/Desktop/nb_genes.csv', col_names = T)
#
nb_gene=read_xlsx('~/Desktop/data/Phylogenomics/0_final_list.xlsx', sheet='OG_busco')
nb_gene=table_to_table(nb_gene, tr)


# Function visualization duplications and looses history ----
visu_OG_history <- function(OG, p, nb=F){
  position=which(colnames(families_changes)==as.character(OG))
  nbgenes=which(colnames(nb_gene)==as.character(OG))
  pv=which(colnames(pval)==as.character(OG))
  tr3+ggtitle(OG)+
    {if(nb)geom_text(aes(label=nb_gene[,nbgenes]), hjust=-.3, size=1.2)}+
    geom_text(aes(label=if_else(as.numeric(pval[,pv])<p, as.character('*'), '')), vjust=-0.25, hjust=-0.35, color='black', size=3)+
    geom_text(aes(label=if_else(as.numeric(families_changes[,position])!=0, families_changes[,position],''), color=ifelse(as.numeric(families_changes[,position])>0, 'positif', 'negatif')),vjust=-0.5, size=1.8)+
    scale_color_manual(values=c(met.brewer('Austria', n=3, type = 'discrete')[1], met.brewer('Austria', n=3, type = 'discrete')[3]))+ 
    theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, face='bold', size=12))
}



p1=visu_OG_history(OG='OG0000128', p=0.05, nb=T)
p2=visu_OG_history(OG='OG0001113', p=0.05, nb=T)
p3=visu_OG_history(OG='OG0001405', p=0.05, nb=T)
p4=visu_OG_history(OG='OG0001557', p=0.05, nb=T)
p5=visu_OG_history(OG='OG0002006', p=0.05, nb=T)
plot_grid(p1,p2,p3,p4,p5, labels = c('A', 'B', 'C', 'D', 'E'), nrow = 3, ncol = 2)


