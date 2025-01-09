###########################################################################
##################### Ancestral state reconstruction #####################
###########################################################################


# Libraries ----
library(readr)     # CRAN v2.1.4
library(ggtree)    # Bioconductor v3.8.2
library(treeio)    # Bioconductor v1.24.3
library(tidyr)     # CRAN v1.3.0
library(stringr)   # CRAN v1.5.1
library(MetBrewer) # CRAN v0.2.0
library(ggplot2)   # CRAN v3.4.4
library(dplyr)     # CRAN v1.1.3
library(readxl)    # CRAN v1.4.3
library(RRphylo)   # CRAN v2.8.0
library(readxl)    # CRAN v1.4.3
library(phytools)  # CRAN v2.0-3
library(geiger)    # CRAN v2.0.11
library(ggtext)


# Working directory ----
setwd('~/Desktop/data/CAFE/astral/')

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
p=ggtree(trait_astral$tree)+geom_tiplab(aes(x=x+3), fontface = "italic")+xlim(c(0,450))
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


for_pie=end[65:126,1:3]
for_pie[,2:3] <- sapply(for_pie[,2:3], as.numeric)
colnames(for_pie)=c('node', 'biting', 'non-biting')
for_pie$node=rownames(end)[65:126]
pies <- nodepie(for_pie, cols = 2:3)
pies <- lapply(pies, function(g) g+scale_fill_manual(values=met.brewer("Egypt")[1:2], labels=c('Blood feeding', 'Non blood feeding'), name='Behaviour'))

offsetlab=125

p2 + geom_inset(pies, width = .05, height = .05) +
  geom_tippoint(aes(col=as.factor(V4)), size=2)+
  scale_color_manual(values=met.brewer("Egypt")[1:2], labels=c('Blood feeding', 'Non blood feeding'), name='Behaviour')+
  geom_highlight(node=121, fill='red', alpha=0.2)+
  geom_highlight(node=70, fill='red', alpha=0.2)+
  geom_cladelabel(node=70, label='Culicidae (mosquitoes)', align=T, color='orange', offset = offsetlab, angle = -90, fontsize = 4, offset.text = 3, hjust='center')+
  geom_cladelabel(node=109, label='Chironomidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=117, label='Ceratopogonidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 2)+
  geom_cladelabel(node=120, label='Psychodidae', align=T, offset = offsetlab, color='orange', fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=106, label='Chaoboridae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_text2(aes(subset= node == 126, label='Outgroups'), hjust=-7.3, size=4, color='grey', vjust = 1.33)+
  geom_cladelabel(node=60, label='Blephariceridae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=39, label='Corethrellidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=40, label='Simuliidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=41, label='Dixidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3)+
  theme(legend.position = 'bottom')








### BUSCO tree ----


# Working directory ----
setwd('~/Desktop/data/CAFE/busco/')

# Import and setup the tree ----
tree = read.newick('~/Desktop/data/CAFE/busco/tree.tre')

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
p=ggtree(trait_astral$tree)+geom_tiplab(aes(x=x+3), fontface = "italic")+xlim(c(0,450))
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


for_pie=end[65:126,1:3]
for_pie[,2:3] <- sapply(for_pie[,2:3], as.numeric)
colnames(for_pie)=c('node', 'biting', 'non-biting')
for_pie$node=rownames(end)[65:126]
pies <- nodepie(for_pie, cols = 2:3)
pies <- lapply(pies, function(g) g+scale_fill_manual(values=met.brewer("Egypt")[1:2], labels=c('Blood feeding', 'Non blood feeding'), name='Behaviour'))

offsetlab=125

p2 + geom_inset(pies, width = .05, height = .05) +
  geom_tippoint(aes(col=as.factor(V4)), size=2)+
  scale_color_manual(values=met.brewer("Egypt")[1:2], labels=c('Blood feeding', 'Non blood feeding'), name='Behaviour')+
  geom_highlight(node=121, fill='red', alpha=0.2)+
  geom_highlight(node=72, fill='red', alpha=0.2)+
  geom_cladelabel(node=72, label='Culicidae (mosquitoes)', align=T, color='orange', offset = offsetlab, angle = -90, fontsize = 4, offset.text = 3, hjust='center')+
  geom_cladelabel(node=109, label='Chironomidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=117, label='Ceratopogonidae', align=T, offset = offsetlab, fontsize = 4, offset.text = 2)+
  geom_cladelabel(node=120, label='Psychodidae', align=T, offset = offsetlab, color='orange', fontsize = 4, offset.text = 3)+
  geom_cladelabel(node=107, label='Chaoboridae', align=T, offset = offsetlab, fontsize = 4, offset.text = 3)+
  geom_text2(aes(subset= node == 126, label='Outgroups'), hjust=-7.3, size=4, color='grey', vjust = 1.33)+
  geom_cladelabel(node=60, label='Blephariceridae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=39, label='Corethrellidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=41, label='Simuliidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  geom_cladelabel(node=40, label='Dixidae', align=T, offset = offsetlab, fontsize = 3, offset.text = 2)+
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3)+
  theme(legend.position = 'bottom')


