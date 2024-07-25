#### Phylogenetic Signal ####

v_BEASTtree <- read.tree("Vellozia_Signal_Phylogeny_BEAST.tre")  ## Phylogenetic trees
v_RAxMLtree <- read.tree("Vellozia_Signal_Phylogeny_RAxML.tre")

## Variables we used were percentage of habitat lost based on our SDMs, and the
## area of suitable habitat based on our SDMs

hab_loss <- read.csv("Vellozia_Species_vs_Percent_Habitat_Lost.csv",
                     row.names = 1)
habitat_loss <- setNames(hab_loss$percentage_of_habitat_lost,
                   rownames(hab_loss))

area <- read.csv("Vellozia_Species_vs_Present_Habitat.csv",
                 row.names = 1)
present_area <- setNames(log(area$present_habitat_.km2.),
                         rownames(area))

#### Calculate Phylogenetic Signal with Blomberg's K ####

hablossK_BEAST <- phylosig(v_BEASTtree,
                 habitat_loss,test=TRUE,nsim=10000)

hablossK_RAxML <- phylosig(v_RAxMLtree,
                           habitat_loss,test=TRUE,nsim=10000)

areaK_BEAST <- phylosig(v_BEASTtree,
                  present_area,test=TRUE,nsim=10000)

areaK_RAxML <- phylosig(v_RAxMLtree,
                        present_area,test=TRUE,nsim=10000)


BEAST_results_K_habloss <- data.frame(hablossK_BEAST$K, hablossK_BEAST$P)
colnames(BEAST_results) <- c("K", "P-value")

RAxML_results_K_habloss <- data.frame(hablossK_RAxML$K, hablossK_RAxML$P)
colnames(RAxML_results) <- c("K", "P-value")

BEAST_results_K_area <- data.frame(areaK_BEAST$K, areaK_BEAST$P)
colnames(BEAST_results) <- c("K", "P-value")

RAxML_results_K_area <- data.frame(areaK_RAxML$K, areaK_RAxML$P)
colnames(RAxML_results) <- c("K", "P-value")


hablossK <- rbind(BEAST_results_K_habloss, RAxML_results_K_habloss)
rownames(hablossK) <- c("BEAST", "RAxML")

areaK <- rbind(BEAST_results_K_area, RAxML_results_K_area)
rownames(areaK) <- c("BEAST", "RAxML")

#### Calculate Phylogenetic Signal with Pagel's Lambda ####

hablossL_BEAST <- phylosig(v_BEASTtree,
                           habitat_loss,test=TRUE, method = "lambda")

hablossL_RAxML <- phylosig(v_RAxMLtree,
                           habitat_loss,test=TRUE, method = "lambda")

areaL_BEAST <- phylosig(v_BEASTtree,
                        present_area,test=TRUE, method = "lambda")

areaL_RAxML <- phylosig(v_RAxMLtree,
                        present_area,test=TRUE, method = "lambda")


BEAST_results_lambda_habloss <- data.frame(hablossL_BEAST$lambda, hablossL_BEAST$P)
colnames(BEAST_results) <- c("lambda", "P-value")

RAxML_results_lambda_habloss <- data.frame(hablossL_RAxML$lambda, hablossL_RAxML$P)
colnames(RAxML_results) <- c("lambda", "P-value")

BEAST_results_lambda_area <- data.frame(areaL_BEAST$lambda, areaL_BEAST$P)
colnames(BEAST_results) <- c("lambda", "P-value")

RAxML_results_lambda_area <- data.frame(areaL_RAxML$lambda, areaL_RAxML$P)
colnames(RAxML_results) <- c("lambda", "P-value")


hablossL <- rbind(BEAST_results_lambda_habloss, RAxML_results_lambda_habloss)
rownames(hablossL) <- c("BEAST", "RAxML")

areaL <- rbind(BEAST_results_lambda_area, RAxML_results_lambda_area)
rownames(areaL) <- c("BEAST", "RAxML")

#### Compile Results ####

area_results <- cbind(areaK, areaL)
habloss_results <- cbind(hablossK, hablossL)

write.csv(area_results, "Vellozia_Present_Area_Phylo_Results.csv")
write.csv(habloss_results, "Vellozia_Percent_Habloss_Phylo_Results.csv")
