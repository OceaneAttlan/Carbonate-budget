
# --------------------------- README ---------------------------
# This script analyses data from the study:
# "Diverging carbonate budgets following tropicalisation of temperate reefs"

# The analysis is divided into three parts:
# 1. Coral calcification rates (Figure 3)
# 2. Community % cover by organism (Figure 2)
# 3. Carbonate production and bioerosion balance (Figure 4)

# --------------------------- Load Required Libraries ---------------------------
library(readxl)        # For reading Excel files
library(ggplot2)         # For plotting
library(DHARMa)          # For diagnostic tests for (G)LMMs
library(glmmTMB)         # For fitting GLMMs, especially zero-inflated models
library(dplyr)           # For data wrangling
library(tidyverse)       # Collection of packages (includes dplyr, ggplot2, etc.)
library(vegan) 
library(mvabund)# For multivariate statistics (PERMANOVA, SIMPER, etc.)
library(emmeans)

# --------------------------- Part 1: Coral Calcification Rates ---------------------------

# Read coral calcification data
calcification_rate <- read_excel("C:/Users/23692811/OneDrive - UWA/Writing/Manuscripts/Coral/Submission/Data_submission_20250623.xlsx", sheet="Figure_3_Coral_Calcification")

# Create Figure 3: Coral calcification rates by species
Figure3 <- ggplot(calcification_rate, aes(x = factor(Species, level=rev(c("Coelestra aspera", "Turbinaria sp", "Acropora sp", "Acropora yongei"))), 
                                          y = Gnet_g.cm2.year, color = Morphotype, group = Species)) + 
  geom_point(size = 10, alpha = 0.2, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) + 
  stat_summary(fun = mean, geom = "point", size = 14, position = position_dodge(0)) + 
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.3, size = 1.5, position = position_dodge(0)) +
  scale_color_manual(values = c("#FF9900","#990000", "#66CCFF")) +
  scale_x_discrete(labels = c(
    expression(italic("Acropora yongei")),
    expression(italic("Acropora")~"sp"),
    expression(italic("Turbinaria")~"sp"),
    expression(italic("Coelestra aspera"))
  )) +
  labs(x = NULL, y = expression ("Calcification rate (g cm"^"-2"~"yr"^"-1"~")")) +
  theme(axis.text.x=element_text(size=25,vjust=1,hjust=0.5,margin=margin(t=10),color = "black"))+
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(size=30,margin=margin(r=20),color = "black"),
        axis.text=element_text(size=25),
        axis.text.y=element_text(margin=margin(r=10),color = "black"),
        strip.text.x = element_text(size = 15,vjust=1,margin=margin(b=10)),
        panel.background = element_rect(fill='transparent'),
        legend.text = element_text(size=25),
        legend.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        # Remove grey background
        strip.background = element_rect(colour = "white", fill = "white", size=20),  # Customize strip background
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))  # Remove grid lines if desired

# Add significance annotations
Figure3 <- Figure3 +
  annotate("text", x = 1, y = max(calcification_rate$Gnet_g.cm2.year, na.rm = TRUE) * 0.60, label = "ab", size = 10) +
  annotate("text", x = 2, y = max(calcification_rate$Gnet_g.cm2.year, na.rm = TRUE) * 1.05, label = "a", size = 10) +
  annotate("text", x = 3, y = max(calcification_rate$Gnet_g.cm2.year, na.rm = TRUE) * 0.8, label = "b", size = 10) +
  annotate("text", x = 4, y = max(calcification_rate$Gnet_g.cm2.year, na.rm = TRUE) * 0.8, label = "ab", size = 10)

# Display the figure
Figure3

# Statistical test: one-way ANOVA for species effect on calcification rate
m0.lm <- lm(Gnet_g.cm2.year ~ Species, data = calcification_rate)
anova(m0.lm)

# Diagnostic plots
testDispersion(m0.lm)
simulationOutput <- simulateResiduals(m0.lm)
plot(simulationOutput)
plotResiduals(simulationOutput, form = m0.lm$Species)

# Model fit diagnostics
AIC(logLik(m0.lm))
plot(m0.lm)

# Post-hoc comparisons using lsmeans
multcomp::cld(lsmeans(m0.lm,  ~ Species), Letters=letters)

# --------------------------- Part 2: Community Composition (% cover) ---------------------------

# Read coral community data
CoralComm_fig <- read_excel("C:/Users/23692811/OneDrive - UWA/Writing/Manuscripts/Coral/Submission/Data_submission_20250623.xlsx", sheet="CoralComm_fig")
Scenario_col <- c("goldenrod1","darkblue","hotpink","gold4")

# Summarize % cover per transect per organism (CCA, Articulated macroalgae, corals)
CoralComm_cover_orga <- CoralComm_fig %>%
  group_by(Scenario, Site, Transect, Organism) %>%
  summarize(Orga_sum_perc_cover = sum(Perc_cover)) %>%
  ungroup()

CoralComm_cover_tot= CoralComm_fig %>%
  group_by(Scenario,Site, Transect)%>%
  summarize(Orga_sum_perc_cover = sum(Perc_cover))%>%
  ungroup()

# Create Figure 4: Community % cover
Figure4B <- ggplot(CoralComm_cover_orga, 
                   aes(x=Site,
                       y=Orga_sum_perc_cover,color=Scenario, shape=Organism, group= Organism)) + 
  facet_wrap(.~Scenario, scales="free",ncol=4)+
  geom_point(size =2, alpha=0.2,position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0)) + 
  stat_summary(aes(group = Organism), fun = mean, geom = "point", size =4,position = position_dodge(0.75)) +
  stat_summary(aes(group = Organism),geom = "errorbar",  fun.args = list(mult = 1), 
               width = 0.3,size=1, position = position_dodge(0.75)) +
  scale_color_manual(values = Scenario_col)+
  labs(x=NULL, y = "Calcifier cover (%)")+
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(size=15,angle=90,margin=margin(r=20)),
        axis.text=element_text(size=10),
        axis.text.y=element_text(margin=margin(r=0)),
        panel.background = element_rect(fill='transparent'),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

Figure4B

Figure4B <- Figure4B +
  guides(
    color = guide_legend(order = 1, title = "Scenario"),
    shape = guide_legend(order = 2, title = "Organism")
  ) 

Figure4A <- ggplot(CoralComm_cover_tot, 
                              aes(x=Scenario,
                                  y=Orga_sum_perc_cover,
                                  fill=Scenario)) + 
  geom_boxplot(alpha = 0.6, width = 0.6, outlier.shape = NA, color = "black") + 
  geom_jitter(aes(color = Scenario), 
              size = 3, width = 0.2, alpha = 0.8, show.legend = FALSE) +
  labs(x = NULL, y = "Total calcifier cover (%)") +
  scale_fill_manual(values=Scenario_col) +
  scale_color_manual(values=Scenario_col) +
  theme(axis.text.x=element_blank())+
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(size=15,angle=90,margin=margin(r=20)),
        axis.text=element_text(size=10),
        axis.text.y=element_text(margin=margin(r=0)),
        strip.text.x = element_text(size = 10,margin=margin(b=10)),
        panel.background = element_rect(fill='transparent'),
        theme(legend.position = "none"),
        strip.background = element_rect(colour = "black", fill = "white"),  
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) 

Figure4A
Figure4A <- Figure4A + theme(legend.position = "none")

Figure4 <-Figure4A/ Figure4B  + 
  plot_layout(heights = c(1, 2), guides = "collect")  


# Display figure
Figure4

# Statistical model: GLMM with zero-inflated Gamma distribution
zfit_gamma <- glmmTMB(Orga_sum_perc_cover ~ Organism * Scenario + (1 | Site),
                      zi = ~1,
                      family = ziGamma(link = "log"),
                      data = CoralComm_cover_orga)

# Model diagnostics
testDispersion(zfit_gamma)
simulationOutput <- simulateResiduals(zfit_gamma)
plot(simulationOutput)
plotResiduals(simulationOutput, form = zfit_gamma$Organism)

# Post-hoc comparisons
multcomp::cld(lsmeans(zfit_gamma, ~ Organism:Scenario), Letters = letters)

# PERMANOVA: community structure differences by scenario
community_matrix <- CoralComm_fig %>%
  select(Scenario, Morph, Site, Transect, Perc_cover) %>%
  pivot_wider(names_from = Morph, values_from = Perc_cover)

response_variable <- as.matrix(community_matrix[, 4:11])
Morph.mat <- sqrt(response_variable)                      # sqrt transform
Morph.dist <- vegdist(Morph.mat, method = 'bray')         # Bray-Curtis dissimilarity

Morph.div <- adonis2(Morph.dist ~ Scenario, data = community_matrix, permutations = 9999)

# Beta-dispersion (homogeneity of dispersion test)
dispersion <- betadisper(Morph.dist, group = community_matrix$Scenario)
permutest(dispersion)

# SIMPER analysis
simper_result <- simper(response_variable, group = community_matrix$Scenario)
summary(simper_result)

# --------------------------- Part 3: Carbonate Production and Bioerosion ---------------------------

# Read summary of bioerosion 
Summary.bioerosion.orga.prod <- read_excel("C:/Users/23692811/OneDrive - UWA/Writing/Manuscripts/Coral/Submission/Data_submission_20250623.xlsx", sheet="Summary.bioerosion.orga.prod")

# Ensure numeric format
CoralComm_fig$`Production_(kg/m2/yr)` <- as.numeric(CoralComm_fig$`Production_(kg/m2/yr)`)

# Step 1: Sum production per transect per Morph
CoralComm_prod_morph <- CoralComm_fig %>%
  group_by(Scenario, Site, Transect, Organism, Morph) %>%
  summarize(tot_prod_morph = sum(`Production_(kg/m2/yr)`), .groups = "drop")

# Step 2: Avg production per site per Morph
Summary.calcifying.morph.prod <- CoralComm_prod_morph %>%
  group_by(Scenario, Site, Organism, Morph) %>%
  summarize(Prod_morph = mean(tot_prod_morph), se = sd(tot_prod_morph)/sqrt(n()), .groups = "drop")

# Step 3-4: Repeat for Organism-level production
CoralComm_prod_orga <- CoralComm_fig %>%
  group_by(Scenario, Site, Transect, Organism) %>%
  summarize(tot_prod_orga = sum(`Production_(kg/m2/yr)`), .groups = "drop")

Summary.calcifying.orga.prod <- CoralComm_prod_orga %>%
  group_by(Scenario, Site, Organism) %>%
  summarize(Prod_orga = mean(tot_prod_orga), se = sd(tot_prod_orga)/sqrt(n()), .groups = "drop")

# Step 5-6: Total production per transect & site
CoralComm_prod_tot <- CoralComm_fig %>%
  group_by(Scenario, Site, Transect) %>%
  summarize(tot_prod_transect = sum(`Production_(kg/m2/yr)`), .groups = "drop")

Summary.calcifying.total.prod <- CoralComm_prod_tot %>%
  group_by(Scenario, Site) %>%
  summarize(Prod_site = mean(tot_prod_transect), se = sd(tot_prod_transect)/sqrt(n()), .groups = "drop")

# Step 7-9: Group by Scenario or Scenario x Morph/Organism
Summary.calcifying.total.prod_scenario <- Summary.calcifying.total.prod %>%
  group_by(Scenario) %>%
  summarize(Prod_scenario = mean(Prod_site), se = sd(Prod_site)/sqrt(n()), .groups = "drop")

Summary.calcifying.total.prod_scenario_morph=Summary.calcifying.morph.prod%>%
  group_by(Scenario,Morph)%>%
  summarize(Prod_scenario_orga = mean(Prod_morph),
            sd=sd(Prod_morph),
            n=n(),
            se = sd(Prod_morph) / sqrt(n()),
            .groups="drop")

Summary.calcifying.total.prod_scenario_orga=Summary.calcifying.orga.prod%>%
  group_by(Scenario,Organism)%>%
  summarize(Prod_scenario_orga = mean(Prod_orga),
            sd=sd(Prod_orga),
            n=n(),
            se = sd(Prod_orga) / sqrt(n()),
            .groups="drop")



# Figure 4: Carbonate budget (production - erosion)
colors <- c("#d24955ff", "#d24955ff", "#5358b0ff", "#5358b0ff", "#5358b0ff", 
            "#5358b0ff", "#5358b0ff", "#5358b0ff", "#e29939ff", "#e29939ff")

Summary.calcifying.morph.prod$Morph <- factor(Summary.calcifying.morph.prod$Morph, 
                                              levels = c("CCA", "Articulated macroalgae", 
                                                         setdiff(unique(Summary.calcifying.morph.prod$Morph), c("CCA", "Articulated macroalgae")))
)

Figure5=ggplot()+
  
  #Production:
  geom_bar(data = Summary.calcifying.morph.prod, stat="identity",
           aes(y=Prod_morph,
               #x=forcats::fct_reorder(Site, Prod_morph, .fun=mean, .desc=T),
               x=Site,
               fill=Morph,
               color=Morph), 
           position="stack", alpha=0.7)+
  ###Errorbars
  geom_errorbar(data = Summary.calcifying.total.prod,  # errorbars from the total 
                aes(x = Site, y=Prod_site,  ymax=Prod_site+se, ymin=Prod_site-se),
                color = "black",width=.05, size=0.75)+
  #Bioerosion:
  geom_bar(data = Summary.bioerosion.orga.prod, stat="identity",
           aes(y=-Prod_orga_site, 
               x=Site,
               fill=Organism,
               color=Organism),
           position="stack", alpha=0.7)+
  ###Errorbars
  geom_errorbar(data = Summary.bioerosion.orga.prod,  # errorbars from the total 
                aes(x = Site, y=-Prod_orga_site,  ymax=-Prod_orga_site+se, ymin=-Prod_orga_site-se),
                color = "black",width=.05, size=0.75)+
  facet_wrap(.~factor(Scenario, levels = c("Temperate kelp", "Transition stage", "Warmer-affinity seaweeds", "Turf/coral")), 
             ncol = 4, scales = "free_x") +
  ##Aesthetics
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  labs(x=NULL, y = "Carbonate production (kg m"^2~"yr"^-1~")")+
  theme(axis.text.x=element_text(angle=90,size=30,vjust=0.5,hjust=1,margin=margin(t=10)))+
  theme(axis.title.x=element_text(),
        axis.title.y=element_text(size=30,angle=90,margin=margin(r=20)),
        axis.text=element_text(size=30),
        axis.text.y=element_text(margin=margin(r=10)),
        strip.text.x = element_text(size = 27,margin=margin(b=10)),
        panel.background = element_rect(fill='transparent'),
        legend.text = element_text(size=30),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add frame around each facet
        # Remove grey background
        strip.background = element_rect(colour = "black", fill = "white"),  # Customize strip background
        panel.grid = element_blank())  # Remove grid lines if desired

# Display figure
Figure5

# Statistical model: GLMM with zero-inflated Gamma distribution
# for total carbonate production
zfit_gamma <- glmmTMB(tot_prod_transect ~ Scenario + (1 | Site),
                      family = ziGamma(link = "log"),
                      data = CoralComm_prod_tot)

# Model checks and inference
testDispersion(zfit_gamma)
simulationOutput <- simulateResiduals(zfit_gamma)
plot(simulationOutput)
plotResiduals(simulationOutput, form = zfit_gamma$Scenario)
Anova(zfit_gamma, type = "III")

# Post-hoc comparisons
multcomp::cld(lsmeans(zfit_gamma,  ~ Scenario), Letters = letters)
lsm <- emmeans(zfit_gamma, ~ Scenario)
pairwise_comparisons <- pairs(lsm, adjust = "tukey")
summary(pairwise_comparisons)


# Statistical community analyses
# Data Preparation
# 

# Reshape community data to wide format
community_matrix_prod <- CoralComm_prod_morph %>%
  select(Scenario, Morph, Site, Transect, tot_prod_morph) %>%
  pivot_wider(names_from = Morph, values_from = tot_prod_morph)

# Extract response matrix (morph productivity)
response_variable_prod <- as.matrix(community_matrix_prod[, 4:11])

# Apply square-root transformation
Morph.mat <- sqrt(response_variable_prod)


# PERMANOVA Analysis

# Calculate Bray-Curtis dissimilarity
Morph.dist <- vegdist(Morph.mat, method = "bray")

# Run PERMANOVA to test for differences between scenarios
permanova_result <- adonis2(
  Morph.dist ~ as.factor(community_matrix_prod$Scenario),
  data = community_matrix_prod,
  permutations = 9999
)
print(permanova_result)

# Homogeneity of Dispersion

# Test for homogeneity of multivariate dispersion
dispersion <- betadisper(Morph.dist, group = community_matrix_prod$Scenario)
dispersion_test <- permutest(dispersion)
print(dispersion_test)

# mvabund: Model-Based Analysis


# Convert to mvabund object
Morph_mvabund <- mvabund(community_matrix_prod[, 4:11])

# Plot abundance distributions
par(mar = c(2, 10, 2, 2))
boxplot(Morph_mvabund, horizontal = TRUE, las = 2, main = "Morph Abundance")

# Mean-variance relationship
meanvar.plot(Morph_mvabund)

# Visualize abundance by scenario
plot(Morph_mvabund ~ community_matrix_prod$Scenario, cex.axis = 0.8, cex = 0.8)

# Fit GLMs (Poisson and Negative Binomial)
mod_poisson <- manyglm(Morph_mvabund ~ community_matrix_prod$Scenario, family = "poisson")
mod_negbin <- manyglm(Morph_mvabund ~ community_matrix_prod$Scenario, family = "negative.binomial")

# Model diagnostics
plot(mod_poisson)
plot(mod_negbin)

# Model comparison and significance tests
anova_result <- anova(mod_negbin)
anova_uni <- anova(mod_negbin, p.uni = "adjusted")
print(anova_result)
print(anova_uni)

# SIMPER Analysis

# SIMPER: Identify morphs contributing to differences between scenarios
simper_result <- simper(response_variable_prod, group = community_matrix_prod$Scenario)

# Summary of contributions
summary(simper_result)

# Full SIMPER output
print(simper_result)

