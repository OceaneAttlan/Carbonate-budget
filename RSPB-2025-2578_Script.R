
# --------------------------- README ---------------------------
# This script analyses data from the study:
# "Diverging carbonate budgets following tropicalisation of temperate reefs"

# The analysis is divided into three parts:
# 1. Tropicalisation states (Figure 2)
# 2. Community % cover by organism (Figure 3)
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
library(patchwork)

# --------------------------- Part 1: Tropicalisation states ---------------------------

# Load Data

# Read Excel data from the 'Figure_2_habitat' sheet
habitat <- read_excel("path[..]/Data_RSPB-2025-2578.xlsx", 
                  sheet = "Figure_2_habitat")

# Calculate Percent Cover per Functional Group per Site 

# Create a wide-format matrix with summed counts per site
habitat_matrix <- habitat %>%
  group_by(Site, Functional_Group) %>%
  summarize(Count = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = Site, values_from = Count, values_fill = 0)

# Store original column names
original_colnames <- colnames(habitat_matrix)

# Calculate percent cover for each functional group at each site
habitat_matrix_percent <- habitat_matrix %>%
  mutate(across(-Functional_Group, ~ . / sum(.) * 100, .names = "Percent_{.col}")) %>%
  select(Functional_Group, starts_with("Percent_"))

# Restore site names as column headers
colnames(habitat_matrix_percent) <- original_colnames

# Convert to long format for plotting
habitat_long <- habitat_matrix_percent %>%
  pivot_longer(cols = -Functional_Group, names_to = "Site", values_to = "Percentage")

# Bar Plot: Functional Group Composition by Site 

# Define functional group color palette
functional_group_colors <- c(
  "Articulated macroalgae"   = "cadetblue3",
  "CCA"                      = "darkcyan",
  "Cooler-affinity seaweed" = "goldenrod1",
  "Hard coral"               = "hotpink",
  "Other invertebrates"      = "gray",
  "Other primary producers"  = "darkolivegreen4",
  "Sand"                     = "burlywood",
  "Sargassum sp"             = "gold4",
  "Turf"                     = "darkblue"
)

# Set factor levels for functional group and site order
habitat_long <- habitat_long %>%
  mutate(
    Functional_Group = factor(Functional_Group, levels = names(functional_group_colors)),
    Site = factor(Site, levels = rev(c("KAL1", "KAL2", "KAL3", 
                                       "PGC1", "PGC2", "PGC3",
                                       "PGK1", "PGK2", "PGK3", 
                                       "PER1", "PER2", "PER3")))
  )

# Create the ggplot barplot
Figure2 <- ggplot(habitat_long, aes(x = Site, y = Percentage, fill = Functional_Group)) +
  geom_bar(stat = "identity", width = 0.95) +
  coord_flip() +
  scale_fill_manual(values = functional_group_colors) +
  labs(x = "Site", y = "Cover (%)") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(hjust = 1),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.size = unit(1, "cm")
  ) +
  guides(fill = guide_legend(nrow = 3))

# Display the figure
Figure2



# --------------------------- Part 2: Community Composition (% cover) ---------------------------

# Read coral community data
CoralComm_fig <- read_excel("path[..]/Data_RSPB-2025-2578.xlsx", sheet="CoralComm_fig")
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

# Create Figure 3: Community % cover
Figure3B <- ggplot(CoralComm_cover_orga, 
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
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

Figure3B

Figure3B <- Figure3B +
  guides(
    color = guide_legend(order = 1, title = "Scenario"),
    shape = guide_legend(order = 2, title = "Organism")
  ) 

Figure3A <- ggplot(CoralComm_cover_tot, 
                              aes(x=Scenario,
                                  y=Orga_sum_perc_cover,
                                  fill=Scenario)) + 
  facet_grid(~Scenario,scales="free")+
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

Figure3A
Figure3A <- Figure3A + theme(legend.position = "none")

Figure3 <-Figure3A/ Figure3B  + 
  plot_layout(heights = c(1, 2), guides = "collect")  


# Display figure
Figure3

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
Summary.bioerosion.orga.prod <- read_excel("path[..]/Data_RSPB-2025-2578.xlsx", sheet="Summary.bioerosion.orga.prod")

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

Figure4=ggplot()+
  
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
Figure4

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

