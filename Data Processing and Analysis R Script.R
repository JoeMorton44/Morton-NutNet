###Nutrient Network Manuscript Code

#Load in basic packages
require(stringr)
require(plyr)
require(dplyr)
require(ggplot2)

###====Loading in and formatting initial data====

#This loads in the Nutrient Network species level cover data,  calculates the 3-year mean percent cover values for each species and incorporates our GS values

#Load in cover data - the percent cover of each species found on each site, in the pre-treatment year and the last 3 years
cover <- read.csv( "27-site-cover-data.csv")

#Calculate 3-year mean cover value
#Split the data by site using unique plot code
cover_list <- split( cover, cover$plot_code )   

#Take the average of the last 3 years cover, keeping the pre-treatment year
cover_means <- lapply( cover_list,              
                       function(plot){
                         plot_mean <- plot[ plot$year_trt != "0", ] %>%
                           dplyr::group_by( .dots = colnames ( plot ) [ c(1:2, 5:12, 14) ] ) %>%
                           summarise_all( mean );
                         plot_data <- bind_rows( plot_mean, 
                                                 plot[ plot$year_trt == "0" , ] );
                         plot_data
                       } )                      

#Bind rows together to make final dataset
cover_data <- data.frame( bind_rows( cover_means ) )  

#Rename year to distinguish 3-year mean from pretreatment year
cover_data <- cover_data %>%                          
  mutate( year = as.factor(
    case_when( year_trt == "0" ~ "initial",
               year_trt != "0" ~ "mean" ) ) )
              

#Load in GS data
Cval <- read.csv( "GS-by-site.csv" )

#Generate unique code for each species on each site
Cval$site_species_code <- paste( Cval$site_code, 
                                 Cval$taxon, 
                                 sep = "." )
cover_data$site_species_code <- paste( cover_data$site_code, 
                                       cover_data$taxon, 
                                       sep = "." )

#Now merge the GS data and cover data

#Generate shorter GS data, then merge with the cover data
Cval_final <- Cval[ , c(10, 6, 3:5, 8:9) ]

cover_GS <- merge( cover_data, 
                   Cval_final, 
                   by = "site_species_code", 
                   all.x = T )

#Clean up data
cover_GS <- cover_GS[ , c( 2:5, 13, 6:11, 17:19, 16, 20:21, 15, 12) ]

#Remove non-angiosperms and "unknowns"
cover_final <- cover_GS[ !is.na( cover_GS$species_GS_code ) , ]       #Data we will use
cover_NA <- cover_GS[ is.na( cover_GS$species_GS_code ) , ]           #To check - should be unknowns, ferns, bryophytes and fungi

write.csv( cover_final, "27-site-species-final-data.csv")

###====Selection of climate variables====
require(ggbiplot)
require(FactoMineR)
require(factoextra)
require(corrplot)

#This will extract the bioclim variables for the sites we are interested in and use a PCA to identify climate variables to include in our analysis
  
climate <- read.csv("27-site-climate-data.csv", row.names = 1)

cor_matrix <- cor( climate[ , c( 10:28 ) ] )

write.csv(cor_matrix, "Climate Correlation Matrix.csv")

#Generate a PCA for all climate variables
pca <- prcomp( climate[ , c( 10:28 ) ],
               center = TRUE, 
               scale. = TRUE )

#Examine loadings/component importance
pca
summary( pca )

#Plot PCA to observe variables
ggbiplot( pca, 
          alpha = 0)+
  geom_point( aes( colour = climate$Country ), 
              size = 3 )+
  labs( colour = "Country" )+
  xlim( -2.3, 2 )+
  theme_bw()+
  theme( plot.margin = margin( 0, 0, 0, 0 ) )+
  theme( axis.text = element_text( size = 16 ),
         axis.title = element_text( size = 16 ),
         legend.text = element_text( size = 16 ),
         legend.title = element_text( size = 16 ) )

#save to files
ggsave("All climate variables PCA.png",
       dpi = 600)

#Isolate the temperature variables and run PCA
pcatemp <- prcomp( climate[ , c( 10:20 ) ], 
                   center = TRUE,
                   scale. = TRUE ) 

pcatemp
summary(pcatemp)

#We can also plot that as above - uses modified version of ggbiplot, see attached file and run
ggbiplot2( pcatemp, alpha = 0,
           varname.size = 6.5,
           varname.adjust = 1.1)+
  geom_point( aes( colour = climate$Country ), size = 4.5)+
  labs( colour = "Country" )+
  xlim( -2.7, 2.5 )+
  theme_bw()+
  theme( plot.margin = margin( 20, 0, 0, 0 ) )+
  theme( axis.text = element_text( size = 25 ),
         axis.title = element_text( size = 25 ),
         legend.position = "none",
         panel.border = element_rect(linewidth = 1.34))

#Save temperature PCA
ggsave( "Temperature variables PCA.png",
       dpi = 600 )

#Look at the contributions of each variable to the first two PCs
temp.pca <- PCA( climate[ , c( 10 : 20 ) ],
                 ncp = 3,
                 graph = TRUE )

var <- get_pca_var( temp.pca ) 

corrplot( var$cos2,
          is.corr=FALSE )

#Isolate the precipitation variables and run PCA
pcaprecip <- prcomp( climate[ , c( 21 : 28 ) ], 
                        center = TRUE, 
                        scale. = TRUE )
pcaprecip 
summary(pcaprecip )

#Plot precipitation variables PCA - uses modified version of ggbiplot, see attached file and run

ggbiplot2( pcaprecip, alpha = 0, 
           varname.size = 6.5,
           varname.adjust = 1.1)+
  geom_point( aes( colour = climate$Country ), size = 4.5)+
  labs( colour = "Country" )+
  theme_bw()+
  xlim( -2, 2.7 )+
  ylim(-2, 2)+
  theme( plot.margin = margin( 20, 0, 0, 0 ) )+
  theme( axis.text = element_text( size = 25 ),
         axis.title = element_text( size = 25 ),
         panel.border = element_rect(linewidth = 1.34))

ggsave( "Precipitation variables PCA.png",
        dpi = 600 )

#Look at the contributions of each variable to the first two PCs
precip.pca <- PCA( climate[ , c( 21 : 28 ) ],
                 ncp = 3,
                 graph = TRUE )

var <- get_pca_var(precip.pca ) 

corrplot( var$cos2,
          is.corr = FALSE )

###====Quantifying variation in genome size and percent cover across functional groups====

#This plots and analyses the variation in genome size and percent cover between functional groups in the data for the 27 sites. 

#Variation in GS across between functional groups in the Nutrient Network
GS_dat <- read.csv("GS-by-site.csv")

#Just get one row for each species
GS <- unique( GS_dat[ , c( 6, 8, 4:5 ) ] )

#Remove NAs
GS_final <- GS[ !is.na( GS$GS_1C ) , ]

#Fit model for variation in GS between life forms
GS_mod <- lm( log(GS_1C) ~ lifeform, 
              data = GS_final )

#Observe model output
anova(GS_mod)
summary(GS_mod)

#We can also use Tukey HSD tests to observe significant differences between them
GS_mod2 <- aov( log( GS_final$GS_1C ) ~ GS_final$lifeform  )

anova( GS_mod2 )

#Fit tukey test
TUKEY_GS <- TukeyHSD( x = GS_mod2, 'GS_final$lifeform', 
                   conf.level = 0.95 )
TUKEY_GS

#Look at pairwise comparisons at the  95% confidence level 
plot( TUKEY_GS , 
      las=1 , 
      col="brown" )


#To compare average GS between C3 and C4 grass species
#Subset the data for just the grasses
Grass_GS <- GS_final[ GS_final$lifeform == "Grass" , ]

#t.test for difference in GS between C3 and C4 grasses
t.test(GS_1C ~ Photosynthetic_Pathway, data = Grass_GS )

## To group species by those that are and are not significantly different

require(multcompView)

#Create a function to label significantly different groups
generate_label_df <- function( TUKEY, variable ){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame( multcompLetters( Tukey.levels )[ 'Letters' ] )
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$treatment = rownames( Tukey.labels )
  Tukey.labels = Tukey.labels[ order( Tukey.labels$treatment ) , ]
  return( Tukey.labels )
}

#Extract the significantly different groups
LABELS_GS <- generate_label_df( TUKEY_GS , 
                             'GS_final$lifeform' )
LABELS_GS


#Variation in percent cover across between functional groups in the Nutrient Network

cover_final <- read.csv("27-site-species-final-data.csv", row.names = 1)

#Extract only the data needed for this analysis
cover_dat <- cover_final[ , c( 2, 19, 5, 13, 18 ) ] 

#Split the data by each plot within each site, and by initial year or the mean of the last 3 years
cover_dat$plotyear_code <- paste( cover_dat$plot_code, cover_dat$year)
plotyear_list <- split( cover_dat, cover_dat$plotyear_code )

#Calculate the total percent cover by angiosperms on each plot in each time bracket
cover_totals <- bind_rows( lapply( plotyear_list, function( plot ) {
  total <- sum( plot$max_cover )
  data <- data.frame( plotyear_code = unique( plot$plotyear_code ),
                      total_cover = total )
} ) )

#Seperately, calculate the total cover taken up by each plant functional group on each plot
cover_sums <- cover_dat %>%
  group_by( plot_code, plotyear_code, year, site_code, lifeform ) %>%
  summarise( cover = sum( max_cover ) )

#Merge the data, keeping only data for the pretreatment year
cover_sums <- merge( cover_sums[ cover_sums$year == "initial", ], cover_totals, 
                     by = "plotyear_code" )

#Calculate the proportion of total cover taken up by each functional group
cover_sums$cover_pct <- ( cover_sums$cover / cover_sums$total_cover ) * 100

#Fit a model for how the proportion of total cover on a plot varies between functional groups
cover_mod <- lmer(cover_pct ~ lifeform + (1|plot_code), data = cover_sums)
anova(cover_mod)
summary(cover_mod)

#Use TUKEY HSD to label significantly different groups (as above)
cover_mod2 <- aov(cover_sums$cover_pct ~ cover_sums$lifeform)

TUKEY_cover <- TukeyHSD( x= cover_mod2, 
                    'cover_sums$lifeform', 
                    conf.level=0.95 )

TUKEY_cover

plot( TUKEY_cover ,
      las=1 , 
      col="brown" )

#Extract the significantly different groups
LABELS_cover <- generate_label_df( TUKEY_cover , 
                             'cover_sums$lifeform' )
LABELS_cover

###====Analysis of cover-weighted genome size====

require(lme4)
require(lmerTest)
require(MuMIn)
require(nlme)

#This calculates the cover-weighted mean 1C value (cwGS) for each plot on each site, calculates the change in cwGS from control and pretreatment conditions using log response ratios, and fits these change in cwGS values against nutrients and the selected climate variables using linear mixed-effects ANOVA

cover_final <- read.csv( "27-site-species-final-data.csv", 
                         row.names = 1 )

#Split site by each plot in each year 
cover_final$plotyear_code <- paste( cover_final$plot_code, 
                                    cover_final$year, 
                                    sep = "_" )
cover_year_list <- split( cover_final, cover_final$plotyear_code )

#Calculate species richness and richness of species with GS values
cwGS1 <- lapply( cover_year_list,   
                   function( plot ){
                     GS_rich  = nrow( plot[ plot$max_cover != 0 & !is.na(plot$GS_1C) , ] );
                     plot$GS_rich = GS_rich;
                     sp_rich = nrow( plot[ plot$max_cover != 0 , ] );
                     plot$sp_rich = sp_rich;
                     plot
                   })
check <- bind_rows( cwGS1 )

#Calculate weighted GS 
cwGS2 <- bind_rows( lapply( cwGS1, 
                             function( plot ){
                               if( plot$GS_rich != 1 ){
                                 wtmean = nlme::gls( log10( GS_1C ) ~ 1,
                                                    weights = nlme::varFixed( ~1/max_cover ),
                                                    data = plot,
                                                    na.action =  na.omit );
                                 wtmeanGS = 10^( wtmean$coeff ); 
                               } else{
                                 wtmeanGS = mean( plot$GS_1C[ plot$max_cover != 0 ], na.rm = T );
                               }
                               row = data.frame( plotyear_code = unique( plot$plotyear_code ),
                                                plot_cwGS = wtmeanGS )
                               row
                             } ) )

cwGS3 <- bind_rows(cwGS1)

cwGS <- merge(cwGS3, cwGS2, by = "plotyear_code", all.x = T)

#Seperate out the plot level data
plot_data <- unique( cwGS[ , c( 1:10, 20:23 ) ] )

#Load in Table S1 to add this to plot data
site_background <- read.csv( "27-site-background-data.csv" )

plot_data <- merge( plot_data, site_background[ , c( 1, 6, 10, 13, 21, 24, 34 ) ] )

#Save the csv file for later
write.csv( plot_data, "27-site-plot-data.csv" )

plot_data <- read.csv( "27-site-plot-data.csv", row.names = 1 )

#Now to calculate the change in weighted genome size - from initial and from control 
#From initial
plot_list <- split( plot_data, plot_data$plot_code ) 

plot_initial <- bind_rows( lapply( plot_list, 
                            function( plot ){
                              plot$initial_cwGS <- plot$plot_cwGS[ plot$year == "initial" ];
                              plot
                            } ) )

plot_list2 <- split( plot_initial, plot_initial$site_code )

#From control 
plot_control <- bind_rows( lapply( plot_list2, 
                            function( plot ){
                              plot$control_cwGS <- mean( plot$plot_cwGS[ plot$trt == "Control" & plot$year == "mean" ] );
                              plot
                            } ) )

#Calculate the LRRs
plot_final <- plot_control[ plot_control$year == "mean" , ]
plot_final$LRR_tp <- log( plot_final$plot_cwGS / plot_final$initial_cwGS )
plot_final$LRR_tc <- log( plot_final$plot_cwGS / plot_final$control_cwGS )

#Put the NP variable in order
plot_final$NP <- factor(plot_final$NP, levels = c("Control", "N", "P", "NP"))

#Preliminary plots
plot( LRR_tc ~ NP, data = plot_final )
plot( LRR_tp ~ NP, data = plot_final )

#Run linear mixed effects models for both LRRs
tc_model <- lmer(LRR_tc ~ N * P + (1|site_code/block), 
                 data = plot_final)

#Get ANOVA outputs and summary
anova(tc_model)
summary(tc_model)
#Look at R-squared
r.squaredGLMM(tc_model)

tp_model <- lmer(LRR_tp ~ N * P + (1|site_code/block), 
                 data = plot_final)

#Get ANOVA outputs and summary
anova(tp_model)
summary(tp_model)
#Look at R-squared
r.squaredGLMM(tp_model)

#Now to incorporate climate variables into the model
#First scale variables as they vary in range/magnitude
plot_final$mat <- scale( plot_final$MAT )
plot_final$temp_var <- scale( plot_final$temp_seasonality )
plot_final$map <- scale( plot_final$MAP )
plot_final$map_var <- scale( plot_final$precip_seasonality )

#Run the model including 2 way interactions betweeen nutrients and each climate variable
LRR_clim_mod <- lmer( LRR_tc ~ N*P*temp_var + N*P*mat + N*P*map + N*P*map_var + (1|site_code/block),
                     data = plot_final )

#Get ANOVA outputs and summary
anova( LRR_clim_mod )
summary( LRR_clim_mod )

#Look at R-squared
r.squaredGLMM( LRR_clim_mod )


###====Analysis of species richness====

#This calculates the species richness for each plot on each site, and the richness of species with GS less than and greater than the mean GS value for that site. It then fits these richness values against nutrients and the selected climate variables using linear mixed-effects ANOVA

cover_final <- read.csv( "27-site-species-final-data.csv", 
                         row.names = 1 )

plot_data <- read.csv( "27-site-plot-data.csv", 
                       row.names = 1 )


#First look at richness independent of genome size

#Ensure NP factor is organised correctly
plot_data$NP <- factor(plot_data$NP,
                       levels = c( "Control", "N", "P", "NP" ) )


#Fit linear mixed effect model
richness_mod <- lmer( sp_rich ~ N * P + (1|site_code/block),
                      data = plot_data )

#View model output
anova( richness_mod )
summary( richness_mod )

#Now to establish richness of species with larger and smaller GS 

#First categorise GS of each species relative to mean GS for plot
#Calculate plot mean GS
site_list <- split( cover_final, cover_final$site_code )          #Split by site

site_GS <- bind_rows( lapply( site_list, function( site ){        #Calculate mean GS for each site
  GS_data = unique( site[ , c( 1, 15:16 ) ] );
  site_meanGS = mean( GS_data$GS_1C, 
                      na.rm = T );
  data = data.frame( site_code = unique( site$site_code ),
                     site_meanGS = site_meanGS );
  data
} ) )

cover_final <- merge( cover_final, site_GS, by = "site_code" )    #Merge with the other data

#Assign GS categories for each species
cover_final <- cover_final %>%
  mutate( GS_Cat = as.factor(
    case_when( GS_1C > site_meanGS ~ "Large",
               GS_1C < site_meanGS ~ "Small",
               is.na(GS_1C) ~ "Unknown" ) ) )

cover_final$GS_Cat <- factor( cover_final$GS_Cat, 
                              levels = c( "Unknown", "Small", "Large" ) )

#Now calculate the richness of species in each GS category
#Select only rows with non-zero cover values (i.e. actually present)
cover_richness <- cover_final[ cover_final$max_cover != 0, ]

#Create a code unique code for each category on each site in each year 
cover_richness$plotyear_code <- paste( cover_richness$plot_code, 
                                       cover_richness$year,
                                       sep = "_" )
cover_richness$GS_Cat_code  <- paste( cover_richness$plotyear_code,
                                      cover_richness$GS_Cat,
                                      sep = "_" )

GS_cat_list <- split(cover_richness, cover_richness$plotyear_code)

#Tally species in each GS category
tally <- bind_rows( lapply( GS_cat_list, function( plot ){
  tally <- plot %>% group_by( GS_Cat_code ) %>% tally()
  tally
}))

#Merge with nutrient and plot information
Richness <- merge( unique( cover_richness[ , c( 1, 3:9, 19:23 ) ] ), tally, 
                   by = "GS_Cat_code", 
                   all.x = T )
Richness <- rename( Richness, GS_richness = n )

#Merge with climate data
richness <- merge( Richness, plot_data[ , c( 2, 11:12, 14:19 ) ],
                         by = "plotyear_code",
                         all.x = T )

#We'll fit models for after treatment only, not pretreatment
richness_final <- richness[ richness$year == "mean" , ]

#Ensure factors are ordered as expected
richness_final$NP <- factor( richness_final$NP,
                             levels = c( "Control", "N", "P", "NP" ) )

#We don't need unknowns 
richness_final$GS_Cat[richness_final$GS_Cat == "Unknown"] <- NA
richness$GS_Cat <- factor( richness$GS_Cat,
                           levels = c( "Small", "Large" ) )


#Fit model of richness against nutrient treatment and GS
GSrich_model <- lmer(GS_richness ~ N * P * GS_Cat + (1|site_code/block),
                     data = richness_final)

#View model output
anova(GSrich_model)
summary(GSrich_model)


#Now to incorporate climate variables into the model
#First scale variables as they vary in range/magnitude
richness_final$mat <- scale( richness_final$MAT )
richness_final$temp_var <- scale( richness_final$temp_seasonality )
richness_final$map <- scale( richness_final$MAP )
richness_final$map_var <- scale( richness_final$precip_seasonality )

#Run the model including 2 way interactions betweeen nutrients and each climate variable
richclim_mod <- lmer(GS_richness ~ N*P*GS_Cat*mat + 
                       N*P*GS_Cat*map + 
                       N*P*GS_Cat*temp_var + 
                       N*P*GS_Cat*map_var + 
                       (1|site_code/block),
                     data = richness_final)

#Get ANOVA outputs and summary
anova(richclim_mod)
summary(richclim_mod)


###====Compiling phylogeny====

#This loads in a phylogeny of species found across the Nutrient Network, pruned from a larger megaphylogeny (Qian & Jim, 2016). It also further prunes the phylogeny to make smaller phylogenies for each of the 6 functional groups used in this study.

# Qian, H. & Jin, Y. (2016). An updated megaphylogeny of plants, 
#   a tool for generating plant phylogenies and an analysis of 
#   phylogenetic community structure. J Plant Ecol, 9, 233–239.

require(ape)
require(phytools)

#Load in tree
NutNet_tree <- read.tree("NutNet Tree.tre")

#Alter names to match up with the syntax of the .tre file
cover_final <- read.csv( "27-site-species-final-data.csv", row.names = 1 )

cover_final$phylo <- sub( "(\\w+\\s+\\w+).*", "\\1",
                          cover_final$taxon )           #Remove subspecies information

cover_final$phylo <- gsub( " ", "_", 
                           cover_final$phylo )          #Replace spaces with _

#Get list of species across NutNet (just the species with GS values)
sp_list <- unique( cover_final$phylo[ !is.na(cover_final$GS_1C) ] , ) 

#Remove some of the oddly named bits
NutNet_species <- data.frame( NutNet_tree$tip.label )

# get the tip labels corresponding to species
tip_labels = lapply( 1:length( sp_list ), 
                     function( x ) grep( paste( sp_list[[x]], collapse='.+' ),
                                         NutNet_tree$tip.label, ignore.case=T ) )

# show the species that haven't been found
not_found = which( sapply( tip_labels, length )==0 )
sp_list[ not_found ]

#Make sure there are no NAs in the tip labels 
tip_labels <- sapply( tip_labels, function( x ) x[ 1 ] )
not_NA <- !is.na( tip_labels )
tip_labels <- tip_labels[ not_NA ]
NutNet_tree$tip.label[ tip_labels ] <- sp_list[ not_NA ]

#Construct the tree
phylogeny <- keep.tip( NutNet_tree, tip=tip_labels )

#Save files
write.tree( phylogeny, "27_site_phylogeny.tre" )
write.csv( cover_final, "27-site-phylo-data.csv" )

#Now prune phylogeny for each separate brms analysis

#1 Prune phylogeny to only include species we measured GS values for

measured <- cover_final[ cover_final$Source == "Measured" , ]

#Compile species list
measured_list <- unique(measured$phylo)

# get the tip labels corresponding to species
tips_measured = lapply(1:length(measured_list), 
                       function(x) grep(paste(measured_list[[x]], collapse='.+'),
                                        phylogeny$tip.label, ignore.case=T))

# show the species that haven't been found
not_found = which(sapply(tips_measured, length)==0)
measured_list[not_found]

#Make sure there are no NAs in the tip labels 
tips_measured <- sapply( tips_measured, function( x ) x[ 1 ] )
not_NA <- !is.na( tips_measured )
tips_measured <- tips_measured[ not_NA ]
NutNet_tree$tip.label[ tips_measured ] <- sp_list[ not_NA ]

#Construct and write the tree
measured_tree <- keep.tip(phylogeny, tip=tips_measured)
write.tree(measured_tree, "measured-GS-tree.tre")


#Now prune for grasses

cover_noNA <- cover_final[!is.na(cover_final$GS_1C),]

grass <- cover_noNA[ cover_noNA$lifeform == "Grass" , ]

#Compile species list
grass_list <- unique(grass$phylo)

# get the tip labels corresponding to species
tips_grass = lapply(1:length(grass_list), 
                    function(x) grep(paste(grass_list[[x]], collapse='.+'),
                                     phylogeny$tip.label, ignore.case=T))

# show the species that haven't been found
not_found = which(sapply(tips_grass, length)==0)
grass_list[not_found]

#Make sure there are no NAs in the tip labels 
tips_grass <- sapply( tips_grass, function( x ) x[ 1 ] )
not_NA <- !is.na( tips_grass )
tips_grass <- tips_grass[ not_NA ]
phylogeny$tip.label[ tips_grass ] <- sp_list[ not_NA ]

#Construct and write the tree
grass_tree <- keep.tip(phylogeny, tip=tips_grass)
write.tree(grass_tree, "grass-GS-tree.tre")

###====Analysis of species cover responses to nutrients====

#This calculates the change in percent cover of each species and fits this against nutrient treatment and GS in a phylogenetic linear mixed effect model, across all species, species for which we directly measured genome size, and for each functional group individually

require(brms)
require(ape)
require(ggplot2)

cover_phylo <- read.csv( "27-site-phylo-data.csv", 
                         row.names = 1 )

#Remove odd duplicates and species with no GS data
cover_phylo <- unique( cover_phylo[ !is.na( cover_phylo$GS_1C ), ] )

#First calculate change in cover from initial value
cover_phylo$plotspecies_code <- paste( cover_phylo$plot_code, 
                                       cover_phylo$taxon, 
                                       sep = "_" )

phylo_list <- split( cover_phylo, cover_phylo$plotspecies_code )

#New column with initial cover 
phylo_data <- bind_rows( lapply( phylo_list, function( species ){
  species$initial_cover = species$max_cover[species$year == "initial" ];
  species
} ) )

#Calculate change in percent cover
phylo_data$delta_cover <- phylo_data$max_cover - phylo_data$initial_cover

#Get together final dataset by removing any species that were not present before or after treatment
phylo_final <- phylo_data[ phylo_data$year == "mean" , ]
phylo_final$presence <- 1
phylo_final$presence[ phylo_final$max_cover == "0" & phylo_final$initial_cover == "0" ] <- 0
phylo_final <- phylo_final [phylo_final$presence == "1" , ]

#Now to build models

#1 Overall model
phylogeny <- read.tree("27_site_phylogeny.tre")

#log transform GS 
phylo_final$log.GS <- log(phylo_final$GS_1C)

#Get the vcv matrix up 
A <- ape:::vcv.phylo(phylogeny)

#Run the model with log.GS
full_brms <- brm(
  delta_cover ~ log.GS * N * P +                    #Fixed effects
    ( 1|gr( phylo, cov = A ) ) +                 #Phylogenetic covariance matrix as random effect
    ( 1|site_code/block ),                       #Block nested within site as random effect
  data = phylo_final,                            #Data source
  family = gaussian(),                           #Distribution of data
  data2 = list( A = A ),                         #Source of variance-covariance matrix for phylogeny
  prior = c( prior( normal( 0, 1 ), "b" ) ),
  sample_prior = TRUE, chains = 3, 
  iter = 20000, warmup = 8000,                   #Define any priors for the model
  cores = getOption( "mc.cores", 3 ),            #Run on 3 cores to speed up models
  control = list( adapt_delta = 0.9,
                  max_treedepth = 15 ) )         #Set based on preliminary runs to ensure convergence

#View model output
summary( full_brms ) 
#View posterior distributions and check for convergence (caterpillar plots)
plot( full_brms, N = 2, ask = TRUE )       
#View model outputs graphically
plot( conditional_effects( full_brms ), points = TRUE )


#2 Measured GS model
measured <- phylo_final[phylo_final$Source == "Measured",]
measured_tree <- read.tree("measured_GS_tree.tre")

#Get the vcv matrix up 
A <- ape:::vcv.phylo(measured_tree)

#Run the model with log.GS
measured_brms <- brm(
  delta_cover ~ log.GS * N * P + 
    ( 1|gr( phylo, cov = A ) ) + 
    ( 1|site_code/block ), 
  data = measured, 
  family = gaussian(), 
  data2 = list( A = A ),
  prior = c( prior( normal( 0, 1 ), "b" ) ),
  sample_prior = TRUE, chains = 3, 
  iter = 20000, warmup = 8000, 
  cores = getOption( "mc.cores", 3 ),
  control = list( adapt_delta = 0.9,
                 max_treedepth = 15 ) )

summary( measured_brms )
plot( measured_brms, N = 2, ask = TRUE )
plot( conditional_effects( measured_brms ), points = TRUE )


#3 Overall model including functional group as a interacting factor
phylogeny <- read.tree("27_site_phylogeny.tre")
                     
#Get the vcv matrix up 
A <- ape:::vcv.phylo(phylogeny)

#Run the model with log.GS
lifeform_brms <- brm(
  delta_cover ~ log.GS * N * P * lifeform +      #Fixed effects
    ( 1|gr( phylo, cov = A ) ) +                 #Phylogenetic covariance matrix as random effect
    ( 1|site_code/block ),                       #Block nested within site as random effect
  data = phylo_final,                            #Data source
  family = gaussian(),                           #Distribution of data
  data2 = list( A = A ),                         #Source of variance-covariance matrix for phylogeny
  prior = c( prior( normal( 0, 1 ), "b" ) ),
  sample_prior = TRUE, chains = 3, 
  iter = 20000, warmup = 8000,                   #Define any priors for the model
  cores = getOption( "mc.cores", 3 ),            #Run on 3 cores to speed up models
  control = list( adapt_delta = 0.9,
                  max_treedepth = 15 ) )         #Set based on preliminary runs to ensure convergence

#View model output
summary( lifeform_brms ) 
#View posterior distributions and check for convergence (caterpillar plots)
plot( lifeform_brms, N = 2, ask = TRUE )       
#View model outputs graphically
plot( conditional_effects( full_brms ), points = TRUE )

                     
#4 Grass photosynthesis model
photosynthesis <- phylo_final[phylo_final$Source == "Grass",]
grass_tree <- read.tree("grass_GS_tree.tre")

#Get the vcv matrix up 
A <- ape:::vcv.phylo(grass_tree)

#Run the model with log.GS
photosynth_brms <- brm(
  delta_cover ~ log.GS * NP * Photosynthetic_Pathway + 
    ( 1|gr( phylo, cov = A ) ) + 
    ( 1|site_code/block ), 
  data = woody, 
  family = gaussian(), 
  data2 = list( A = A ),
  prior = c( prior( normal( 0, 1 ), "b" ) ),
  sample_prior = TRUE, chains = 3, 
  iter = 20000, warmup = 8000, 
  cores = getOption( "mc.cores", 3 ),
  control = list(adapt_delta = 0.9,
                 max_treedepth = 15))

summary( photosynth_brms )
plot( photosynth_brms, N = 2, ask = TRUE )
plot( conditional_effects( photosynth_brms ), points = TRUE )

