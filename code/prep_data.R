# Description: Clean and tabulate input survey microdata for analysis
# Inputs: DHS and PMA 2020 microdata (prepped for Global Burden of Disease)
# Outputs: Dataset with admin 1-level tabulations for analysis
rm(list = ls())
if(!require(pacman)) install.packages('pacman')
pacman::p_load(data.table, magrittr, haven, survey, sp, maptools, rgdal, ggplot2, stringr)

# Read in dataset(s) - start with most recent dhs

data_df <- read_dta('../data/MACRO_DHS_WN_ETH_2016_2016_218568.dta') %>% as.data.table

# Using gbd criteria to determine need of contraception: using contraception or fecund/sexually active/no wish to become pregnant in next two years
data_df[, not_fecund := (pregnant == 1|
                           (never_used_contra == 1 & years_since_first_cohabit > 5 & is.na(months_since_last_birth))|
                           no_menses == 1|
                           desire_children_infertile == 1)][is.na(not_fecund), not_fecund := F] #

data_df[, contra_need := (current_contra != 'not using') | (!not_fecund & sex_in_last_month == 1 & desire_later == 1)][is.na(contra_need), contra_need := F]

appendix_dhs_1 <- copy(data_df[, .N, by = .(5*age_year%/%5, admin_1)])[, source := 'dhs_total']
appendix_dhs_2 <- copy(data_df[, .N, by = .(5*age_year%/%5, admin_1, contra_need)])[, source := 'dhs_need_status']

data_df <- data_df[contra_need == T]
data_df[, modern := !current_contra %like% 'periodic abstinence|not using|standard days|lactational|withdrawal' %>% as.integer]

dhs_design <- svydesign(id = ~psu + hh_id, strata = ~strata, weights = ~pweight, data = data_df)
dhs_collapse <- svyby(~modern, ~admin_1, dhs_design, svymean)
dhs_collapse <- data.table(source = 'dhs', year = 2016, admin_1 = dhs_collapse$admin_1, prev_met_mod_need = dhs_collapse$modernTRUE, prev_met_mod_need_se = dhs_collapse$se.modernTRUE)

dhs_sample_sizes <- data_df[, .N, by = admin_1]
dhs_collapse <- merge(dhs_collapse, dhs_sample_sizes, by = 'admin_1')
dhs_collapse[admin_1 %like% 'addis adaba', admin_1 := 'addis ababa']

data_df <- read_dta('../data/JHSPH_PERFORMANCE_MONITORING_ACCOUNTABILITY_SURVEY_PMA2020_WN_ETH_2016_2016_285891.dta') %>% as.data.table

# Filter down to correct ages
data_df <- data_df[age_year >= 15 & age_year <= 49]

# Clean dates
data_df[, interview_date := str_extract(interview_date, "^.*2016") %>% as.Date(., '%B %d, %Y')]
data_df[, recent_cohabit_start_date := str_extract(recent_cohabit_start_date, "^.*2016") %>% as.Date(., '%B %d, %Y')]

# Clean last sexual activity variable - express all in terms of months
data_df[, sexually_active := (last_sex_unit == 1 & last_sex/30 < 1)|
          (last_sex_unit == 2 & last_sex/4 < 1)|
          (last_sex_unit == 3 & last_sex < 1)|
          (last_sex_unit == 4 & last_sex*12 < 1)][is.na(sexually_active), sexually_active := F]

# Clean desire timing variable
data_df[, desire_later := (desire_unit == 1 & desire_timing > 24)|
          (desire_unit == 2 & desire_timing > 2)|
          desire_timing %in% c(-88, 5)|
          desire_gate != 1][is.na(desire_later), desire_later := F]

# Using gbd criteria to determine need of contraception: using contraception or fecund/sexually active/no wish to become pregnant in next two years
data_df[, not_fecund := (pregnant == 1|
                           (never_used_contra == 1 & is.na(last_birth_date))|
                           no_menses == 1|
                           desire_children_infertile == 1)][is.na(not_fecund), not_fecund := F] #

data_df[, contra_need := (current_contra != '') | (!not_fecund & sexually_active & desire_later)][is.na(contra_need), contra_need := F]

appendix_pma_1 <- copy(data_df[, .N, by = .(5*age_year%/%5, admin_1)])[, source := 'pma_total']
appendix_pma_2 <- copy(data_df[, .N, by = .(5*age_year%/%5, admin_1, contra_need)])[, source := 'pma_need_status']


data_df <- data_df[contra_need == T]
data_df[, modern := !current_contra %like% '^other_traditional|^rhythm|^withdrawal|^LAM|beads' %>% as.integer][current_contra == "", modern := F]
data_df <- data_df[!is.na(hhweight) & !is.na(pweight)]

pma_design <- svydesign(id = ~psu + ~hh_id, strata = ~strata, weights = ~hhweight + ~pweight, data = data_df)
pma_collapse <- svyby(~modern, ~admin_1, pma_design, svymean)
pma_collapse <- data.table(source = 'pma2020', year = 2016, admin_1 = pma_collapse$admin_1, prev_met_mod_need = pma_collapse$modernTRUE, prev_met_mod_need_se = pma_collapse$se.modernTRUE)

pma_sample_size <- data_df[, .N, by = .(admin_1)]
pma_collapse <- merge(pma_collapse, pma_sample_size, by = 'admin_1')

pma_collapse[, admin_1 := gsub('[^[:alpha:] ]', '', admin_1) %>% trimws]
pma_collapse[admin_1 %like% 'benishangul', admin_1 := 'benishangul']
pma_collapse[admin_1 %like% 'somali', admin_1 := 'somali']
pma_collapse[admin_1 %like% 'oromiya', admin_1 := 'oromia']
pma_collapse[admin_1 %like% 'gambella', admin_1 := 'gambela']

all_prepped <- rbind(pma_collapse, dhs_collapse, use.names = T, fill = T)
all_prepped <- all_prepped[order(source, admin_1)]

## Read in shape file from lbd
# shape_file <- readRDS('../data/lbd_standard_admin_1.rds') #too large to push to github but prepped file is still in repo
shape_file <- readRDS('../data/eth_shp_prepped.rds')
ethiopia_codes <- shape_file@data %>% as.data.table %>% .[ADM0_NAME == 'Ethiopia']
shape_file <- subset(shape_file, ADM1_CODE %in% ethiopia_codes$ADM1_CODE)

all_prepped[, ADM1_CODE := rep(ethiopia_codes$ADM1_CODE[c(1:8, 10, 9, 11)], 2)]

write.csv(all_prepped, '../data/prepped.csv', row.names = F)
saveRDS(shape_file, '../data/eth_shp_prepped.rds')

all_appendix <- rbind(appendix_dhs_1, appendix_dhs_2, appendix_pma_1, appendix_pma_2, use.names = T, fill = T)
all_appendix[, admin_1 := gsub('[^[:alpha:] ]', '', admin_1) %>% trimws]
all_appendix[admin_1 %like% 'benishangul', admin_1 := 'benishangul']
all_appendix[admin_1 %like% 'somali', admin_1 := 'somali']
all_appendix[admin_1 %like% 'oromiya', admin_1 := 'oromia']
all_appendix[admin_1 %like% 'gambella', admin_1 := 'gambela']
all_appendix[admin_1 %like% 'addis', admin_1 := 'addis ababa']

all_appendix <- all_appendix[order(admin_1, age_year, -source)]
setnames(all_appendix, c('age_year', 'source', 'contra_need'), c('Age', 'Source subset', 'Contraception Need'))
all_appendix[, `Source subset`:= paste0(`Source subset`, ': ', `Contraception Need`)][, `Contraception Need`:= NULL]
all_appendix <- dcast(all_appendix, ...~admin_1, value.var = 'N')

write.csv(all_appendix, '../data/appendix_survey_pops.csv', row.names = F)
