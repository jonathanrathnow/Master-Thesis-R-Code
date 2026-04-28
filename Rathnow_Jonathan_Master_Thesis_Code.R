# ============================================================
# MASTER THESIS ANALYSIS SCRIPT
# ============================================================
# A Zeitenwende in Public Opinion?
# Party Politics and the New Structure of Support 
# for Defence Strengthening in Germany since 2022
#
# Hertie School
# Master of Public Policy
# Jonathan Rathnow
# Year of Graduation: 2026
# Thesis Advisor: Prof. Dr. Tobias Bunde
#
# Data: German Longitudinal Election Study (GLES) Panel
#
# This script prepares the GLES data and reproduces the
# descriptive statistics, regression models, figures, and
# robustness checks used in the thesis.
# ============================================================



# ============================================================
# 0. SETUP
# ============================================================

library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(estimatr)
library(broom)
library(stargazer)
library(sandwich)
library(car)
library(MASS)

# Set path to the folder containing the GLES .sav files.
data_path <- "~/Library/CloudStorage/OneDrive-HertieSchool/Anlagen/Master Thesis Data/GLES Data/Datasets/GLES_waves"


# ============================================================
# 1. LOAD RAW DATA
# ============================================================
# The analysis uses eight GLES panel waves:
# W8, W10, W11, and W17 before the 2022 invasion;
# W22, W24, W26, and W32 after the 2022 invasion.
# Three additional profile files contain the demographic variables
w1to9 <- read_sav(file.path(data_path, "ZA6838_w1to9_sA_v6-0-0.sav"))
w10   <- read_sav(file.path(data_path, "ZA6838_w10_sA_v6-0-0.sav"))
w11   <- read_sav(file.path(data_path, "ZA6838_w11_sA_v6-0-0.sav"))
w17   <- read_sav(file.path(data_path, "ZA6838_w17_sA_v6-0-0.sav"))

w22   <- read_sav(file.path(data_path, "ZA7728_w22_v1-0-0.sav"))
w24   <- read_sav(file.path(data_path, "ZA7730_w24_v1-0-0.sav"))
w26   <- read_sav(file.path(data_path, "ZA7732_w26_sA_v2-0-0.sav"))
w32   <- read_sav(file.path(data_path, "ZA10121_w32_sA_v1-0-0.sav"))

profile_a1 <- read_sav(file.path(data_path, "ZA6804_de_v7-0-0.sav"))
profile_a2 <- read_sav(file.path(data_path, "ZA6838_wa2_sA_v6-0-0.sav"))
profile_a5 <- read_sav(file.path(data_path, "ZA7961_wa5_sA_v1-0-0.sav"))

cat("Raw wave and profile files loaded.\n")
cat("Profile rows: A1 =", nrow(profile_a1),
    "| A2 =", nrow(profile_a2),
    "| A5 =", nrow(profile_a5), "\n")


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================
# The functions below are used repeatedly when building the
# cleaned analysis dataset. They keep the recoding rules in one
# place, which makes the rest of the script easier to read.


# Keep valid responses to the defence-spending item.
# The item is measured on a 1-5 scale. Other values are treated
# as missing.

clean_defence <- function(x) {
  x_num <- as.numeric(x)
  ifelse(x_num %in% 1:5, x_num, NA_real_)
}


# Recode party preference into the party categories used in the
# thesis. CDU and CSU are combined into one CDU/CSU category.
# Other parties and non-valid responses are set to missing.

clean_party <- function(x) {
  x_num <- as.numeric(x)
  
  case_when(
    x_num %in% c(1, 2, 3) ~ 1,    # CDU/CSU
    x_num == 4            ~ 4,    # SPD
    x_num == 5            ~ 5,    # FDP
    x_num == 6            ~ 6,    # Greens
    x_num == 7            ~ 7,    # Left
    x_num == 322          ~ 322,  # AfD
    x_num == 392          ~ 392,  # BSW
    TRUE                  ~ NA_real_
  )
}


# Attach readable party names to the recoded party variable.

party_name <- function(x) {
  case_when(
    x == 1   ~ "CDU/CSU",
    x == 4   ~ "SPD",
    x == 5   ~ "FDP",
    x == 6   ~ "Greens",
    x == 7   ~ "Left",
    x == 322 ~ "AfD",
    x == 392 ~ "BSW",
    TRUE     ~ NA_character_
  )
}


# Recode political interest so that higher values mean greater
# political interest. In the original GLES coding, lower values
# indicate stronger interest.

clean_pol_interest <- function(x) {
  x_num <- as.numeric(x)
  ifelse(x_num %in% 1:5, 6 - x_num, NA_real_)
}


# Build a cleaned wave-level dataset with the same variable names
# across waves. This allows all selected waves to be combined later.
#
# Arguments:
# df          = raw GLES wave dataset
# id_var      = respondent ID variable
# wave_label  = wave name used in the thesis
# year        = survey year
# dv_var      = defence-spending variable in that wave
# iv_var      = party-preference variable in that wave
# weight_var  = optional survey weight, where available

make_wave <- function(df, id_var = "lfdn", wave_label, year,
                      dv_var, iv_var, weight_var = NULL) {
  
  df %>%
    transmute(
      lfdn = .data[[id_var]],
      wave = wave_label,
      year = year,
      defence = clean_defence(.data[[dv_var]]),
      party_code = clean_party(.data[[iv_var]]),
      party = party_name(party_code),
      weight = if (!is.null(weight_var) && weight_var %in% names(df)) {
        as.numeric(df[[weight_var]])
      } else {
        NA_real_
      }
    )
}


# Extract robust standard errors from estimatr models for use in
# custom tables.

robust_se <- function(model) {
  sqrt(diag(model$vcov))
}


# Coalesce variables only if they exist in the respective dataset.
# This is useful for profile variables because variable names differ
# across GLES profile files.

coalesce_existing <- function(df, vars) {
  present <- vars[vars %in% names(df)]
  
  if (length(present) == 0) {
    return(rep(NA, nrow(df)))
  }
  
  out <- df[[present[1]]]
  
  if (length(present) > 1) {
    for (v in present[-1]) {
      out <- dplyr::coalesce(out, df[[v]])
    }
  }
  
  out
}

cat("Setup, data loading, and helper functions completed.\n")

# ============================================================
# 3. BUILD ANALYSIS DATASET
# ============================================================
# The eight selected GLES waves are harmonised into the same
# basic structure. W8 uses reported vote in the 2017 election;
# later waves use current vote intention.


# ------------------------------------------------------------
# 3.1 Clean wave-level datasets
# ------------------------------------------------------------

d8 <- make_wave(
  w1to9,
  wave_label = "W8",
  year = 2017,
  dv_var = "kp8_2880y",
  iv_var = "kp8_200ba"
)

d10 <- make_wave(
  w10,
  wave_label = "W10",
  year = 2018,
  dv_var = "kp10_2880y",
  iv_var = "kp10_190ba",
  weight_var = "wei5_mz"
)

d11 <- make_wave(
  w11,
  wave_label = "W11",
  year = 2019,
  dv_var = "kp11_2880y",
  iv_var = "kp11_190ba",
  weight_var = "wei5_mz"
)

d17 <- make_wave(
  w17,
  wave_label = "W17",
  year = 2021,
  dv_var = "kp17_2880y",
  iv_var = "kp17_190ba",
  weight_var = "wei17_mz"
)

d22 <- make_wave(
  w22,
  wave_label = "W22",
  year = 2022,
  dv_var = "kp22_2880y",
  iv_var = "kp22_190ba"
)

d24 <- make_wave(
  w24,
  wave_label = "W24",
  year = 2023,
  dv_var = "kp24_2880y",
  iv_var = "kp24_190ba"
)

d26 <- make_wave(
  w26,
  wave_label = "W26",
  year = 2024,
  dv_var = "kp26_2880y",
  iv_var = "kp26_190ba",
  weight_var = "kp26_wei_mz"
)

d32 <- make_wave(
  w32,
  wave_label = "W32",
  year = 2025,
  dv_var = "kp32_2880y",
  iv_var = "kp32_190ba",
  weight_var = "kp32_wei_mz"
)

cat("Cleaned wave datasets created.\n")
cat("Rows per cleaned wave:\n")
cat("d8 :", nrow(d8), "\n")
cat("d10:", nrow(d10), "\n")
cat("d11:", nrow(d11), "\n")
cat("d17:", nrow(d17), "\n")
cat("d22:", nrow(d22), "\n")
cat("d24:", nrow(d24), "\n")
cat("d26:", nrow(d26), "\n")
cat("d32:", nrow(d32), "\n")


# ------------------------------------------------------------
# 3.2 Identify respondents in the analysed waves
# ------------------------------------------------------------

analysis_ids <- bind_rows(
  d8  %>% dplyr::select(lfdn),
  d10 %>% dplyr::select(lfdn),
  d11 %>% dplyr::select(lfdn),
  d17 %>% dplyr::select(lfdn),
  d22 %>% dplyr::select(lfdn),
  d24 %>% dplyr::select(lfdn),
  d26 %>% dplyr::select(lfdn),
  d32 %>% dplyr::select(lfdn)
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(lfdn = as.numeric(lfdn))

cat("Unique respondent IDs across analysed waves:", nrow(analysis_ids), "\n")


# ------------------------------------------------------------
# 3.3 Extract profile information
# ------------------------------------------------------------
# Demographic controls are stored in profile files. Because the
# variable names differ across profile files, coalesce_existing()
# keeps the first available version of each variable.

profile_a1_core <- profile_a1 %>%
  dplyr::semi_join(analysis_ids, by = "lfdn") %>%
  dplyr::transmute(
    lfdn = as.numeric(lfdn),
    gender_raw     = coalesce_existing(., c("kpx_2280")),
    birth_year_raw = coalesce_existing(., c("kpx_2290")),
    education_raw  = coalesce_existing(., c("kp1_2320", "kpa1_2320")),
    bundesland_raw = coalesce_existing(., c("kp1_2601"))
  )

profile_a2_core <- profile_a2 %>%
  dplyr::semi_join(analysis_ids, by = "lfdn") %>%
  dplyr::transmute(
    lfdn = as.numeric(lfdn),
    gender_raw     = coalesce_existing(., c("kpx_2280")),
    birth_year_raw = coalesce_existing(., c("kpx_2290s")),
    education_raw  = coalesce_existing(., c("kpa2_2320")),
    bundesland_raw = coalesce_existing(., c("kpa2_2601"))
  )

profile_a5_core <- profile_a5 %>%
  dplyr::semi_join(analysis_ids, by = "lfdn") %>%
  dplyr::transmute(
    lfdn = as.numeric(lfdn),
    gender_raw     = coalesce_existing(., c("kpa5_2280")),
    birth_year_raw = coalesce_existing(., c("kpa5_2290s")),
    education_raw  = coalesce_existing(., c("kpa5_2320")),
    bundesland_raw = coalesce_existing(., c("kpa5_2601"))
  )

cat("Core profile tables created.\n")
cat("profile_a1_core rows:", nrow(profile_a1_core), "\n")
cat("profile_a2_core rows:", nrow(profile_a2_core), "\n")
cat("profile_a5_core rows:", nrow(profile_a5_core), "\n")


# ------------------------------------------------------------
# 3.4 Merge profile information
# ------------------------------------------------------------
# The three profile files are combined into one respondent-level
# table before the demographic controls are recoded.

profile_demographics <- analysis_ids %>%
  dplyr::left_join(
    profile_a1_core %>%
      dplyr::rename(
        gender_a1     = gender_raw,
        birth_year_a1 = birth_year_raw,
        education_a1  = education_raw,
        bundesland_a1 = bundesland_raw
      ),
    by = "lfdn"
  ) %>%
  dplyr::left_join(
    profile_a2_core %>%
      dplyr::rename(
        gender_a2     = gender_raw,
        birth_year_a2 = birth_year_raw,
        education_a2  = education_raw,
        bundesland_a2 = bundesland_raw
      ),
    by = "lfdn"
  ) %>%
  dplyr::left_join(
    profile_a5_core %>%
      dplyr::rename(
        gender_a5     = gender_raw,
        birth_year_a5 = birth_year_raw,
        education_a5  = education_raw,
        bundesland_a5 = bundesland_raw
      ),
    by = "lfdn"
  ) %>%
  dplyr::transmute(
    lfdn,
    gender_raw = dplyr::coalesce(
      as.numeric(gender_a1),
      as.numeric(gender_a2),
      as.numeric(gender_a5)
    ),
    birth_year_raw = dplyr::coalesce(
      as.numeric(birth_year_a1),
      as.numeric(birth_year_a2),
      as.numeric(birth_year_a5)
    ),
    education_raw = dplyr::coalesce(
      as.numeric(education_a1),
      as.numeric(education_a2),
      as.numeric(education_a5)
    ),
    bundesland_raw = dplyr::coalesce(
      as.numeric(bundesland_a1),
      as.numeric(bundesland_a2),
      as.numeric(bundesland_a5)
    )
  )

cat("Merged respondent-level profile table created.\n")


# ------------------------------------------------------------
# 3.5 Check demographic coverage
# ------------------------------------------------------------

profile_coverage <- profile_demographics %>%
  dplyr::summarise(
    n_ids = dplyr::n(),
    gender_nonmissing = sum(!is.na(gender_raw)),
    birth_year_nonmissing = sum(!is.na(birth_year_raw)),
    education_nonmissing = sum(!is.na(education_raw)),
    bundesland_nonmissing = sum(!is.na(bundesland_raw))
  )

print(profile_coverage)


# ------------------------------------------------------------
# 3.6 Clean demographic controls
# ------------------------------------------------------------
# These controls are later merged into gles_analysis by respondent ID.

profile_demographics_clean <- profile_demographics %>%
  dplyr::transmute(
    lfdn,
    gender = dplyr::case_when(
      gender_raw == 1 ~ "Male",
      gender_raw == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    birth_year = dplyr::case_when(
      birth_year_raw >= 1900 & birth_year_raw <= 2010 ~ birth_year_raw,
      TRUE ~ NA_real_
    ),
    education = dplyr::case_when(
      education_raw %in% c(1, 2) ~ "Low",
      education_raw %in% c(3, 4) ~ "Medium",
      education_raw == 5 ~ "High",
      TRUE ~ NA_character_
    ),
    east_west = dplyr::case_when(
      bundesland_raw %in% 1:10 ~ "West",
      bundesland_raw %in% 11:16 ~ "East",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::mutate(
    gender = factor(gender, levels = c("Male", "Female")),
    education = factor(education, levels = c("Low", "Medium", "High")),
    east_west = factor(east_west, levels = c("West", "East"))
  )

cat("Demographic controls cleaned.\n")

cat("\nGender:\n")
print(table(profile_demographics_clean$gender, useNA = "ifany"))

cat("\nEducation:\n")
print(table(profile_demographics_clean$education, useNA = "ifany"))

cat("\nEast/West:\n")
print(table(profile_demographics_clean$east_west, useNA = "ifany"))

cat("\nBirth year summary:\n")
print(summary(profile_demographics_clean$birth_year))


# ------------------------------------------------------------
# 3.7 Build gles_analysis
# ------------------------------------------------------------
# This is the main long-format dataset before controls are added.

gles_long <- bind_rows(d8, d10, d11, d17, d22, d24, d26, d32)

gles_analysis <- gles_long %>%
  filter(!is.na(defence), !is.na(party)) %>%
  mutate(
    party = factor(
      party,
      levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD", "BSW")
    )
  )

cat("gles_analysis created.\n")

cat("Usable analysis cases by wave:\n")
gles_analysis %>%
  group_by(wave) %>%
  summarise(n = n(), .groups = "drop") %>%
  print()

saveRDS(gles_analysis, file.path(data_path, "gles_analysis_ready.rds"))
write.csv(gles_analysis, file.path(data_path, "gles_analysis_ready.csv"), row.names = FALSE)

cat("gles_analysis saved.\n")


# ------------------------------------------------------------
# 3.8 Add political interest and post-2022 indicator
# ------------------------------------------------------------
# Political interest is wave-specific and reverse-coded earlier
# so that higher values indicate greater political interest.

pol_interest_all <- bind_rows(
  w1to9 %>% transmute(lfdn, wave = "W8",  pol_interest = clean_pol_interest(kp8_010)),
  w10   %>% transmute(lfdn, wave = "W10", pol_interest = clean_pol_interest(kp10_010)),
  w11   %>% transmute(lfdn, wave = "W11", pol_interest = clean_pol_interest(kp11_010)),
  w17   %>% transmute(lfdn, wave = "W17", pol_interest = clean_pol_interest(kp17_010)),
  w22   %>% transmute(lfdn, wave = "W22", pol_interest = clean_pol_interest(kp22_010)),
  w24   %>% transmute(lfdn, wave = "W24", pol_interest = clean_pol_interest(kp24_010)),
  w26   %>% transmute(lfdn, wave = "W26", pol_interest = clean_pol_interest(kp26_010)),
  w32   %>% transmute(lfdn, wave = "W32", pol_interest = clean_pol_interest(kp32_010))
)

gles_analysis <- gles_analysis %>%
  left_join(pol_interest_all, by = c("lfdn", "wave")) %>%
  mutate(
    post2022 = ifelse(wave %in% c("W22", "W24", "W26", "W32"), 1, 0)
  )

cat("Political interest and post-2022 indicator added.\n")


# ------------------------------------------------------------
# 3.9 Add demographic controls to gles_analysis
# ------------------------------------------------------------
# Demographic controls are respondent-level variables from the
# profile files. Age is calculated relative to the survey year.

gles_analysis <- gles_analysis %>%
  dplyr::left_join(profile_demographics_clean, by = "lfdn") %>%
  dplyr::mutate(
    age = year - birth_year,
    age = dplyr::if_else(age >= 18 & age <= 110, age, NA_real_)
  )

cat("Demographic controls added to gles_analysis.\n")

cat("Political interest missingness by wave:\n")
gles_analysis %>%
  group_by(wave) %>%
  summarise(
    n_total = n(),
    n_pol_interest_valid = sum(!is.na(pol_interest)),
    .groups = "drop"
  ) %>%
  print()

cat("Post-2022 indicator check:\n")
print(table(gles_analysis$wave, gles_analysis$post2022))


# ============================================================
# 4. BUILD PANEL DATASETS
# ============================================================
# The panel analysis compares respondents' support for increased
# defence spending in W17 with their later responses in W22, W24,
# W26, and W32.
#
# Party preference and political interest are measured at W17.
# The outcome is the change in defence spending support between
# W17 and the respective later wave.


# ------------------------------------------------------------
# 4.1 Identify respondents observed in both W17 and W22
# ------------------------------------------------------------
# The W17 baseline is first restricted to respondents who are
# also observed in W22. This reproduces the matched baseline
# used in the thesis analysis.

matched_w17_w22 <- w17$lfdn[w17$lfdn %in% w22$lfdn]


# ------------------------------------------------------------
# 4.2 Prepare the W17 baseline
# ------------------------------------------------------------
# This contains the pre-invasion values used in the panel models:
# defence support, party preference, and political interest.

w17_valid <- w17 %>%
  filter(lfdn %in% matched_w17_w22) %>%
  transmute(
    lfdn = lfdn,
    defence_w17 = clean_defence(kp17_2880y),
    party_w17 = clean_party(kp17_190ba),
    pol_interest_w17 = clean_pol_interest(kp17_010)
  ) %>%
  filter(!is.na(defence_w17), !is.na(party_w17))


# ------------------------------------------------------------
# 4.3 Build W17-W22 change dataset
# ------------------------------------------------------------

w22_valid <- w22 %>%
  filter(lfdn %in% w17_valid$lfdn) %>%
  transmute(
    lfdn = lfdn,
    defence_w22 = clean_defence(kp22_2880y)
  ) %>%
  filter(!is.na(defence_w22))

panel_w17_w22 <- inner_join(w17_valid, w22_valid, by = "lfdn") %>%
  mutate(change_w17_w22 = defence_w22 - defence_w17)


# ------------------------------------------------------------
# 4.4 Build W17-W24 change dataset
# ------------------------------------------------------------

w24_valid <- w24 %>%
  filter(lfdn %in% w17_valid$lfdn) %>%
  transmute(
    lfdn = lfdn,
    defence_w24 = clean_defence(kp24_2880y)
  ) %>%
  filter(!is.na(defence_w24))

panel_w17_w24 <- inner_join(w17_valid, w24_valid, by = "lfdn") %>%
  mutate(change_w17_w24 = defence_w24 - defence_w17)


# ------------------------------------------------------------
# 4.5 Build W17-W26 change dataset
# ------------------------------------------------------------

w26_valid <- w26 %>%
  filter(lfdn %in% w17_valid$lfdn) %>%
  transmute(
    lfdn = lfdn,
    defence_w26 = clean_defence(kp26_2880y)
  ) %>%
  filter(!is.na(defence_w26))

panel_w17_w26 <- inner_join(w17_valid, w26_valid, by = "lfdn") %>%
  mutate(change_w17_w26 = defence_w26 - defence_w17)


# ------------------------------------------------------------
# 4.6 Build W17-W32 change dataset
# ------------------------------------------------------------

w32_valid <- w32 %>%
  filter(lfdn %in% w17_valid$lfdn) %>%
  transmute(
    lfdn = lfdn,
    defence_w32 = clean_defence(kp32_2880y)
  ) %>%
  filter(!is.na(defence_w32))

panel_w17_w32 <- inner_join(w17_valid, w32_valid, by = "lfdn") %>%
  mutate(change_w17_w32 = defence_w32 - defence_w17)


# ------------------------------------------------------------
# 4.7 Check panel sample sizes
# ------------------------------------------------------------

cat("Panel datasets created.\n")

cat("Matched panel counts:\n")
cat("W17-W22:", nrow(panel_w17_w22), "\n")
cat("W17-W24:", nrow(panel_w17_w24), "\n")
cat("W17-W26:", nrow(panel_w17_w26), "\n")
cat("W17-W32:", nrow(panel_w17_w32), "\n")

cat("Valid W17 political interest counts:\n")
cat("W17-W22:", sum(!is.na(panel_w17_w22$pol_interest_w17)), "\n")
cat("W17-W24:", sum(!is.na(panel_w17_w24$pol_interest_w17)), "\n")
cat("W17-W26:", sum(!is.na(panel_w17_w26$pol_interest_w17)), "\n")
cat("W17-W32:", sum(!is.na(panel_w17_w32$pol_interest_w17)), "\n")

# ============================================================
# 5. DESCRIPTIVE STATISTICS — THESIS CHAPTER 5.1
# ============================================================
# This section creates the descriptive outputs used in Section 5.1:
# - Table 2: overall support for increased defence spending by wave
# - Table 3: mean support by party and wave
# - Figure 1: party-level trends in mean defence support, 2017-2025


# ------------------------------------------------------------
# 5.1 Overall support for increased defence spending by wave
# ------------------------------------------------------------
# Table 2 is based on all respondents with a valid response to
# the defence-spending item in each wave.

desc_sample_overview <- gles_long %>%
  dplyr::filter(!is.na(defence)) %>%
  dplyr::group_by(wave, year) %>%
  dplyr::summarise(
    N = dplyr::n(),
    mean_defence = round(mean(defence, na.rm = TRUE), 2),
    sd_defence = round(sd(defence, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year)

cat("Table 2 object created.\n")
print(desc_sample_overview)

table2_export <- desc_sample_overview %>%
  dplyr::mutate(year = as.character(year)) %>%
  dplyr::rename(
    Wave = wave,
    Year = year,
    N = N,
    `Mean defence support` = mean_defence,
    SD = sd_defence
  ) %>%
  as.data.frame()

stargazer(
  table2_export,
  summary = FALSE,
  type = "html",
  out = file.path(data_path, "table2_desc_sample_overview.html"),
  title = "Table 2. Overall Support for Increased Defence Spending by Wave",
  rownames = FALSE,
  digits = 2
)

cat("Table 2 saved as table2_desc_sample_overview.html\n")


# ------------------------------------------------------------
# 5.2 Mean support by party and wave
# ------------------------------------------------------------
# Table 3 reports the party-level means underlying the descriptive
# trend figure.

desc_by_party_wave <- gles_analysis %>%
  dplyr::filter(!is.na(party), !is.na(defence)) %>%
  dplyr::group_by(wave, year, party) %>%
  dplyr::summarise(
    mean_defence = round(mean(defence, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(year, party)

desc_mean_defence_by_party <- desc_by_party_wave %>%
  dplyr::select(wave, year, party, mean_defence) %>%
  tidyr::pivot_wider(
    names_from = party,
    values_from = mean_defence
  ) %>%
  dplyr::arrange(year) %>%
  dplyr::select(wave, year, `CDU/CSU`, SPD, Greens, FDP, Left, AfD, BSW)

cat("Table 3 object created.\n")
print(desc_mean_defence_by_party)

write.csv(
  desc_mean_defence_by_party,
  file.path(data_path, "table3_party_means.csv"),
  row.names = FALSE
)

table3_export <- desc_mean_defence_by_party %>%
  dplyr::mutate(year = as.character(year)) %>%
  dplyr::rename(
    Wave = wave,
    Year = year
  ) %>%
  as.data.frame()

stargazer(
  table3_export,
  summary = FALSE,
  type = "html",
  out = file.path(data_path, "table3_party_means.html"),
  title = "Table 3. Mean Support for Increased Defence Spending by Party and Wave",
  rownames = FALSE,
  digits = 2
)

cat("Table 3 saved as table3_party_means.csv and table3_party_means.html\n")


# ------------------------------------------------------------
# 5.3 Figure 1: Party trends in mean defence support
# ------------------------------------------------------------
# The x-axis uses approximate survey timing instead of only the
# calendar year. The dashed vertical line marks February 2022.

survey_x_lookup <- data.frame(
  wave = c("W8", "W10", "W11", "W17", "W22", "W24", "W26", "W32"),
  survey_x = c(
    2017 + 10/12,  # Sep/Oct 2017
    2018 + 11/12,  # Nov 2018
    2019 + 6/12,   # Jun 2019
    2021 + 7/12,   # Jul 2021
    2022 + 5/12,   # May 2022
    2023 + 5/12,   # May 2023
    2024 + 6/12,   # Jun 2024
    2025 + 5/12    # May 2025
  ),
  survey_label = c(
    "Sep 2017", "Nov 2018", "Jun 2019", "Jul 2021",
    "May 2022", "May 2023", "Jun 2024", "May 2025"
  )
)

desc_plot <- desc_by_party_wave %>%
  left_join(survey_x_lookup, by = "wave")

labels_plot <- desc_plot %>%
  filter(wave %in% c("W17", "W22"))

party_colours <- c(
  "CDU/CSU" = "#000000",
  "SPD"     = "#E3000F",
  "Greens"  = "#1AA037",
  "FDP"     = "#FFD700",
  "Left"    = "#D90199",
  "AfD"     = "#009EE0",
  "BSW"     = "#800080"
)

fig_descriptive <- ggplot(
  desc_plot,
  aes(x = survey_x, y = mean_defence, colour = party, group = party)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = 2022 + 2/12,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2022 + 2/12 + 0.05,
    y = 4.42,
    label = "February 2022",
    hjust = 0,
    size = 2.8,
    colour = "grey40"
  ) +
  geom_text_repel(
    data = labels_plot,
    aes(label = sprintf("%.2f", mean_defence)),
    size = 2.2,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.6,
    force = 3,
    force_pull = 0.1,
    segment.size = 0.25,
    segment.colour = "grey50",
    min.segment.length = 0.1,
    direction = "both"
  ) +
  scale_colour_manual(values = party_colours) +
  scale_y_continuous(
    limits = c(1.5, 4.6),
    breaks = c(1.5, 2, 2.5, 3, 3.5, 4, 4.5),
    name = "Mean support for increased defence spending (1-5)"
  ) +
  scale_x_continuous(
    breaks = survey_x_lookup$survey_x,
    labels = c(
      "Sep 2017", "Nov 2018", "Jun 2019", "Jul 2021",
      "May 2022", "May 2023", "Jun 2024", "May 2025"
    ),
    name = "Survey Wave"
  ) +
  labs(colour = "Party preference") +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(fig_descriptive)

ggsave(
  filename = file.path(data_path, "fig1_descriptive_party_trends.png"),
  plot = fig_descriptive,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Descriptive statistics section completed.\n")
cat("Figure 1 saved as fig1_descriptive_party_trends.png\n")


# ============================================================
# 6. WAVE-BY-WAVE OLS MODELS — THESIS CHAPTER 5.2
# ============================================================
# This section estimates the separate wave-by-wave OLS models
# reported in Section 5.2 and in Appendix Tables 4 and 5.
#
# The dependent variable is support for increased defence spending.
# CDU/CSU is the reference party category.
#
# Controls:
# - political interest
# - age
# - gender
# - education
# - East/West residence
#
# Outputs:
# - Figure 2: party coefficients relative to CDU/CSU across waves
# - Table 4: pre-invasion wave-specific OLS models
# - Table 5: post-invasion wave-specific OLS models


# ------------------------------------------------------------
# 6.1 Define wave-specific estimation samples
# ------------------------------------------------------------
# The wave_model_data() function keeps the observations used in
# the controlled wave-by-wave models. Respondents must have valid
# values on the controls added to the model.

wave_model_data <- function(w) {
  gles_analysis %>%
    dplyr::filter(
      wave == w,
      !is.na(pol_interest),
      !is.na(age),
      !is.na(gender),
      !is.na(education),
      !is.na(east_west)
    )
}


# ------------------------------------------------------------
# 6.2 Estimate wave-by-wave OLS models
# ------------------------------------------------------------
# lm_robust() estimates the same linear model as lm(), but directly
# reports heteroskedasticity-robust standard errors. HC2 is used
# throughout the wave-specific models.

ols_wave8 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W8"),
  se_type = "HC2"
)

ols_wave10 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W10"),
  se_type = "HC2"
)

ols_wave11 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W11"),
  se_type = "HC2"
)

ols_wave17 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W17"),
  se_type = "HC2"
)

ols_wave22 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W22"),
  se_type = "HC2"
)

ols_wave24 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W24"),
  se_type = "HC2"
)

ols_wave26 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W26"),
  se_type = "HC2"
)

ols_wave32 <- lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = wave_model_data("W32"),
  se_type = "HC2"
)

cat("Wave-by-wave OLS models estimated.\n")
cat("Complete-case counts by wave:\n")
cat("W8 :", nrow(wave_model_data("W8")), "\n")
cat("W10:", nrow(wave_model_data("W10")), "\n")
cat("W11:", nrow(wave_model_data("W11")), "\n")
cat("W17:", nrow(wave_model_data("W17")), "\n")
cat("W22:", nrow(wave_model_data("W22")), "\n")
cat("W24:", nrow(wave_model_data("W24")), "\n")
cat("W26:", nrow(wave_model_data("W26")), "\n")
cat("W32:", nrow(wave_model_data("W32")), "\n")

summary(ols_wave10)
summary(ols_wave8)
summary(ols_wave11)
summary(ols_wave17)


# ------------------------------------------------------------
# 6.3 Prepare party coefficients for Figure 2
# ------------------------------------------------------------
# Figure 2 plots only the party coefficients from the wave-specific
# OLS models. These coefficients show each party's difference from
# CDU/CSU in the same wave, net of the controls.

extract_party_coefficients <- function(model, wave_label, survey_x) {
  broom::tidy(model, conf.int = TRUE) %>%
    dplyr::filter(grepl("^party", term)) %>%
    dplyr::mutate(
      party = gsub("^party", "", term),
      wave = wave_label,
      survey_x = survey_x
    ) %>%
    dplyr::select(wave, survey_x, party, estimate, conf.low, conf.high)
}

all_ols_coefficients <- bind_rows(
  extract_party_coefficients(ols_wave8,  "W8",  2017 + 10/12),
  extract_party_coefficients(ols_wave10, "W10", 2018 + 11/12),
  extract_party_coefficients(ols_wave11, "W11", 2019 +  6/12),
  extract_party_coefficients(ols_wave17, "W17", 2021 +  7/12),
  extract_party_coefficients(ols_wave22, "W22", 2022 +  5/12),
  extract_party_coefficients(ols_wave24, "W24", 2023 +  5/12),
  extract_party_coefficients(ols_wave26, "W26", 2024 +  6/12),
  extract_party_coefficients(ols_wave32, "W32", 2025 +  5/12)
)

all_ols_coefficients <- all_ols_coefficients %>%
  mutate(
    survey_label = case_when(
      wave == "W8"  ~ "Sep 2017",
      wave == "W10" ~ "Nov 2018",
      wave == "W11" ~ "Jun 2019",
      wave == "W17" ~ "Jul 2021",
      wave == "W22" ~ "May 2022",
      wave == "W24" ~ "May 2023",
      wave == "W26" ~ "Jun 2024",
      wave == "W32" ~ "May 2025"
    ),
    party = factor(
      party,
      levels = c("SPD", "Greens", "FDP", "Left", "AfD", "BSW")
    )
  )


# ------------------------------------------------------------
# 6.4 Figure 2: OLS coefficients relative to CDU/CSU
# ------------------------------------------------------------

party_colours_no_cdu <- c(
  "SPD"    = "#E3000F",
  "Greens" = "#1AA037",
  "FDP"    = "#FFD700",
  "Left"   = "#D90199",
  "AfD"    = "#009EE0",
  "BSW"    = "#800080"
)

fig_ols_coefficients <- ggplot(
  all_ols_coefficients,
  aes(x = survey_x, y = estimate, colour = party, group = party)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = party),
    alpha = 0.12,
    colour = NA
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 2022 + 2/12,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = 2022 + 2/12 + 0.05,
    y = 0.15,
    label = "February 2022",
    hjust = 0,
    size = 2.8,
    colour = "grey40"
  ) +
  scale_colour_manual(values = party_colours_no_cdu) +
  scale_fill_manual(values = party_colours_no_cdu) +
  scale_y_continuous(
    name = "Coefficient relative to CDU/CSU (defence support, 1-5)",
    breaks = seq(-1.6, 0.2, by = 0.2)
  ) +
  scale_x_continuous(
    breaks = c(
      2017 + 10/12, 2018 + 11/12, 2019 + 6/12,
      2021 + 7/12, 2022 + 5/12, 2023 + 5/12,
      2024 + 6/12, 2025 + 5/12
    ),
    labels = c(
      "Sep 2017", "Nov 2018", "Jun 2019", "Jul 2021",
      "May 2022", "May 2023", "Jun 2024", "May 2025"
    ),
    name = "Survey Wave"
  ) +
  labs(
    colour = "Party preference",
    fill = "Party preference"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(fig_ols_coefficients)

ggsave(
  filename = file.path(data_path, "fig2_ols_coefficients.png"),
  plot = fig_ols_coefficients,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Figure 2 saved as fig2_ols_coefficients.png\n")

all_ols_coefficients_df <- as.data.frame(all_ols_coefficients)
print(all_ols_coefficients_df, row.names = FALSE)


# ------------------------------------------------------------
# 6.5 Prepare Table 4 and Table 5 exports
# ------------------------------------------------------------
# stargazer works more smoothly with regular lm objects.
# The lm models below use the same formulas and samples as the
# lm_robust models above. The robust HC2 standard errors from
# lm_robust are then supplied manually to stargazer.

table_wave_data <- function(w) {
  wave_model_data(w) %>%
    dplyr::filter(
      !is.na(defence),
      !is.na(party),
      !is.na(pol_interest)
    ) %>%
    droplevels() %>%
    as.data.frame()
}

lm_wave8 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W8")
)

lm_wave10 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W10")
)

lm_wave11 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W11")
)

lm_wave17 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W17")
)

lm_wave22 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W22")
)

lm_wave24 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W24")
)

lm_wave26 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W26")
)

lm_wave32 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = table_wave_data("W32")
)

se_w8  <- robust_se(ols_wave8)
se_w10 <- robust_se(ols_wave10)
se_w11 <- robust_se(ols_wave11)
se_w17 <- robust_se(ols_wave17)
se_w22 <- robust_se(ols_wave22)
se_w24 <- robust_se(ols_wave24)
se_w26 <- robust_se(ols_wave26)
se_w32 <- robust_se(ols_wave32)


# ------------------------------------------------------------
# 6.6 Check that table models match the main wave models
# ------------------------------------------------------------

cat("Table 4/5 export sample sizes:\n")
cat("W8 :", nrow(table_wave_data("W8")),  " | model N =", nobs(lm_wave8),  "\n")
cat("W10:", nrow(table_wave_data("W10")), " | model N =", nobs(lm_wave10), "\n")
cat("W11:", nrow(table_wave_data("W11")), " | model N =", nobs(lm_wave11), "\n")
cat("W17:", nrow(table_wave_data("W17")), " | model N =", nobs(lm_wave17), "\n")
cat("W22:", nrow(table_wave_data("W22")), " | model N =", nobs(lm_wave22), "\n")
cat("W24:", nrow(table_wave_data("W24")), " | model N =", nobs(lm_wave24), "\n")
cat("W26:", nrow(table_wave_data("W26")), " | model N =", nobs(lm_wave26), "\n")
cat("W32:", nrow(table_wave_data("W32")), " | model N =", nobs(lm_wave32), "\n")

stopifnot(
  nrow(table_wave_data("W8"))  == nobs(lm_wave8),
  nrow(table_wave_data("W10")) == nobs(lm_wave10),
  nrow(table_wave_data("W11")) == nobs(lm_wave11),
  nrow(table_wave_data("W17")) == nobs(lm_wave17),
  nrow(table_wave_data("W22")) == nobs(lm_wave22),
  nrow(table_wave_data("W24")) == nobs(lm_wave24),
  nrow(table_wave_data("W26")) == nobs(lm_wave26),
  nrow(table_wave_data("W32")) == nobs(lm_wave32)
)

cat("Coefficient equivalence checks:\n")
print(all.equal(unname(coef(lm_wave8)),  unname(coef(ols_wave8)),  tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave10)), unname(coef(ols_wave10)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave11)), unname(coef(ols_wave11)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave17)), unname(coef(ols_wave17)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave22)), unname(coef(ols_wave22)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave24)), unname(coef(ols_wave24)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave26)), unname(coef(ols_wave26)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_wave32)), unname(coef(ols_wave32)), tolerance = 1e-10))


# ------------------------------------------------------------
# 6.7 Table 4: Pre-invasion wave-specific OLS models
# ------------------------------------------------------------

stargazer(
  lm_wave8, lm_wave10, lm_wave11, lm_wave17,
  se = list(se_w8, se_w10, se_w11, se_w17),
  type = "text",
  title = "Table 4. OLS Regression: Defence Spending Support (Pre-Invasion Waves) with Demographic Controls",
  column.labels = c("W8 Sep 2017", "W10 Nov 2018", "W11 Jun 2019", "W17 Jul 2021"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the reference party category. Male, low education, and West Germany are the demographic reference categories. Political interest is reverse-coded so that higher values indicate greater political interest. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)

stargazer(
  lm_wave8, lm_wave10, lm_wave11, lm_wave17,
  se = list(se_w8, se_w10, se_w11, se_w17),
  type = "html",
  out = file.path(data_path, "table4_ols_pre_invasion_controls.html"),
  title = "Table 4. OLS Regression: Defence Spending Support (Pre-Invasion Waves) with Demographic Controls",
  column.labels = c("W8 Sep 2017", "W10 Nov 2018", "W11 Jun 2019", "W17 Jul 2021"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the reference party category. Male, low education, and West Germany are the demographic reference categories. Political interest is reverse-coded so that higher values indicate greater political interest. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)


# ------------------------------------------------------------
# 6.8 Table 5: Post-invasion wave-specific OLS models
# ------------------------------------------------------------

stargazer(
  lm_wave22, lm_wave24, lm_wave26, lm_wave32,
  se = list(se_w22, se_w24, se_w26, se_w32),
  type = "text",
  title = "Table 5. OLS Regression: Defence Spending Support (Post-Invasion Waves) with Demographic Controls",
  column.labels = c("W22 May 2022", "W24 May 2023", "W26 Jun 2024", "W32 May 2025"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the reference party category. Male, low education, and West Germany are the demographic reference categories. Political interest is reverse-coded so that higher values indicate greater political interest. BSW was not available prior to 2024. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)

stargazer(
  lm_wave22, lm_wave24, lm_wave26, lm_wave32,
  se = list(se_w22, se_w24, se_w26, se_w32),
  type = "html",
  out = file.path(data_path, "table5_ols_post_invasion_controls.html"),
  title = "Table 5. OLS Regression: Defence Spending Support (Post-Invasion Waves) with Demographic Controls",
  column.labels = c("W22 May 2022", "W24 May 2023", "W26 Jun 2024", "W32 May 2025"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the reference party category. Male, low education, and West Germany are the demographic reference categories. Political interest is reverse-coded so that higher values indicate greater political interest. BSW was not available prior to 2024. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)

cat("Wave-by-wave OLS section completed.\n")
cat("Table 4 saved as table4_ols_pre_invasion_controls.html\n")
cat("Table 5 saved as table5_ols_post_invasion_controls.html\n")


# ============================================================
# 7. POOLED INTERACTION MODEL — THESIS CHAPTER 5.3
# ============================================================
# This section tests whether partisan differences in support for
# increased defence spending changed after February 2022.
#
# The data are pooled across all eight waves. BSW is excluded here
# because it has no pre-2022 observations and therefore cannot be
# compared before and after the invasion.
#
# CDU/CSU is the reference party category.
# Standard errors are clustered by respondent ID because some
# respondents appear in more than one wave.

# ------------------------------------------------------------
# 7.1 Build pooled dataset for H2
# ------------------------------------------------------------

gles_h2 <- gles_analysis %>%
  dplyr::filter(
    party %in% c("CDU/CSU", "SPD", "FDP", "Greens", "Left", "AfD"),
    !is.na(pol_interest),
    !is.na(age),
    !is.na(gender),
    !is.na(education),
    !is.na(east_west)
  ) %>%
  dplyr::mutate(
    party = factor(
      party,
      levels = c("CDU/CSU", "SPD", "FDP", "Greens", "Left", "AfD")
    ),
    post2022 = ifelse(wave %in% c("W22", "W24", "W26", "W32"), 1, 0),
    lfdn_int = as.integer(lfdn)
  )

cat("Pooled H2 dataset created successfully.\n")
cat("Observations in pooled H2 sample:", nrow(gles_h2), "\n")
cat("Unique respondents in pooled H2 sample:", length(unique(gles_h2$lfdn_int)), "\n")


# ------------------------------------------------------------
# 7.2 Estimate pooled OLS models with clustered standard errors
# ------------------------------------------------------------

# Model 1: party preference, post-2022 period, and controls
# Model 2: adds party x post-2022 interactions

ols_h2_main_effects_clustered <- estimatr::lm_robust(
  defence ~ party + post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2,
  clusters = lfdn_int,
  se_type = "CR0"
)

ols_h2_pooled_clustered <- estimatr::lm_robust(
  defence ~ party * post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2,
  clusters = lfdn_int,
  se_type = "CR0"
)

cat("Pooled clustered OLS models estimated successfully.\n")

# ------------------------------------------------------------
# 7.3 Estimate lm versions for stargazer and joint F-test
# ------------------------------------------------------------
# The lm objects are used for table export and for the joint F-test.
# The clustered standard errors are supplied separately.

lm_h2_main_effects <- lm(
  defence ~ party + post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2
)

pooled_interaction_lm <- lm(
  defence ~ party * post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2
)


# ------------------------------------------------------------
# 7.4 Joint F-test of the interaction terms
# ------------------------------------------------------------
# This tests whether the party x post-2022 interaction terms are
# jointly different from zero.

vcov_clustered <- sandwich::vcovCL(
  pooled_interaction_lm,
  cluster = ~ lfdn_int,
  data = gles_h2
)

joint_test_h2_clustered <- car::linearHypothesis(
  pooled_interaction_lm,
  c(
    "partySPD:post2022 = 0",
    "partyFDP:post2022 = 0",
    "partyGreens:post2022 = 0",
    "partyLeft:post2022 = 0",
    "partyAfD:post2022 = 0"
  ),
  vcov. = vcov_clustered,
  test = "F"
)

cat("Joint F-test for pooled interaction terms completed.\n")
print(joint_test_h2_clustered)


# ------------------------------------------------------------
# 7.5 Export Table 6
# ------------------------------------------------------------

se_h2_main_effects_clustered <- robust_se(ols_h2_main_effects_clustered)
se_h2_pooled_clustered <- robust_se(ols_h2_pooled_clustered)

stargazer(
  lm_h2_main_effects,
  pooled_interaction_lm,
  se = list(se_h2_main_effects_clustered, se_h2_pooled_clustered),
  type = "text",
  title = "Table 6. Pooled Model: Party Preference, Post-2022 Shift, and Demographic Controls",
  column.labels = c("Main effects + controls", "With interactions + controls"),
  dep.var.labels = "Defence spending support (1-5)",
  covariate.labels = c(
    "SPD",
    "FDP",
    "Greens",
    "Left",
    "AfD",
    "Post-2022",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany",
    "SPD x Post-2022",
    "FDP x Post-2022",
    "Greens x Post-2022",
    "Left x Post-2022",
    "AfD x Post-2022"
  ),
  omit = "Constant",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Clustered standard errors (CR0) by respondent ID in parentheses.",
    "Models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the party reference category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "BSW is excluded because no pre-2022 observations are available.",
    "Positive interaction = convergence toward CDU/CSU post-2022.",
    "Negative interaction = divergence from CDU/CSU post-2022."
  ),
  notes.append = FALSE
)

stargazer(
  lm_h2_main_effects,
  pooled_interaction_lm,
  se = list(se_h2_main_effects_clustered, se_h2_pooled_clustered),
  type = "html",
  out = file.path(data_path, "table6_pooled_interaction_h2_clustered.html"),
  title = "Table 6. Pooled Model: Party Preference, Post-2022 Shift, and Demographic Controls",
  column.labels = c("Main effects + controls", "With interactions + controls"),
  dep.var.labels = "Defence spending support (1-5)",
  covariate.labels = c(
    "SPD",
    "FDP",
    "Greens",
    "Left",
    "AfD",
    "Post-2022",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany",
    "SPD x Post-2022",
    "FDP x Post-2022",
    "Greens x Post-2022",
    "Left x Post-2022",
    "AfD x Post-2022"
  ),
  omit = "Constant",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Clustered standard errors (CR0) by respondent ID in parentheses.",
    "Models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the party reference category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "BSW is excluded because no pre-2022 observations are available.",
    "Positive interaction = convergence toward CDU/CSU post-2022.",
    "Negative interaction = divergence from CDU/CSU post-2022."
  ),
  notes.append = FALSE
)

cat("Table 6 saved as table6_pooled_interaction_h2_clustered.html\n")


# ------------------------------------------------------------
# 7.6 Figure 3: Post-2022 changes in party gaps
# ------------------------------------------------------------

interaction_coefficients_clustered <- broom::tidy(
  ols_h2_pooled_clustered,
  conf.int = TRUE
) %>%
  dplyr::filter(grepl(":post2022", term)) %>%
  dplyr::mutate(
    party = dplyr::case_when(
      term == "partySPD:post2022"    ~ "SPD",
      term == "partyFDP:post2022"    ~ "FDP",
      term == "partyGreens:post2022" ~ "Greens",
      term == "partyLeft:post2022"   ~ "Left",
      term == "partyAfD:post2022"    ~ "AfD"
    ),
    party = factor(
      party,
      levels = c("AfD", "Left", "FDP", "Greens", "SPD")
    )
  )

fig_interaction_effects <- ggplot(
  interaction_coefficients_clustered,
  aes(x = estimate, y = party, colour = party)
) +
  geom_segment(
    aes(x = conf.low, xend = conf.high, y = party, yend = party),
    linewidth = 1.2
  ) +
  geom_point(size = 4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  scale_colour_manual(values = c(
    "SPD"    = "#E3000F",
    "Greens" = "#1AA037",
    "FDP"    = "#FFD700",
    "Left"   = "#D90199",
    "AfD"    = "#009EE0"
  )) +
  scale_x_continuous(
    limits = c(-0.7, 0.7),
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6),
    name = "Post-2022 change in gap relative to CDU/CSU\n(positive = convergence, negative = divergence)"
  ) +
  labs(y = NULL) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.text.y = element_text(size = 11)
  )

print(fig_interaction_effects)

ggsave(
  filename = file.path(data_path, "fig3_interaction_effects_clustered.png"),
  plot = fig_interaction_effects,
  width = 7,
  height = 3,
  dpi = 300
)


cat("Figure 3 saved as fig3_interaction_effects_clustered.png\n")
cat("Pooled interaction model section completed.\n")



# ============================================================
# 8. PANEL CHANGE-SCORE ANALYSIS — THESIS SECTION 5.4
# ============================================================
# This section estimates the within-person change models reported
# in Table 7.
#
# The baseline is W17, the last pre-invasion wave. The dependent
# variable is the change in support for increased defence spending
# between W17 and each later wave.
#
# Party preference, political interest, and demographic controls
# are measured at baseline. CDU/CSU is the reference category.


# ------------------------------------------------------------
# 8.1 Add W17 baseline demographic controls
# ------------------------------------------------------------
# The panel datasets were built earlier in Section 4. Here they are
# supplemented with baseline demographic controls and a cleaned
# W17 party factor for the regression models.

make_panel_party_factor <- function(x) {
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  
  x_lab <- case_when(
    x_num == 1   ~ "CDU/CSU",
    x_num == 4   ~ "SPD",
    x_num == 6   ~ "Greens",
    x_num == 5   ~ "FDP",
    x_num == 7   ~ "Left",
    x_num == 322 ~ "AfD",
    TRUE         ~ as.character(x)
  )
  
  factor(
    x_lab,
    levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
  )
}

panel_baseline_controls <- profile_demographics_clean %>%
  dplyr::transmute(
    lfdn,
    gender,
    education,
    east_west,
    age_w17 = dplyr::case_when(
      !is.na(birth_year) ~ 2021 - birth_year,
      TRUE ~ NA_real_
    )
  )

panel_w17_w22 <- panel_w17_w22 %>%
  dplyr::select(-dplyr::any_of(c("gender", "education", "east_west", "age_w17"))) %>%
  dplyr::left_join(panel_baseline_controls, by = "lfdn") %>%
  dplyr::mutate(party_w17 = make_panel_party_factor(party_w17))

panel_w17_w24 <- panel_w17_w24 %>%
  dplyr::select(-dplyr::any_of(c("gender", "education", "east_west", "age_w17"))) %>%
  dplyr::left_join(panel_baseline_controls, by = "lfdn") %>%
  dplyr::mutate(party_w17 = make_panel_party_factor(party_w17))

panel_w17_w26 <- panel_w17_w26 %>%
  dplyr::select(-dplyr::any_of(c("gender", "education", "east_west", "age_w17"))) %>%
  dplyr::left_join(panel_baseline_controls, by = "lfdn") %>%
  dplyr::mutate(party_w17 = make_panel_party_factor(party_w17))

panel_w17_w32 <- panel_w17_w32 %>%
  dplyr::select(-dplyr::any_of(c("gender", "education", "east_west", "age_w17"))) %>%
  dplyr::left_join(panel_baseline_controls, by = "lfdn") %>%
  dplyr::mutate(party_w17 = make_panel_party_factor(party_w17))

cat("Panel datasets prepared with W17 baseline demographic controls.\n")
cat("Panel sample sizes after party recoding:\n")
cat("W17-W22:", nrow(panel_w17_w22), "\n")
cat("W17-W24:", nrow(panel_w17_w24), "\n")
cat("W17-W26:", nrow(panel_w17_w26), "\n")
cat("W17-W32:", nrow(panel_w17_w32), "\n")

cat("Baseline-control coverage in panel datasets:\n")
cat("W17-W22 complete on controls:",
    sum(!is.na(panel_w17_w22$age_w17) & !is.na(panel_w17_w22$gender) &
          !is.na(panel_w17_w22$education) & !is.na(panel_w17_w22$east_west)), "\n")
cat("W17-W24 complete on controls:",
    sum(!is.na(panel_w17_w24$age_w17) & !is.na(panel_w17_w24$gender) &
          !is.na(panel_w17_w24$education) & !is.na(panel_w17_w24$east_west)), "\n")
cat("W17-W26 complete on controls:",
    sum(!is.na(panel_w17_w26$age_w17) & !is.na(panel_w17_w26$gender) &
          !is.na(panel_w17_w26$education) & !is.na(panel_w17_w26$east_west)), "\n")
cat("W17-W32 complete on controls:",
    sum(!is.na(panel_w17_w32$age_w17) & !is.na(panel_w17_w32$gender) &
          !is.na(panel_w17_w32$education) & !is.na(panel_w17_w32$east_west)), "\n")


# ------------------------------------------------------------
# 8.2 Estimate change-score models
# ------------------------------------------------------------
# Outcome: change in defence-spending support since W17.
# Predictors: W17 party preference, W17 political interest, and
# W17 baseline demographic controls.

panel_model_data <- function(df) {
  df %>%
    dplyr::filter(
      !is.na(pol_interest_w17),
      !is.na(age_w17),
      !is.na(gender),
      !is.na(education),
      !is.na(east_west)
    )
}

ols_change_w17_w22 <- lm_robust(
  change_w17_w22 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = panel_model_data(panel_w17_w22),
  se_type = "HC2"
)

ols_change_w17_w24 <- lm_robust(
  change_w17_w24 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = panel_model_data(panel_w17_w24),
  se_type = "HC2"
)

ols_change_w17_w26 <- lm_robust(
  change_w17_w26 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = panel_model_data(panel_w17_w26),
  se_type = "HC2"
)

ols_change_w17_w32 <- lm_robust(
  change_w17_w32 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = panel_model_data(panel_w17_w32),
  se_type = "HC2"
)

cat("Change-score models estimated with W17 baseline controls.\n")
cat("Model observation counts:\n")
cat("W17-W22:", nobs(ols_change_w17_w22), "\n")
cat("W17-W24:", nobs(ols_change_w17_w24), "\n")
cat("W17-W26:", nobs(ols_change_w17_w26), "\n")
cat("W17-W32:", nobs(ols_change_w17_w32), "\n")


# ------------------------------------------------------------
# 8.3 Descriptive change tables
# ------------------------------------------------------------
# Mean within-person change by W17 party preference.
# These descriptive tables use the same complete-case samples as
# the controlled panel regression models.

desc_change_by_party_long <- bind_rows(
  panel_model_data(panel_w17_w22) %>%
    transmute(
      comparison = "W17-W22",
      party_w17,
      change = change_w17_w22
    ),
  
  panel_model_data(panel_w17_w24) %>%
    transmute(
      comparison = "W17-W24",
      party_w17,
      change = change_w17_w24
    ),
  
  panel_model_data(panel_w17_w26) %>%
    transmute(
      comparison = "W17-W26",
      party_w17,
      change = change_w17_w26
    ),
  
  panel_model_data(panel_w17_w32) %>%
    transmute(
      comparison = "W17-W32",
      party_w17,
      change = change_w17_w32
    )
) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c("W17-W22", "W17-W24", "W17-W26", "W17-W32")
    ),
    party_w17 = factor(
      party_w17,
      levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
    )
  )

desc_change_by_party_summary <- desc_change_by_party_long %>%
  group_by(comparison, party_w17) %>%
  summarise(
    n = n(),
    mean_change = round(mean(change, na.rm = TRUE), 3),
    sd_change = round(sd(change, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(comparison, party_w17)

desc_change_by_party_wide <- desc_change_by_party_summary %>%
  dplyr::select(comparison, party_w17, mean_change) %>%
  tidyr::pivot_wider(
    names_from = party_w17,
    values_from = mean_change
  )

cat("Descriptive change tables created.\n")

write.csv(
  desc_change_by_party_summary,
  file.path(data_path, "desc_change_by_party_summary.csv"),
  row.names = FALSE
)

write.csv(
  desc_change_by_party_wide,
  file.path(data_path, "desc_change_by_party_wide.csv"),
  row.names = FALSE
)

cat("Descriptive change tables saved as CSV files.\n")


# ------------------------------------------------------------
# 8.4 Rebuild plain lm models for Table 7 export
# ------------------------------------------------------------
# Plain lm objects are used for the stargazer table. They use the
# same formulas and complete-case samples as the lm_robust models.

table7_data_w22 <- panel_model_data(panel_w17_w22) %>%
  dplyr::select(
    change_w17_w22, party_w17, pol_interest_w17,
    age_w17, gender, education, east_west
  ) %>%
  as.data.frame()

table7_data_w24 <- panel_model_data(panel_w17_w24) %>%
  dplyr::select(
    change_w17_w24, party_w17, pol_interest_w17,
    age_w17, gender, education, east_west
  ) %>%
  as.data.frame()

table7_data_w26 <- panel_model_data(panel_w17_w26) %>%
  dplyr::select(
    change_w17_w26, party_w17, pol_interest_w17,
    age_w17, gender, education, east_west
  ) %>%
  as.data.frame()

table7_data_w32 <- panel_model_data(panel_w17_w32) %>%
  dplyr::select(
    change_w17_w32, party_w17, pol_interest_w17,
    age_w17, gender, education, east_west
  ) %>%
  as.data.frame()

table7_data_w22$party_w17 <- factor(
  table7_data_w22$party_w17,
  levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
)

table7_data_w24$party_w17 <- factor(
  table7_data_w24$party_w17,
  levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
)

table7_data_w26$party_w17 <- factor(
  table7_data_w26$party_w17,
  levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
)

table7_data_w32$party_w17 <- factor(
  table7_data_w32$party_w17,
  levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
)

lm_table7_w22 <- lm(
  change_w17_w22 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = table7_data_w22
)

lm_table7_w24 <- lm(
  change_w17_w24 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = table7_data_w24
)

lm_table7_w26 <- lm(
  change_w17_w26 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = table7_data_w26
)

lm_table7_w32 <- lm(
  change_w17_w32 ~ party_w17 + pol_interest_w17 + age_w17 + gender + education + east_west,
  data = table7_data_w32
)

cat("Table 7 lm models rebuilt successfully.\n")
cat("N W17-W22:", nobs(lm_table7_w22), "\n")
cat("N W17-W24:", nobs(lm_table7_w24), "\n")
cat("N W17-W26:", nobs(lm_table7_w26), "\n")
cat("N W17-W32:", nobs(lm_table7_w32), "\n")


# ------------------------------------------------------------
# 8.5 Compute HC2 robust standard errors for Table 7
# ------------------------------------------------------------

se_table7_w22 <- sqrt(diag(vcovHC(lm_table7_w22, type = "HC2")))
se_table7_w24 <- sqrt(diag(vcovHC(lm_table7_w24, type = "HC2")))
se_table7_w26 <- sqrt(diag(vcovHC(lm_table7_w26, type = "HC2")))
se_table7_w32 <- sqrt(diag(vcovHC(lm_table7_w32, type = "HC2")))

cat("Robust HC2 standard errors for Table 7 created.\n")


# ------------------------------------------------------------
# 8.6 Export Table 7
# ------------------------------------------------------------

stargazer(
  lm_table7_w22, lm_table7_w24, lm_table7_w26, lm_table7_w32,
  se = list(se_table7_w22, se_table7_w24, se_table7_w26, se_table7_w32),
  type = "html",
  out = file.path(data_path, "table7_change_models.html"),
  title = "Table 7. Change in Defence Spending Support with Baseline Controls",
  column.labels = c("W17-W22", "W17-W24", "W17-W26", "W17-W32"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest (W17)",
    "Age (W17)",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels.include = FALSE,
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include baseline controls for age, gender, education, and East/West residence. CDU/CSU is the reference party category. Male, low education, and West Germany are the demographic reference categories. Coefficients indicate how much more or less each party group's support changed relative to CDU/CSU supporters between W17 and the respective follow-up wave. Political interest is measured in W17.",
  notes.append = FALSE
)

cat("Table 7 saved as table7_change_models.html\n")
cat("Panel change-score analysis section completed.\n")

# ============================================================
# 9. ROBUSTNESS CHECKS — THESIS CHAPTER 5.5 
# ============================================================
# These checks assess whether the main findings depend on specific
# modelling or measurement choices.


# ------------------------------------------------------------
# 9.1 Ordered-logit robustness check
# ------------------------------------------------------------
# The defence-spending item is ordinal. The main analysis uses OLS
# for interpretability, while this robustness check re-estimates
# the same controlled specifications using ordered logit.


# ------------------------------------------------------------
# 9.1.1 Build ordered-logit datasets
# ------------------------------------------------------------

ologit_wave_data <- function(w) {
  table_wave_data(w) %>%
    dplyr::mutate(
      defence_ord = ordered(defence, levels = c(1, 2, 3, 4, 5))
    ) %>%
    droplevels() %>%
    as.data.frame()
}

ologit_h2_data <- gles_analysis %>%
  dplyr::filter(
    wave %in% c("W8", "W10", "W11", "W17", "W22", "W24", "W26", "W32"),
    !is.na(defence),
    !is.na(party),
    !is.na(pol_interest),
    !is.na(age),
    !is.na(gender),
    !is.na(education),
    !is.na(east_west),
    party %in% c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
  ) %>%
  dplyr::mutate(
    defence_ord = ordered(defence, levels = c(1, 2, 3, 4, 5)),
    post2022 = ifelse(wave %in% c("W22", "W24", "W26", "W32"), 1, 0),
    party = factor(
      party,
      levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD")
    )
  ) %>%
  droplevels() %>%
  as.data.frame()

cat("Ordered-logit datasets created.\n")
cat("Pooled ordered-logit H2 sample size:", nrow(ologit_h2_data), "\n")


# ------------------------------------------------------------
# 9.1.2 Estimate wave-by-wave ordered-logit models
# ------------------------------------------------------------

ologit_wave8 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W8"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave10 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W10"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave11 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W11"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave17 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W17"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave22 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W22"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave24 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W24"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave26 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W26"),
  Hess = TRUE,
  method = "logistic"
)

ologit_wave32 <- MASS::polr(
  defence_ord ~ party + pol_interest + age + gender + education + east_west,
  data = ologit_wave_data("W32"),
  Hess = TRUE,
  method = "logistic"
)

cat("Wave-by-wave ordered-logit models estimated.\n")


# ------------------------------------------------------------
# 9.1.3 Estimate pooled ordered-logit interaction model
# ------------------------------------------------------------

ologit_h2_pooled <- MASS::polr(
  defence_ord ~ party * post2022 + pol_interest + age + gender + education + east_west,
  data = ologit_h2_data,
  Hess = TRUE,
  method = "logistic"
)

cat("Pooled ordered-logit interaction model estimated.\n")
cat("Pooled ordered-logit N:", nobs(ologit_h2_pooled), "\n")


# ------------------------------------------------------------
# 9.1.4 Extract ordered-logit results
# ------------------------------------------------------------
# polr() does not report p-values by default. Approximate two-sided
# p-values are calculated from the reported t values.

extract_polr_results <- function(model, model_name) {
  ct <- coef(summary(model))
  p_vals <- 2 * pnorm(abs(ct[, "t value"]), lower.tail = FALSE)
  
  out <- data.frame(
    model = model_name,
    term = rownames(ct),
    estimate = ct[, "Value"],
    std.error = ct[, "Std. Error"],
    statistic = ct[, "t value"],
    p.value = p_vals,
    stringsAsFactors = FALSE
  )
  
  out$stars <- ifelse(out$p.value < 0.01, "***",
                      ifelse(out$p.value < 0.05, "**",
                             ifelse(out$p.value < 0.10, "*", "")))
  
  out$estimate_display <- paste0(sprintf("%.3f", out$estimate), out$stars)
  out$se_display <- paste0("(", sprintf("%.3f", out$std.error), ")")
  
  return(out)
}

ologit_wave_results <- dplyr::bind_rows(
  extract_polr_results(ologit_wave8,  "W8"),
  extract_polr_results(ologit_wave10, "W10"),
  extract_polr_results(ologit_wave11, "W11"),
  extract_polr_results(ologit_wave17, "W17"),
  extract_polr_results(ologit_wave22, "W22"),
  extract_polr_results(ologit_wave24, "W24"),
  extract_polr_results(ologit_wave26, "W26"),
  extract_polr_results(ologit_wave32, "W32")
)

ologit_wave_results_core <- ologit_wave_results %>%
  dplyr::filter(!grepl("\\|", term))

ologit_h2_results <- extract_polr_results(ologit_h2_pooled, "Pooled H2") %>%
  dplyr::filter(!grepl("\\|", term))

ologit_h2_results_or <- ologit_h2_results %>%
  dplyr::mutate(
    odds_ratio = exp(estimate),
    odds_ratio_display = sprintf("%.3f", odds_ratio)
  )


# ------------------------------------------------------------
# 9.1.5 Check ordered-logit samples
# ------------------------------------------------------------

cat("Wave-specific ordered-logit observation counts:\n")
cat("W8 :", nobs(ologit_wave8),  "\n")
cat("W10:", nobs(ologit_wave10), "\n")
cat("W11:", nobs(ologit_wave11), "\n")
cat("W17:", nobs(ologit_wave17), "\n")
cat("W22:", nobs(ologit_wave22), "\n")
cat("W24:", nobs(ologit_wave24), "\n")
cat("W26:", nobs(ologit_wave26), "\n")
cat("W32:", nobs(ologit_wave32), "\n")
cat("Pooled H2:", nobs(ologit_h2_pooled), "\n")

stopifnot(
  nrow(ologit_wave_data("W8"))  == nobs(ologit_wave8),
  nrow(ologit_wave_data("W10")) == nobs(ologit_wave10),
  nrow(ologit_wave_data("W11")) == nobs(ologit_wave11),
  nrow(ologit_wave_data("W17")) == nobs(ologit_wave17),
  nrow(ologit_wave_data("W22")) == nobs(ologit_wave22),
  nrow(ologit_wave_data("W24")) == nobs(ologit_wave24),
  nrow(ologit_wave_data("W26")) == nobs(ologit_wave26),
  nrow(ologit_wave_data("W32")) == nobs(ologit_wave32)
)

stopifnot(
  nobs(ologit_wave8)  == nobs(lm_wave8),
  nobs(ologit_wave10) == nobs(lm_wave10),
  nobs(ologit_wave11) == nobs(lm_wave11),
  nobs(ologit_wave17) == nobs(lm_wave17),
  nobs(ologit_wave22) == nobs(lm_wave22),
  nobs(ologit_wave24) == nobs(lm_wave24),
  nobs(ologit_wave26) == nobs(lm_wave26),
  nobs(ologit_wave32) == nobs(lm_wave32)
)

cat("Ordered-logit sample checks passed.\n")


# ------------------------------------------------------------
# 9.1.6 Save ordered-logit output files
# ------------------------------------------------------------

write.csv(
  ologit_wave_results_core,
  file.path(data_path, "table_ologit_wave_results_core.csv"),
  row.names = FALSE
)

write.csv(
  ologit_h2_results,
  file.path(data_path, "table_ologit_h2_results_core.csv"),
  row.names = FALSE
)

write.csv(
  ologit_h2_results_or,
  file.path(data_path, "table_ologit_h2_results_odds_ratios.csv"),
  row.names = FALSE
)

saveRDS(
  ologit_h2_pooled,
  file.path(data_path, "ologit_h2_pooled.rds")
)

cat("Ordered-logit robustness outputs saved.\n")


# ------------------------------------------------------------
# 9.1.7 Export Table 8
# ------------------------------------------------------------

table8_ologit_pooled <- ologit_h2_results_or %>%
  dplyr::filter(term %in% c(
    "partySPD",
    "partyFDP",
    "partyGreens",
    "partyLeft",
    "partyAfD",
    "post2022",
    "partySPD:post2022",
    "partyFDP:post2022",
    "partyGreens:post2022",
    "partyLeft:post2022",
    "partyAfD:post2022",
    "pol_interest",
    "age",
    "genderFemale",
    "educationMedium",
    "educationHigh",
    "east_westEast"
  )) %>%
  dplyr::mutate(
    term_order = dplyr::case_when(
      term == "partySPD" ~ 1,
      term == "partyFDP" ~ 2,
      term == "partyGreens" ~ 3,
      term == "partyLeft" ~ 4,
      term == "partyAfD" ~ 5,
      term == "post2022" ~ 6,
      term == "partySPD:post2022" ~ 7,
      term == "partyFDP:post2022" ~ 8,
      term == "partyGreens:post2022" ~ 9,
      term == "partyLeft:post2022" ~ 10,
      term == "partyAfD:post2022" ~ 11,
      term == "pol_interest" ~ 12,
      term == "age" ~ 13,
      term == "genderFemale" ~ 14,
      term == "educationMedium" ~ 15,
      term == "educationHigh" ~ 16,
      term == "east_westEast" ~ 17
    ),
    term_label = dplyr::case_when(
      term == "partySPD" ~ "SPD",
      term == "partyFDP" ~ "FDP",
      term == "partyGreens" ~ "Greens",
      term == "partyLeft" ~ "Left",
      term == "partyAfD" ~ "AfD",
      term == "post2022" ~ "Post-2022",
      term == "partySPD:post2022" ~ "SPD x Post-2022",
      term == "partyFDP:post2022" ~ "FDP x Post-2022",
      term == "partyGreens:post2022" ~ "Greens x Post-2022",
      term == "partyLeft:post2022" ~ "Left x Post-2022",
      term == "partyAfD:post2022" ~ "AfD x Post-2022",
      term == "pol_interest" ~ "Political interest",
      term == "age" ~ "Age",
      term == "genderFemale" ~ "Female",
      term == "educationMedium" ~ "Education: Medium",
      term == "educationHigh" ~ "Education: High",
      term == "east_westEast" ~ "East Germany",
      TRUE ~ term
    )
  ) %>%
  dplyr::arrange(term_order) %>%
  dplyr::select(
    term_label,
    estimate_display,
    se_display,
    odds_ratio_display
  )

cat("Table 8 display object created.\n")
print(table8_ologit_pooled)

html_table8 <- c(
  "<html>",
  "<head><meta charset='UTF-8'></head>",
  "<body>",
  "<table style='border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; font-size: 12px;'>",
  "<caption style='caption-side: top; text-align: center; font-weight: bold; padding: 8px;'>",
  "Table 8. Ordered-Logit Robustness Check: Pooled Post-2022 Interaction Model with Demographic Controls",
  "</caption>",
  "<tr><td colspan='4' style='border-top: 1px solid black;'></td></tr>",
  "<tr>",
  "<th style='text-align: left; padding: 6px;'>Term</th>",
  "<th style='padding: 6px;'>Coefficient</th>",
  "<th style='padding: 6px;'>Std. Error</th>",
  "<th style='padding: 6px;'>Odds Ratio</th>",
  "</tr>",
  "<tr><td colspan='4' style='border-bottom: 1px solid black;'></td></tr>"
)

for (i in 1:nrow(table8_ologit_pooled)) {
  html_table8 <- c(
    html_table8,
    paste0(
      "<tr>",
      "<td style='text-align: left; padding: 6px;'>", table8_ologit_pooled$term_label[i], "</td>",
      "<td style='padding: 6px;'>", table8_ologit_pooled$estimate_display[i], "</td>",
      "<td style='padding: 6px;'>", table8_ologit_pooled$se_display[i], "</td>",
      "<td style='padding: 6px;'>", table8_ologit_pooled$odds_ratio_display[i], "</td>",
      "</tr>"
    )
  )
}

html_table8 <- c(
  html_table8,
  paste0(
    "<tr><td colspan='4' style='border-top: 1px solid black; text-align: left; padding: 8px;'>",
    "<strong>Note:</strong> Entries are ordered-logit coefficients, with model-based standard errors in parentheses and odds ratios in the final column. ",
    "The model includes controls for political interest, age, gender, education, and East/West residence. ",
    "CDU/CSU is the reference party category; male, low education, and West Germany are the demographic reference categories. ",
    "BSW is excluded because no pre-2022 observations are available. ",
    "Positive interaction coefficients indicate convergence toward the CDU/CSU baseline after 2022; negative interaction coefficients indicate divergence.",
    "</td></tr>"
  ),
  "</table>",
  "</body>",
  "</html>"
)

writeLines(
  html_table8,
  file.path(data_path, "table8_ologit_pooled_h2.html")
)

cat("Table 8 saved as table8_ologit_pooled_h2.html\n")



# ------------------------------------------------------------
# 9.2 Excluding W8
# ------------------------------------------------------------
# W8 is the only analysed wave that measures party preference
# using reported vote choice in the 2017 federal election rather
# than current vote intention.
#
# This robustness check re-estimates the pooled H2 models without
# W8 to assess whether the pooled interaction results depend on
# this measurement difference.
#
# The specification is the same as in the pooled model in Thesis
# Chapter 5.3.


# ------------------------------------------------------------
# 9.2.1 Build pooled dataset without W8
# ------------------------------------------------------------

gles_h2_no_w8 <- gles_h2 %>%
  dplyr::filter(wave != "W8")

cat("Reduced pooled dataset without W8 created.\n")
cat("Observations in original gles_h2:", nrow(gles_h2), "\n")
cat("Observations in gles_h2_no_w8:", nrow(gles_h2_no_w8), "\n")
cat("Remaining waves:\n")
print(table(gles_h2_no_w8$wave))


# ------------------------------------------------------------
# 9.2.2 Re-estimate pooled main-effects model without W8
# ------------------------------------------------------------

lm_h2_main_effects_no_w8 <- lm(
  defence ~ party + post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2_no_w8
)

ols_h2_main_effects_no_w8_clustered <- lm_robust(
  defence ~ party + post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2_no_w8,
  clusters = lfdn_int,
  se_type = "CR0"
)

se_h2_main_effects_no_w8_clustered <- robust_se(
  ols_h2_main_effects_no_w8_clustered
)


# ------------------------------------------------------------
# 9.2.3 Re-estimate pooled interaction model without W8
# ------------------------------------------------------------

pooled_interaction_lm_no_w8 <- lm(
  defence ~ party * post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2_no_w8
)

ols_h2_pooled_no_w8_clustered <- lm_robust(
  defence ~ party * post2022 + pol_interest + age + gender + education + east_west,
  data = gles_h2_no_w8,
  clusters = lfdn_int,
  se_type = "CR0"
)

se_h2_pooled_no_w8_clustered <- robust_se(
  ols_h2_pooled_no_w8_clustered
)


# ------------------------------------------------------------
# 9.2.4 Joint F-test of interaction terms without W8
# ------------------------------------------------------------

vcov_clustered_no_w8 <- vcovCL(
  pooled_interaction_lm_no_w8,
  cluster = ~ lfdn_int,
  data = gles_h2_no_w8
)

joint_test_h2_no_w8_clustered <- linearHypothesis(
  pooled_interaction_lm_no_w8,
  c(
    "partySPD:post2022 = 0",
    "partyFDP:post2022 = 0",
    "partyGreens:post2022 = 0",
    "partyLeft:post2022 = 0",
    "partyAfD:post2022 = 0"
  ),
  vcov. = vcov_clustered_no_w8,
  test = "F"
)

joint_f_value_no_w8 <- round(joint_test_h2_no_w8_clustered[2, "F"], 2)
joint_numdf_no_w8 <- joint_test_h2_no_w8_clustered[2, "Df"]
joint_dendf_no_w8 <- joint_test_h2_no_w8_clustered[2, "Res.Df"]
joint_p_value_no_w8 <- joint_test_h2_no_w8_clustered[2, "Pr(>F)"]

joint_p_display_no_w8 <- ifelse(
  joint_p_value_no_w8 < 0.001,
  "p < .001",
  paste0("p = ", sprintf("%.3f", joint_p_value_no_w8))
)

table9_joint_test_note <- paste0(
  "Joint F-test of interaction terms: F(",
  joint_numdf_no_w8, ", ",
  joint_dendf_no_w8, ") = ",
  sprintf("%.2f", joint_f_value_no_w8),
  ", ",
  joint_p_display_no_w8,
  "."
)

cat("Joint F-test without W8 completed.\n")
print(joint_test_h2_no_w8_clustered)


# ------------------------------------------------------------
# 9.2.5 Check model samples and coefficient equivalence
# ------------------------------------------------------------

cat("No-W8 pooled observation counts:\n")
cat("Main-effects model N:", nobs(lm_h2_main_effects_no_w8), "\n")
cat("Interaction model N:", nobs(pooled_interaction_lm_no_w8), "\n")
cat("Clustered interaction model N:", nobs(ols_h2_pooled_no_w8_clustered), "\n")

cat("Coefficient equivalence check for interaction model:\n")
print(all.equal(
  unname(coef(pooled_interaction_lm_no_w8)),
  unname(coef(ols_h2_pooled_no_w8_clustered)),
  tolerance = 1e-10
))


# ------------------------------------------------------------
# 9.2.6 Export Table 9
# ------------------------------------------------------------

stargazer(
  lm_h2_main_effects_no_w8,
  pooled_interaction_lm_no_w8,
  se = list(se_h2_main_effects_no_w8_clustered, se_h2_pooled_no_w8_clustered),
  type = "text",
  title = "Table 9. Robustness Check: Pooled Model Excluding Wave 8, with Demographic Controls",
  column.labels = c("Main effects + controls", "With interactions + controls"),
  dep.var.labels = "Defence spending support (1-5)",
  covariate.labels = c(
    "SPD",
    "FDP",
    "Greens",
    "Left",
    "AfD",
    "Post-2022",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany",
    "SPD x Post-2022",
    "FDP x Post-2022",
    "Greens x Post-2022",
    "Left x Post-2022",
    "AfD x Post-2022"
  ),
  omit = "Constant",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Clustered standard errors (CR0) by respondent ID in parentheses.",
    "Models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the party reference category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Wave 8 is excluded because it uses reported actual vote rather than vote intention to measure party preference.",
    "BSW is excluded because no pre-2022 observations are available.",
    "Positive interaction = convergence toward CDU/CSU post-2022.",
    "Negative interaction = divergence from CDU/CSU post-2022.",
    "* p < .10, ** p < .05, *** p < .01.",
    table9_joint_test_note
  ),
  notes.append = FALSE
)

stargazer(
  lm_h2_main_effects_no_w8,
  pooled_interaction_lm_no_w8,
  se = list(se_h2_main_effects_no_w8_clustered, se_h2_pooled_no_w8_clustered),
  type = "html",
  out = file.path(data_path, "table9_pooled_no_w8.html"),
  title = "Table 9. Robustness Check: Pooled Model Excluding Wave 8, with Demographic Controls",
  column.labels = c("Main effects + controls", "With interactions + controls"),
  dep.var.labels = "Defence spending support (1-5)",
  covariate.labels = c(
    "SPD",
    "FDP",
    "Greens",
    "Left",
    "AfD",
    "Post-2022",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany",
    "SPD x Post-2022",
    "FDP x Post-2022",
    "Greens x Post-2022",
    "Left x Post-2022",
    "AfD x Post-2022"
  ),
  omit = "Constant",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Clustered standard errors (CR0) by respondent ID in parentheses.",
    "Models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the party reference category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Wave 8 is excluded because it uses reported actual vote rather than vote intention to measure party preference.",
    "BSW is excluded because no pre-2022 observations are available.",
    "Positive interaction = convergence toward CDU/CSU post-2022.",
    "Negative interaction = divergence from CDU/CSU post-2022.",
    "* p < .10, ** p < .05, *** p < .01.",
    table9_joint_test_note
  ),
  notes.append = FALSE
)

cat("Excluding-W8 robustness check completed.\n")
cat("Table 9 saved as table9_pooled_no_w8.html\n")
print(summary(ols_h2_pooled_no_w8_clustered))


# ------------------------------------------------------------
# 9.3 Selective weighting check
# ------------------------------------------------------------
# This robustness check re-estimates selected wave-specific OLS
# models with survey weights for the waves in which suitable
# weight variables are available.
#
# Weighted and unweighted estimates are compared on the same
# weighted-eligible complete-case samples. This keeps the comparison
# focused on the effect of weighting rather than on changes in the
# sample composition.
#
# The specification mirrors the controlled wave-by-wave models from
# Thesis Chapter 5.2.


# ------------------------------------------------------------
# 9.3.1 Build weighted complete-case samples
# ------------------------------------------------------------

weighted_wave_data <- function(w) {
  gles_analysis %>%
    dplyr::filter(
      wave == w,
      !is.na(defence),
      !is.na(party),
      !is.na(pol_interest),
      !is.na(age),
      !is.na(gender),
      !is.na(education),
      !is.na(east_west),
      !is.na(weight),
      weight > 0
    ) %>%
    droplevels() %>%
    as.data.frame()
}

weighted_waves <- c("W10", "W11", "W17", "W26", "W32")

cat("Weighted robustness samples created.\n")
cat("N W10:", nrow(weighted_wave_data("W10")), "\n")
cat("N W11:", nrow(weighted_wave_data("W11")), "\n")
cat("N W17:", nrow(weighted_wave_data("W17")), "\n")
cat("N W26:", nrow(weighted_wave_data("W26")), "\n")
cat("N W32:", nrow(weighted_wave_data("W32")), "\n")


# ------------------------------------------------------------
# 9.3.2 Estimate weighted controlled models
# ------------------------------------------------------------

model_w10_weighted_robustness <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W10"),
  weights = weight,
  se_type = "HC2"
)

model_w11_weighted_robustness <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W11"),
  weights = weight,
  se_type = "HC2"
)

model_w17_weighted_robustness <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W17"),
  weights = weight,
  se_type = "HC2"
)

model_w26_weighted_robustness <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W26"),
  weights = weight,
  se_type = "HC2"
)

model_w32_weighted_robustness <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W32"),
  weights = weight,
  se_type = "HC2"
)

cat("Weighted controlled models estimated.\n")


# ------------------------------------------------------------
# 9.3.3 Estimate unweighted models on the same samples
# ------------------------------------------------------------
# These models use the same observations as the weighted models.
# This allows a direct weighted-versus-unweighted comparison.

model_w10_unweighted_same_sample <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W10"),
  se_type = "HC2"
)

model_w11_unweighted_same_sample <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W11"),
  se_type = "HC2"
)

model_w17_unweighted_same_sample <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W17"),
  se_type = "HC2"
)

model_w26_unweighted_same_sample <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W26"),
  se_type = "HC2"
)

model_w32_unweighted_same_sample <- estimatr::lm_robust(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W32"),
  se_type = "HC2"
)

cat("Unweighted same-sample comparison models estimated.\n")


# ------------------------------------------------------------
# 9.3.4 Compare weighted and unweighted coefficients
# ------------------------------------------------------------

extract_core_terms_weight_check <- function(model, wave_label) {
  tibble::tibble(
    wave = wave_label,
    spd = unname(coef(model)["partySPD"]),
    greens = unname(coef(model)["partyGreens"]),
    afd = unname(coef(model)["partyAfD"]),
    pol_interest = unname(coef(model)["pol_interest"]),
    age = unname(coef(model)["age"]),
    female = unname(coef(model)["genderFemale"]),
    east = unname(coef(model)["east_westEast"])
  )
}

weighted_core <- dplyr::bind_rows(
  extract_core_terms_weight_check(model_w10_weighted_robustness, "W10"),
  extract_core_terms_weight_check(model_w11_weighted_robustness, "W11"),
  extract_core_terms_weight_check(model_w17_weighted_robustness, "W17"),
  extract_core_terms_weight_check(model_w26_weighted_robustness, "W26"),
  extract_core_terms_weight_check(model_w32_weighted_robustness, "W32")
) %>%
  dplyr::rename(
    spd_weighted = spd,
    greens_weighted = greens,
    afd_weighted = afd,
    pol_interest_weighted = pol_interest,
    age_weighted = age,
    female_weighted = female,
    east_weighted = east
  )

unweighted_core <- dplyr::bind_rows(
  extract_core_terms_weight_check(model_w10_unweighted_same_sample, "W10"),
  extract_core_terms_weight_check(model_w11_unweighted_same_sample, "W11"),
  extract_core_terms_weight_check(model_w17_unweighted_same_sample, "W17"),
  extract_core_terms_weight_check(model_w26_unweighted_same_sample, "W26"),
  extract_core_terms_weight_check(model_w32_unweighted_same_sample, "W32")
) %>%
  dplyr::rename(
    spd_unweighted = spd,
    greens_unweighted = greens,
    afd_unweighted = afd,
    pol_interest_unweighted = pol_interest,
    age_unweighted = age,
    female_unweighted = female,
    east_unweighted = east
  )

table_weighted_robustness_comparison <- unweighted_core %>%
  dplyr::left_join(weighted_core, by = "wave") %>%
  dplyr::mutate(
    spd_difference = spd_weighted - spd_unweighted,
    greens_difference = greens_weighted - greens_unweighted,
    afd_difference = afd_weighted - afd_unweighted,
    pol_interest_difference = pol_interest_weighted - pol_interest_unweighted,
    age_difference = age_weighted - age_unweighted,
    female_difference = female_weighted - female_unweighted,
    east_difference = east_weighted - east_unweighted
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = -wave,
      .fns = ~ round(.x, 3)
    )
  )

cat("Weighted/unweighted same-sample comparison table created.\n")
print(table_weighted_robustness_comparison)

write.csv(
  table_weighted_robustness_comparison,
  file.path(data_path, "table_weighted_robustness_comparison.csv"),
  row.names = FALSE
)

cat("Comparison table saved as table_weighted_robustness_comparison.csv\n")


# ------------------------------------------------------------
# 9.3.5 Rebuild plain weighted lm models for table export
# ------------------------------------------------------------
# The lm objects are used for stargazer table export. The same
# weighted samples and formulas are used as in the weighted
# lm_robust models above.
#
# The older object names are kept to avoid breaking later code.

lm_table_5_5_3_w10 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W10"),
  weights = weight
)

lm_table_5_5_3_w11 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W11"),
  weights = weight
)

lm_table_5_5_3_w17 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W17"),
  weights = weight
)

lm_table_5_5_3_w26 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W26"),
  weights = weight
)

lm_table_5_5_3_w32 <- lm(
  defence ~ party + pol_interest + age + gender + education + east_west,
  data = weighted_wave_data("W32"),
  weights = weight
)

cat("Plain weighted lm models for table export created.\n")
cat("N W10:", nobs(lm_table_5_5_3_w10), "\n")
cat("N W11:", nobs(lm_table_5_5_3_w11), "\n")
cat("N W17:", nobs(lm_table_5_5_3_w17), "\n")
cat("N W26:", nobs(lm_table_5_5_3_w26), "\n")
cat("N W32:", nobs(lm_table_5_5_3_w32), "\n")


# ------------------------------------------------------------
# 9.3.6 Compute HC2 robust standard errors
# ------------------------------------------------------------

se_table_5_5_3_w10 <- sqrt(diag(sandwich::vcovHC(lm_table_5_5_3_w10, type = "HC2")))
se_table_5_5_3_w11 <- sqrt(diag(sandwich::vcovHC(lm_table_5_5_3_w11, type = "HC2")))
se_table_5_5_3_w17 <- sqrt(diag(sandwich::vcovHC(lm_table_5_5_3_w17, type = "HC2")))
se_table_5_5_3_w26 <- sqrt(diag(sandwich::vcovHC(lm_table_5_5_3_w26, type = "HC2")))
se_table_5_5_3_w32 <- sqrt(diag(sandwich::vcovHC(lm_table_5_5_3_w32, type = "HC2")))

cat("HC2 robust standard errors for weighted tables created.\n")


# ------------------------------------------------------------
# 9.3.7 Check weighted table samples and coefficients
# ------------------------------------------------------------

cat("Weighted export sample sizes:\n")
cat("W10:", nrow(weighted_wave_data("W10")), " | model N =", nobs(lm_table_5_5_3_w10), "\n")
cat("W11:", nrow(weighted_wave_data("W11")), " | model N =", nobs(lm_table_5_5_3_w11), "\n")
cat("W17:", nrow(weighted_wave_data("W17")), " | model N =", nobs(lm_table_5_5_3_w17), "\n")
cat("W26:", nrow(weighted_wave_data("W26")), " | model N =", nobs(lm_table_5_5_3_w26), "\n")
cat("W32:", nrow(weighted_wave_data("W32")), " | model N =", nobs(lm_table_5_5_3_w32), "\n")

stopifnot(
  nrow(weighted_wave_data("W10")) == nobs(lm_table_5_5_3_w10),
  nrow(weighted_wave_data("W11")) == nobs(lm_table_5_5_3_w11),
  nrow(weighted_wave_data("W17")) == nobs(lm_table_5_5_3_w17),
  nrow(weighted_wave_data("W26")) == nobs(lm_table_5_5_3_w26),
  nrow(weighted_wave_data("W32")) == nobs(lm_table_5_5_3_w32)
)

cat("Coefficient equality checks: weighted lm vs weighted lm_robust\n")
print(all.equal(unname(coef(lm_table_5_5_3_w10)), unname(coef(model_w10_weighted_robustness)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_table_5_5_3_w11)), unname(coef(model_w11_weighted_robustness)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_table_5_5_3_w17)), unname(coef(model_w17_weighted_robustness)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_table_5_5_3_w26)), unname(coef(model_w26_weighted_robustness)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_table_5_5_3_w32)), unname(coef(model_w32_weighted_robustness)), tolerance = 1e-10))


# ------------------------------------------------------------
# 9.3.8 Preview Table 10: Weighted pre-2022 models
# ------------------------------------------------------------

stargazer(
  lm_table_5_5_3_w10,
  lm_table_5_5_3_w11,
  lm_table_5_5_3_w17,
  se = list(
    se_table_5_5_3_w10,
    se_table_5_5_3_w11,
    se_table_5_5_3_w17
  ),
  type = "text",
  title = "Table 10. Selective Weighted OLS Robustness Check (Pre-2022 Waves) with Demographic Controls",
  column.labels = c("W10", "W11", "W17"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Robust standard errors (HC2) in parentheses.",
    "All models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the reference party category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Models are estimated with survey weights.",
    "Only waves with suitable weight variables are included.",
    "* p < .10, ** p < .05, *** p < .01."
  ),
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.3.9 Preview Table 11: Weighted post-2022 models
# ------------------------------------------------------------

stargazer(
  lm_table_5_5_3_w26,
  lm_table_5_5_3_w32,
  se = list(
    se_table_5_5_3_w26,
    se_table_5_5_3_w32
  ),
  type = "text",
  title = "Table 11. Selective Weighted OLS Robustness Check (Post-2022 Waves) with Demographic Controls",
  column.labels = c("W26", "W32"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Robust standard errors (HC2) in parentheses.",
    "All models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the reference party category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Models are estimated with survey weights.",
    "Only waves with suitable weight variables are included.",
    "BSW was not available prior to 2024.",
    "* p < .10, ** p < .05, *** p < .01."
  ),
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.3.10 Save Table 10: Weighted pre-2022 models
# ------------------------------------------------------------

stargazer(
  lm_table_5_5_3_w10,
  lm_table_5_5_3_w11,
  lm_table_5_5_3_w17,
  se = list(
    se_table_5_5_3_w10,
    se_table_5_5_3_w11,
    se_table_5_5_3_w17
  ),
  type = "html",
  out = file.path(data_path, "table10_weighted_robustness_pre.html"),
  title = "Table 10. Selective Weighted OLS Robustness Check (Pre-2022 Waves) with Demographic Controls",
  column.labels = c("W10", "W11", "W17"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Robust standard errors (HC2) in parentheses.",
    "All models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the reference party category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Models are estimated with survey weights.",
    "Only waves with suitable weight variables are included.",
    "* p < .10, ** p < .05, *** p < .01."
  ),
  notes.append = FALSE
)

cat("Table 10 saved as table10_weighted_robustness_pre.html\n")


# ------------------------------------------------------------
# 9.3.11 Save Table 11: Weighted post-2022 models
# ------------------------------------------------------------

stargazer(
  lm_table_5_5_3_w26,
  lm_table_5_5_3_w32,
  se = list(
    se_table_5_5_3_w26,
    se_table_5_5_3_w32
  ),
  type = "html",
  out = file.path(data_path, "table11_weighted_robustness_post.html"),
  title = "Table 11. Selective Weighted OLS Robustness Check (Post-2022 Waves) with Demographic Controls",
  column.labels = c("W26", "W32"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = c(
    "Robust standard errors (HC2) in parentheses.",
    "All models include age, gender, education, and East/West residence controls.",
    "CDU/CSU is the reference party category.",
    "Male, low education, and West Germany are the demographic reference categories.",
    "Models are estimated with survey weights.",
    "Only waves with suitable weight variables are included.",
    "BSW was not available prior to 2024.",
    "* p < .10, ** p < .05, *** p < .01."
  ),
  notes.append = FALSE
)

cat("Table 11 saved as table11_weighted_robustness_post.html\n")
cat("Selective weighting robustness check completed.\n")



# ------------------------------------------------------------
# 9.4 Party identification instead of vote intention
# ------------------------------------------------------------
# This robustness check re-estimates the wave-specific OLS models
# using party identification instead of vote intention as the
# partisan predictor.
#
# The specification mirrors the controlled wave-by-wave models from
# Thesis Chapter 5.2. The purpose is to check whether the main
# partisan patterns also appear when partisanship is measured as
# party identification rather than electoral preference.


# ------------------------------------------------------------
# 9.4.1 Build party-identification variable across waves
# ------------------------------------------------------------

party_id_all <- bind_rows(
  w1to9 %>% transmute(
    lfdn = lfdn,
    wave = "W8",
    party_id = party_name(clean_party(kp8_2090a))
  ),
  w10 %>% transmute(
    lfdn = lfdn,
    wave = "W10",
    party_id = party_name(clean_party(kp10_2090a))
  ),
  w11 %>% transmute(
    lfdn = lfdn,
    wave = "W11",
    party_id = party_name(clean_party(kp11_2090a))
  ),
  w17 %>% transmute(
    lfdn = lfdn,
    wave = "W17",
    party_id = party_name(clean_party(kp17_2090a))
  ),
  w22 %>% transmute(
    lfdn = lfdn,
    wave = "W22",
    party_id = party_name(clean_party(kp22_2090a))
  ),
  w24 %>% transmute(
    lfdn = lfdn,
    wave = "W24",
    party_id = party_name(clean_party(kp24_2090a))
  ),
  w26 %>% transmute(
    lfdn = lfdn,
    wave = "W26",
    party_id = party_name(clean_party(kp26_2090a))
  ),
  w32 %>% transmute(
    lfdn = lfdn,
    wave = "W32",
    party_id = party_name(clean_party(kp32_2090a))
  )
)


# ------------------------------------------------------------
# 9.4.2 Build controlled party-ID robustness dataset
# ------------------------------------------------------------

gles_partyid_robustness <- gles_analysis %>%
  dplyr::select(
    lfdn, wave, defence, pol_interest,
    age, gender, education, east_west
  ) %>%
  dplyr::left_join(party_id_all, by = c("lfdn", "wave")) %>%
  dplyr::mutate(
    party_id = factor(
      party_id,
      levels = c("CDU/CSU", "SPD", "Greens", "FDP", "Left", "AfD", "BSW")
    )
  ) %>%
  dplyr::filter(
    !is.na(defence),
    !is.na(party_id),
    !is.na(pol_interest),
    !is.na(age),
    !is.na(gender),
    !is.na(education),
    !is.na(east_west)
  )

cat("Party-identification robustness dataset with demographic controls created.\n")


# ------------------------------------------------------------
# 9.4.3 Compare party-ID sample sizes to main model samples
# ------------------------------------------------------------

partyid_n_check <- gles_analysis %>%
  dplyr::filter(
    !is.na(defence),
    !is.na(party),
    !is.na(pol_interest),
    !is.na(age),
    !is.na(gender),
    !is.na(education),
    !is.na(east_west)
  ) %>%
  dplyr::count(wave, name = "n_main_vote_intention_controls") %>%
  dplyr::left_join(
    gles_partyid_robustness %>%
      dplyr::count(wave, name = "n_party_identification_controls"),
    by = "wave"
  ) %>%
  dplyr::mutate(
    n_party_identification_controls =
      dplyr::coalesce(n_party_identification_controls, 0L),
    n_lost = n_main_vote_intention_controls - n_party_identification_controls
  )

print(partyid_n_check)

cat("Party-identification distribution by wave:\n")
print(gles_partyid_robustness %>% dplyr::count(wave, party_id), n = 100)


# ------------------------------------------------------------
# 9.4.4 Estimate controlled party-ID models
# ------------------------------------------------------------

partyid_wave_data <- function(w) {
  gles_partyid_robustness %>%
    dplyr::filter(wave == w)
}

partyid_formula <- defence ~ party_id + pol_interest + age + gender + education + east_west

ols_partyid_w8 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W8"),
  se_type = "HC2"
)

ols_partyid_w10 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W10"),
  se_type = "HC2"
)

ols_partyid_w11 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W11"),
  se_type = "HC2"
)

ols_partyid_w17 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W17"),
  se_type = "HC2"
)

ols_partyid_w22 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W22"),
  se_type = "HC2"
)

ols_partyid_w24 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W24"),
  se_type = "HC2"
)

ols_partyid_w26 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W26"),
  se_type = "HC2"
)

ols_partyid_w32 <- estimatr::lm_robust(
  partyid_formula,
  data = partyid_wave_data("W32"),
  se_type = "HC2"
)

cat("All eight controlled party-identification robustness models estimated.\n")

cat("Model Ns:\n")
cat("W8 :", nobs(ols_partyid_w8), "\n")
cat("W10:", nobs(ols_partyid_w10), "\n")
cat("W11:", nobs(ols_partyid_w11), "\n")
cat("W17:", nobs(ols_partyid_w17), "\n")
cat("W22:", nobs(ols_partyid_w22), "\n")
cat("W24:", nobs(ols_partyid_w24), "\n")
cat("W26:", nobs(ols_partyid_w26), "\n")
cat("W32:", nobs(ols_partyid_w32), "\n")


# ------------------------------------------------------------
# 9.4.5 Rebuild plain lm models for table export
# ------------------------------------------------------------
# stargazer works more smoothly with regular lm objects.
# These lm models use the same formulas and samples as the
# lm_robust models above. The robust HC2 standard errors are then
# supplied manually.

lm_partyid_w8 <- lm(
  partyid_formula,
  data = partyid_wave_data("W8")
)

lm_partyid_w10 <- lm(
  partyid_formula,
  data = partyid_wave_data("W10")
)

lm_partyid_w11 <- lm(
  partyid_formula,
  data = partyid_wave_data("W11")
)

lm_partyid_w17 <- lm(
  partyid_formula,
  data = partyid_wave_data("W17")
)

lm_partyid_w22 <- lm(
  partyid_formula,
  data = partyid_wave_data("W22")
)

lm_partyid_w24 <- lm(
  partyid_formula,
  data = partyid_wave_data("W24")
)

lm_partyid_w26 <- lm(
  partyid_formula,
  data = partyid_wave_data("W26")
)

lm_partyid_w32 <- lm(
  partyid_formula,
  data = partyid_wave_data("W32")
)

se_partyid_w8  <- robust_se(ols_partyid_w8)
se_partyid_w10 <- robust_se(ols_partyid_w10)
se_partyid_w11 <- robust_se(ols_partyid_w11)
se_partyid_w17 <- robust_se(ols_partyid_w17)
se_partyid_w22 <- robust_se(ols_partyid_w22)
se_partyid_w24 <- robust_se(ols_partyid_w24)
se_partyid_w26 <- robust_se(ols_partyid_w26)
se_partyid_w32 <- robust_se(ols_partyid_w32)

cat("Plain lm objects for party-ID table export created.\n")


# ------------------------------------------------------------
# 9.4.6 Preview Table 12: Pre-2022 party-ID models
# ------------------------------------------------------------

stargazer(
  lm_partyid_w8, lm_partyid_w10, lm_partyid_w11, lm_partyid_w17,
  se = list(se_partyid_w8, se_partyid_w10, se_partyid_w11, se_partyid_w17),
  type = "text",
  title = "Table 12. Party-Identification Robustness Check (Pre-2022 Waves) with Demographic Controls",
  column.labels = c("W8 Sep 2017", "W10 Nov 2018", "W11 Jun 2019", "W17 Jul 2021"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the party reference category. Male, low education, and West Germany are the demographic reference categories. Party identification is used instead of vote intention. Political interest is reverse-coded so that higher values indicate greater political interest. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.4.7 Save Table 12: Pre-2022 party-ID models
# ------------------------------------------------------------

stargazer(
  lm_partyid_w8, lm_partyid_w10, lm_partyid_w11, lm_partyid_w17,
  se = list(se_partyid_w8, se_partyid_w10, se_partyid_w11, se_partyid_w17),
  type = "html",
  out = file.path(data_path, "table12_partyid_robustness_pre.html"),
  title = "Table 12. Party-Identification Robustness Check (Pre-2022 Waves) with Demographic Controls",
  column.labels = c("W8 Sep 2017", "W10 Nov 2018", "W11 Jun 2019", "W17 Jul 2021"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the party reference category. Male, low education, and West Germany are the demographic reference categories. Party identification is used instead of vote intention. Political interest is reverse-coded so that higher values indicate greater political interest. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.4.8 Preview Table 13: Post-2022 party-ID models
# ------------------------------------------------------------

stargazer(
  lm_partyid_w22, lm_partyid_w24, lm_partyid_w26, lm_partyid_w32,
  se = list(se_partyid_w22, se_partyid_w24, se_partyid_w26, se_partyid_w32),
  type = "text",
  title = "Table 13. Party-Identification Robustness Check (Post-2022 Waves) with Demographic Controls",
  column.labels = c("W22 May 2022", "W24 May 2023", "W26 Jun 2024", "W32 May 2025"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the party reference category. Male, low education, and West Germany are the demographic reference categories. Party identification is used instead of vote intention. Political interest is reverse-coded so that higher values indicate greater political interest. BSW was not available prior to 2024. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.4.9 Save Table 13: Post-2022 party-ID models
# ------------------------------------------------------------

stargazer(
  lm_partyid_w22, lm_partyid_w24, lm_partyid_w26, lm_partyid_w32,
  se = list(se_partyid_w22, se_partyid_w24, se_partyid_w26, se_partyid_w32),
  type = "html",
  out = file.path(data_path, "table13_partyid_robustness_post.html"),
  title = "Table 13. Party-Identification Robustness Check (Post-2022 Waves) with Demographic Controls",
  column.labels = c("W22 May 2022", "W24 May 2023", "W26 Jun 2024", "W32 May 2025"),
  covariate.labels = c(
    "SPD",
    "Greens",
    "FDP",
    "Left",
    "AfD",
    "BSW",
    "Political interest",
    "Age",
    "Female",
    "Education: Medium",
    "Education: High",
    "East Germany"
  ),
  dep.var.labels = "Defence spending support (1-5)",
  omit.stat = c("f", "ser"),
  star.cutoffs = c(0.1, 0.05, 0.01),
  notes = "Robust standard errors (HC2) in parentheses. All models include age, gender, education, and East/West residence controls. CDU/CSU is the party reference category. Male, low education, and West Germany are the demographic reference categories. Party identification is used instead of vote intention. Political interest is reverse-coded so that higher values indicate greater political interest. BSW was not available prior to 2024. * p < .10, ** p < .05, *** p < .01.",
  notes.append = FALSE
)


# ------------------------------------------------------------
# 9.4.10 Check coefficient equivalence
# ------------------------------------------------------------
# These checks confirm that the lm objects used for table export
# match the lm_robust objects used for the actual robustness models.

cat("Coefficient equality checks: party-ID lm vs lm_robust\n")
print(all.equal(unname(coef(lm_partyid_w8)),  unname(coef(ols_partyid_w8)),  tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w10)), unname(coef(ols_partyid_w10)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w11)), unname(coef(ols_partyid_w11)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w17)), unname(coef(ols_partyid_w17)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w22)), unname(coef(ols_partyid_w22)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w24)), unname(coef(ols_partyid_w24)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w26)), unname(coef(ols_partyid_w26)), tolerance = 1e-10))
print(all.equal(unname(coef(lm_partyid_w32)), unname(coef(ols_partyid_w32)), tolerance = 1e-10))

cat("Party-identification robustness check completed.\n")
cat("Table 12 saved as table12_partyid_robustness_pre.html\n")
cat("Table 13 saved as table13_partyid_robustness_post.html\n")



cat("Supplementary Russia-validation check completed.\n")

# ============================================================
# END OF SCRIPT
# ============================================================

cat("\nAnalysis script completed.\n")