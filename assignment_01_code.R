##***************************
## Assignment 1 Code
##
## Group 8
##
## 2026-10-27
##
##***************************
#Contributions:
## - Original author: Kristian (GitHub: kgonzales86). Wrote the originalscript.
## - Secondary contributor: Shagufta (GitHub: Shagufta-Shaikh).
##   Edit 1: Added PLE and PLE/Naiad to the retrotransposon lists, merged the
##           duplicated non-TE list, and added class counts to check the filters.
##   Edit 2: Replaced repeated mutate() calls with across().
##   Edit 3: Replaced five separate plots with one faceted histogram and
##           frequency polygon figure.
##   Minor: Updated file paths and added explanatory comments (marked "SS").
##***************************
rm(list = ls())

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

## Reading in data -----
#(SS) - fixed file paths
df_trich_sp = read_csv("trichinella_spiralis_te.csv")
df_trich_sp

df_trich_pap = read_csv("GCA_001447755.1.csv")
df_trich_pap

## Filtering the data -----
# EDIT 1 (SS): (a) The non-TE class list was defined twice (once per species) with identical contents; there is now one list for both species. (b) Added "PLE" and "PLE/Naiad" to both retrotransposon lists. They were previously dropped, excluding 460 spiralis and 406 papuae entries from the analysis. Counts of matching_class were added to check which classes of te and non-te are present.

df_trich_sp |>  count(matching_class, sort = TRUE)
df_trich_pap |> count(matching_class, sort = TRUE) 

ls_non_te = c("Simple_repeat", "Low_complexity", "Unknown",
              "simple_repeat", "Satellite", "rRNA", "snRNA")

df_actual_te_sp = df_trich_sp |> 
  filter(matching_class %notin% ls_non_te)

df_actual_te_sp #(SS)- to take a quick look at the df

#Filtering out all non-TE's from the Trichinella spiralis data frame

#(SS) added "PLE", "PLE/Naiad"
ls_retransposons_sp = c("LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", 
                        "SINE/tRNA-Deu-L2", "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro","LTR/DIRS","PLE", "PLE/Naiad") 

df_trich_sp_retro = df_actual_te_sp |> 
  filter(matching_class %in% ls_retransposons_sp)

df_trich_sp_retro

#Filtering so that only retrotransposons are present in the Trichinella spiralis
#data frame

df_actual_te_pap = df_trich_pap |> 
  filter(matching_class %notin% ls_non_te)

df_actual_te_pap #SS - to take a quick look at the df

#Filtering out all non-TE's from the Trichinella papuae data frame

#(SS) added "PLE", "PLE/Naiad"
ls_retransposons_pap = c("LINE", "LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro", "PLE", "PLE/Naiad")

df_trich_pap_retro = df_actual_te_pap |> 
  filter(matching_class %in% ls_retransposons_pap)

df_trich_pap_retro

#Filtering so that only retrotransposons are present in the Trichinella papuae
#data frame

## Combining the data frames -----

# EDIT 2 (SS): The cleaning step repeated the same str_remove_all() call in three separate mutate() lines for each species (six lines in total). These are replaced with one across() per species that applies the same cleaning to every column in ls_clean_cols. 

ls_clean_cols = c("qry_left", "no_bp_in_complement", "in_repeat_end")

df_trich_sp_final = df_trich_sp_retro |> 
  mutate(species = "spiralis",
         across(all_of(ls_clean_cols), ~ str_remove_all(.x, "-")))

df_trich_pap_final = df_trich_pap_retro |> 
  mutate(species = "papuae",
         across(all_of(ls_clean_cols), ~ str_remove_all(.x, "[()]")))

#Adding a column indicating which species before combining
#Formatting so that the bind_rows function can run without error.
#Formatting done was removing "()" and "-" from values in those 3 columns

df_trich_combined = bind_rows(df_trich_sp_final, df_trich_pap_final)
view(df_trich_combined)
#Combined data frames for more analysis options

## Figures -----

# EDIT 3 (SS): The original script made five separate plots (a histogram and a frequency polygon for each species, plus a combined frequency polygon), so the species had to be compared across different figures. They are replaced by one faceted figure that shows the histogram and the frequency polygon together, with the species stacked on a shared x-axis and the same y-axis scale, so they can be compared directly.

per_insert_fig <- df_trich_combined |> 
  ggplot(aes(x = perc_insert, fill = species, colour = species)) +
  geom_histogram(bins = 30, alpha = 0.4, show.legend = FALSE) +
  geom_freqpoly(bins = 30, linewidth = 0.8, show.legend = FALSE) +
  facet_wrap(~ species, ncol = 1) +
  labs(x = "Percent Insertion", y = "Frequency",
       title = "Percent Insertion of Retrotransposons in Two Trichinella Species")
per_insert_fig
ggsave("per_insert_fig.png")

## Potential novel code -----
  
summary(df_trich_sp_retro$perc_insert)
summary(df_trich_pap_retro$perc_insert)

df_trich_combined_lm = 
  lm(df_trich_combined$perc_insert~df_trich_combined$species)
summary(df_trich_combined_lm)

t.test(df_trich_pap_retro$perc_insert, df_trich_sp_retro$perc_insert)
