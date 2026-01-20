library(egg)
library(ggbeeswarm)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(ggsidekick); theme_set(theme_sleek())

root_dir <- here::here()

# plot
n_sim <- 300

res_df <- read_csv(file.path(root_dir, paste0("2026-01-05/simulation_selfcross_out_", n_sim, "_reps.csv")))
#res_df <- read_csv(file.path(root_dir, paste0("2025-11-09/simulation_selfcross_out_", n_sim, "_reps.csv")))
#res_df <- read_csv(file.path(root_dir, paste0("2025-11-09/simulation_selfcross_out_", n_sim, "_reps.csv")))

# specify true model for each species
correct_models <- tibble(
  species = c("capelin", "pacific cod", "pacific halibut"),
  correct_type = c("ST", "null", "T")
)

res_df |> 
  distinct(AIC, type, iter, species) |> 
  mutate(type = case_when(
    type == "1-000" ~ "null",
    type == "1-100" ~ "S",
    type == "1-010" ~ "T",
    type == "1-110" ~ "ST",
    TRUE ~ type
  )) |>
  mutate(
    type = factor(type, levels = c("null", "S", "T", "ST"))
  ) |> 
  mutate(delta_AIC = AIC - min(AIC), .by = c(iter, species)) |> 
  summarise(fraction_best = mean(delta_AIC == 0), .by = c(type, species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type) |>
  ggplot(aes(type, fraction_best, fill = true_model)) +
  facet_wrap(~species, ncol = 3) +
  guides(fill = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.1, 0.8),
  panel.spacing = unit(1, "lines"),
        plot.margin = margin(5, 10, 5, 5)) +
  coord_cartesian(expand = 0, ylim = c(0,1)) +
  geom_col(width = 0.9) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "model",
    y = "fraction of iterations with best AIC"
  ) +
  scale_fill_brewer(palette = "Set2", name = "operating model",
                    direction = -1)

ggsave(paste0(root_dir, "/figs/self_test_aic.png"), width = 20, height = 8, unit = "cm")


# Now wrangle the data a lot for plotting...
pdat <- res_df |> 
  # keep logit rhoE for now...
  filter(par_name %in% c("log_kappaS", "kappaT", "logit_rhoE", "RMSD")) |> 
  filter(species %in% c("capelin", "pacific halibut")) |> 
  drop_na(par_true) |>
  dplyr::select(-seed, -AIC) |> 
  # make wide so I can calculate AR1 and RMSD
  pivot_wider(names_from = "par_name", values_from = c("par_true", "par_hat")) |> 
  mutate(`par_true_AR(1)` = par_true_kappaT / (par_true_kappaT + 1),
         `par_hat_AR(1)` = par_hat_kappaT / (par_hat_kappaT + 1),
         # par_true_RMSD = ifelse(type == "1-100", 
         #                        sqrt(4 * exp(-2 * par_true_log_kappaS) * (1 - 0 / (1 + 0))),
         #                        sqrt(4 * exp(-2 * par_true_log_kappaS) * (1 - par_true_kappaT / (1 + par_true_kappaT)))),
         par_true_RMSD = sqrt(4 * exp(-2 * par_true_log_kappaS) * (1 - par_true_kappaT / (1 + par_true_kappaT))),
         par_hat_RMSD = ifelse(type == "1-100", 
                               sqrt(4 * exp(-2 * par_hat_log_kappaS) * (1 - 0 / (1 + 0))),
                               sqrt(4 * exp(-2 * par_hat_log_kappaS) * (1 - par_hat_kappaT / (1 + par_hat_kappaT)))),
         ) |> 
  # make long again for plotting easily
  pivot_longer(c("par_true_logit_rhoE", "par_true_log_kappaS", "par_true_kappaT",
                 "par_hat_logit_rhoE", "par_hat_log_kappaS", "par_hat_kappaT",
                 "par_true_AR(1)", "par_hat_AR(1)", "par_true_RMSD", "par_hat_RMSD"),
               names_to = "par_name") |> 
  drop_na(value) |> 
  # clean up the par names and split true vs estimated
  mutate(
    type_of_par = if_else(str_detect(par_name, "^par_true"), "par_true", "par_hat"),
    par_name = str_replace(par_name, "^par_(true|hat)_", "")
  ) |> 
  # rename the model types
  mutate(type = case_when(
    type == "1-000" ~ "null",
    type == "1-100" ~ "S",
    type == "1-010" ~ "T",
    type == "1-110" ~ "ST",
    TRUE ~ type
  )) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type) |> 
  # filter the params we work with
  tidylog::filter(par_name %in% c("log_kappaS", "kappaT", "AR(1)", "RMSD"))

# Filter true values (for horizontal lines)
pdattrue <- pdat |> 
  filter(type_of_par == "par_true") |> 
  distinct(value, par_name, species, type)

# Filter estimated parameters
pdathat <- pdat |> 
  filter(type_of_par == "par_hat") |> 
  mutate(true_model = ifelse(true_model == TRUE, "TRUE", "FALSE")) 

# Add in rmsd (for making different symbols)
rmsd <- pdathat |> 
  filter(par_name == "RMSD" & type %in% c("ST", "S")) |> 
  dplyr::select(iter, type, species, value) |> 
  rename(RMSD = value)
# 
# t <- rmsd |> tidylog::filter(RMSD < 1)
# t

# Join the rmsd to the main par df
pdathat <- pdathat |> 
  tidylog::left_join(rmsd, by = c("iter", "type", "species")) |> 
  mutate(shape_value = ifelse(RMSD < 1 & par_name %in% c("log_kappaS", "RMSD"),
                              "RMSD<1", "RMSD>1"),
         shape_value = replace_na(shape_value, "RMSD>1"))


pdathat <- pdathat |> 
  mutate(
    # Create a combined label for custom ordering
    facet_order = paste(par_name, species, sep = "_"),
    facet_order = factor(facet_order, 
                         levels = c("log_kappaS_capelin", 
                                    "kappaT_capelin", 
                                    "kappaT_pacific halibut",
                                    "RMSD_capelin", 
                                    "AR(1)_capelin", 
                                    "AR(1)_pacific halibut")))

# Do the same for pdattrue
pdattrue <- pdattrue |> 
  mutate(
    # Create a combined label for custom ordering
    facet_order = paste(par_name, species, sep = "_"),
    facet_order = factor(facet_order, 
                         levels = c("log_kappaS_capelin", 
                                    "kappaT_capelin", 
                                    "kappaT_pacific halibut",
                                    "RMSD_capelin", 
                                    "AR(1)_capelin", 
                                    "AR(1)_pacific halibut"))
  )

legend_col <- brewer.pal(3, "Dark2")[1:2]

p1 <- pdathat |> 
  tidylog::filter(!(species == "capelin" & par_name == "kappaT" & value > 200 )) |> 
  ggplot(aes(type, value)) +
  facet_wrap(~facet_order, 
             scales = "free", 
             ncol = 3,
             labeller = labeller(facet_order = as_labeller(
               c("log_kappaS_capelin" = "log(kappa[S])",
                 "kappaT_capelin" = "kappa[T]",
                 "kappaT_pacific halibut" = "kappa[T]",
                 #"RMSD_capelin" = "atop(Root-mean~squared, displacement~(km))",
                 "RMSD_capelin" = "Root*'\u2013'*mean*'\u2013'*squared~displacement~(km)",
                 "AR(1)_capelin" = "rho[t]",
                 "AR(1)_pacific halibut" = "rho[t]"),
               default = label_parsed
             ))
             ) +
  geom_hline(data = pdattrue, aes(yintercept = value), linetype = 2,
             color = "tomato", alpha = 0.7, linewidth = 0.6) +
  geom_quasirandom(aes(fill = true_model, shape = shape_value, color = shape_value), stroke = 0.5, 
                   alpha = 0.3, size = 1.5) +
  scale_shape_manual(values = c(4, 21), name = "") +
  scale_fill_brewer(palette = "Dark2", name = "operating model", direction = -1) +
  scale_color_manual(values = c("grey20", "grey95")) + 
  geom_boxplot(fill = NA, width = 0.2, size = 0.4,
               outlier.shape = NA, color = "gray20") +
  geom_text(data = . %>% distinct(species, facet_order),
            aes(label = species, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.5, size = 3, color = "gray30") +
  labs(y = "estimated value",
       x = "estimation model") +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.5, "cm"),
        strip.text.x.top = element_text(size = 9.5)) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(color = rev(legend_col))),
         color = "none",
         shape = guide_legend(ncol = 2, override.aes = list(color = "gray10")),
         linetype = guide_legend(title.position = "top", title.hjust = 0.5))

tag_facet(p1, fontface = 1, open = "", close = ")", color = "gray20", label = LETTERS[1:6])

ggsave(paste0(root_dir, "/figs/kappa_recovery.png"), width = 22, height = 17, unit = "cm")

