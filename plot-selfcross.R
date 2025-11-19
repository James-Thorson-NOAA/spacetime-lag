library(egg)
library(ggbeeswarm)
library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(ggsidekick); theme_set(theme_sleek())

root_dir <- here::here()

# plot
n_sim <- 100

res_df <- read_csv(file.path(root_dir, paste0("2025-11-09/simulation_selfcross_out_", n_sim, "_reps.csv")))

# specify true model for each species
correct_models <- tibble(
  species = c("capelin", "pacific cod", "pacific halibut"),
  correct_type = c("space+time", "base", "time")
)

res_df |> 
  distinct(AIC, type, iter, species) |> 
  mutate(type = case_when(
    type == "1-000" ~ "base",
    type == "1-100" ~ "space",
    type == "1-010" ~ "time",
    type == "1-110" ~ "space+time",
    TRUE ~ type
  )) |>
  mutate(delta_AIC = AIC - min(AIC), .by = c(iter, species)) |> 
  summarise(fraction_best = mean(delta_AIC == 0), .by = c(type, species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type) |>
  ggplot(aes(fraction_best, type, fill = true_model)) +
  facet_wrap(~str_to_sentence(species), ncol = 3) +
  theme(legend.position = "bottom") +
  geom_col(width = 0.9) +
  labs(
    y = "Model",
    x = "Fraction of iterations with best AIC"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Operating model",
                    direction = -1)

ggsave(paste0(root_dir, "/figs/self_test_aic.png"), width = 20, height = 8, unit = "cm")


# Same plot but add AR(1) [rho_E] and RMSD (from log_kappaS)
pdat <- res_df |> 
  filter(par_name %in% c("log_kappaS", "kappaT", "logit_rhoE", "RMSD")) |> 
  filter(species %in% c("capelin", "pacific halibut")) |> 
  drop_na(par_true) |> 
  mutate(type = case_when(
    type == "1-000" ~ "base",
    type == "1-100" ~ "space",
    type == "1-010" ~ "time",
    type == "1-110" ~ "space+time",
    TRUE ~ type
  )) |>
  mutate(par_label_parsed = case_when(
    par_name == "kappaT" ~ "kappa[T]",
    par_name == "log_kappaS" ~ "log(kappa[S])",
    par_name == "logit_rhoE" ~ "rho[E]",
    TRUE ~ par_name
  ),
  species_label = str_to_sentence(species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type,
         par_hat = ifelse(par_name == "logit_rhoE",
                          2 * plogis(par_hat) - 1,
                          par_hat),
         par_true = ifelse(par_name == "logit_rhoE",
                           2 * plogis(par_true) - 1,
                           par_true))

#mutate(RMSD = (4 / exp( 2.0 * log_kappaS))^2) |> 
RMSD <- res_df |> 
  filter(par_name == "log_kappaS") |> 
  filter(species %in% c("capelin", "pacific halibut")) |> 
  mutate(par_hat = sqrt(4 / exp( 2.0 * par_hat)),
         par_true = sqrt(4 / exp( 2.0 * par_true)),
         par_name = "RMSD") |> 
  drop_na(par_true) |> 
  mutate(type = case_when(
    type == "1-000" ~ "base",
    type == "1-100" ~ "space",
    type == "1-010" ~ "time",
    type == "1-110" ~ "space+time",
    TRUE ~ type
  )) |>
  mutate(par_label_parsed = case_when(
    par_name == "RMSD" ~ "RMSD",
    TRUE ~ par_name
  ),
  species_label = str_to_sentence(species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type)

# Now add in RMSD
pdat2 <- pdat |>
  bind_rows(RMSD)

pdat2 <- pdat2 |> 
  mutate(par_label_parsed =
           factor(par_label_parsed, 
                  levels = c("log(kappa[S])", "kappa[T]", "rho[E]", "RMSD"))) |>
  mutate(true_model = ifelse(true_model == TRUE, "AIC-selected", "Not AIC-selected")) |> 
  mutate(
    rmsd_value = ifelse(any(par_name == "RMSD"), 
                        par_hat[par_name == "RMSD"][1], 
                        NA_real_),
    log_kappaS_value = ifelse(any(par_name == "log_kappaS"), 
                              par_hat[par_name == "log_kappaS"][1], 
                              NA_real_),
    .by = c(species, type, iter, seed)
  ) |>
  ungroup() |>
  mutate(color_value = case_when(
    par_name == "log_kappaS" ~ scale(rmsd_value)[,1],
    par_name == "RMSD" ~ scale(log_kappaS_value)[,1],
    TRUE ~ NA_real_
  ))

p1 <- pdat2 |> 
  ggplot(aes(type, par_hat)) +
  facet_wrap(~ species_label + par_label_parsed, 
             scales = "free", 
             ncol = 3,
             labeller = labeller(par_label_parsed = label_parsed, 
                                 species_label = function(x) "")) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(data = . %>% filter(is.na(color_value)),
                   aes(shape = true_model), alpha = 0.3, size = 1.1, color = "gray50") +
  geom_quasirandom(data = . %>% filter(!is.na(color_value)), color = "gray50", stroke = 0.1,
                   aes(fill = color_value, shape = true_model), alpha = 1, size = 1.1) +
  # scale_color_viridis_c(name = expression(log(kappa[S])~or~RMSD (scaled)), 
  #                       option = "plasma") +
  scale_fill_gradient2(name = expression(log(kappa[S])~or~RMSD (scaled))) +
  scale_color_brewer(palette = "Dark2", name = "") +
  scale_shape_manual(values = c("AIC-selected" = 21, "Not AIC-selected" = 23), name = "") +
  # geom_boxplot(fill = NA, width = 0.2, size = 0.4,
  #              outlier.shape = NA, color = "gray20") +
  geom_boxplot(
    aes(color = true_model),
    fill = NA,
    width = 0.15,
    size = 0.3,
    outlier.shape = NA
  ) +
  geom_text(data = . %>% distinct(species_label, par_label_parsed),
            aes(label = species_label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.2, size = 4, color = "gray30") +
  labs(y = "Estimated value",
       x = "Estimation model") +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.2, "cm")) +
  guides(shape = "none",
         #shape = guide_legend(ncol = 1, override.aes = list(alpha = 0.9, size = 3)),
         color = guide_legend(ncol = 1, override.aes = list(width = 0)))

tag_facet(p1, fontface = 1, open = "", close = "", color = "gray20")
ggsave(paste0(root_dir, "/figs/kappa_recovery.png"), width = 22, height = 17, unit = "cm")



# Splitting the plot to make life easier 
pdat_kappa <- pdat2 |> 
  filter(par_name %in% c("log_kappaS", "RMSD"))

pdat_other <- pdat2 |> 
  filter(par_name %in% c("kappaT", "logit_rhoE"))

# log_kappaS + RMSD
p2 <- pdat_kappa |> 
  ggplot(aes(type, par_hat)) +
  facet_wrap(~ species_label + par_label_parsed, 
             scales = "free", ncol = 2,
             labeller = labeller(par_label_parsed = label_parsed,
                                 species_label = function(x) "")) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(aes(shape = true_model, fill = color_value), 
                   alpha = 1, size = 1.1, color = "gray85") +
  scale_shape_manual(values = c("AIC-selected" = 21, "Not AIC-selected" = 23), name = "") +
  scale_fill_gradient2(name = expression(log(kappa[S])~or~RMSD (scaled))) +
  geom_boxplot(fill = NA, color = "black", width = 0.15, size = 0.3, outlier.shape = NA) +
  labs(y = "Estimated value", x = "Estimation model") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
         shape = guide_legend(ncol = 1, override.aes = list(alpha = 1, color = "grey50", size = 3))) + 
  theme(legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.2, "cm"))

tag_facet(p2, fontface = 1, open = "", close = "", color = "gray20")

ggsave(paste0(root_dir, "/figs/kappa_recovery_logkappa_RMSD.png"), 
       width = 20, height = 12, units = "cm")


# kappaT + rhoE
p_other <- pdat_other |> 
  mutate(true_model = ifelse(true_model, "AIC-selected", "Not AIC-selected")) |>
  ggplot(aes(type, par_hat)) +
  facet_wrap(~ species_label + par_label_parsed, 
             scales = "free", ncol = 2,
             labeller = labeller(par_label_parsed = label_parsed,
                                 species_label = function(x) "")) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(aes(shape = true_model, alpha = true_model), size = 0.6, color = "gray10") +
  geom_boxplot(fill = NA, color = "black", width = 0.15, size = 0.3, outlier.shape = NA) +
  scale_shape_manual(values = c("AIC-selected" = 19, "Not AIC-selected" = 1), name = "") +
  scale_alpha_manual(values = c(0.2, 0.15), name = "") +
  labs(y = "Estimated value", x = "Estimation model") +
  #coord_flip() +
  theme(legend.position = "bottom")

ggsave(paste0(root_dir, "/figs/kappa_recovery_kappaT_rhoE.png"),
       width = 17, height = 14, units = "cm")

