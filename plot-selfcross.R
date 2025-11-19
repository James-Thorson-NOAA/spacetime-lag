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
  theme(panel.spacing = unit(1, "lines"),
        plot.margin = margin(5, 10, 5, 5)) +
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
  mutate(shape_value = ifelse(rmsd_value < 1 & par_name %in% c("log_kappaS", "RMSD"),
                              "RMSD<1", "RMSD>1"))

p1 <- pdat2 |> 
  ggplot(aes(type, par_hat)) +
  facet_wrap(~ species_label + par_label_parsed, 
             scales = "free", 
             ncol = 3,
             labeller = labeller(par_label_parsed = label_parsed, 
                                 species_label = function(x) "")) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(aes(fill = true_model, shape = shape_value),
                   alpha = 0.8, size = 1, color = "gray70") +
  scale_shape_manual(values = c(4, 21), name = "") +
  scale_fill_manual(values = c("gray50", "gray99"), name = "") +
  geom_boxplot(fill = NA, width = 0.2, size = 0.4,
               outlier.shape = NA, color = "gray20") +
  geom_text(data = . %>% distinct(species_label, par_label_parsed),
            aes(label = species_label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.2, size = 3, color = "gray30") +
  labs(y = "Estimated value",
       x = "Estimation model") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(color = c("gray20", "gray90"))))

tag_facet(p1, fontface = 1, open = "", close = ")", color = "gray20")

ggsave(paste0(root_dir, "/figs/kappa_recovery.png"), width = 22, height = 17, unit = "cm")

