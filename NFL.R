library(nflreadr)
library(dplyr)
library(ggplot2)
library(tibble)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ==================================================================
# SETTINGS
# ==================================================================
seasons          <- 2025:2026       # more seasons are fine (e.g. 2023:2026)
include_playoffs <- FALSE
st_types         <- c("extra_point", "kickoff", "punt", "field_goal")
use_garbage_filter <- TRUE
wp_cutoff        <- 0.001

# Sampler settings. For a first test run use chains = 2, iter = 1000.
n_chains <- 4
n_iter   <- 2000

# ==================================================================
# DATA: one row per game (home perspective)
# ==================================================================
pbp <- load_pbp(seasons = seasons) %>% mutate(game_id = as.character(game_id))
if (!include_playoffs) pbp <- pbp %>% filter(season_type == "REG")
pbp$neutral_flag <- if ("location" %in% names(pbp)) pbp$location == "Neutral" else FALSE

sched <- pbp %>%
  group_by(game_id) %>%
  summarise(season       = first(season),
            home_team    = first(home_team),
            away_team    = first(away_team),
            neutral_site = isTRUE(first(neutral_flag)),
            home_points  = max(home_score, na.rm = TRUE),
            away_points  = max(away_score, na.rm = TRUE),
            .groups = "drop") %>%
  filter(is.finite(home_points), is.finite(away_points))

# Offense/defense: scrimmage plays only, special teams excluded
plays <- pbp %>%
  filter((rush == 1 | pass == 1), !play_type %in% st_types, !is.na(epa), !is.na(posteam))
st_plays <- pbp %>%
  filter(play_type %in% st_types, !is.na(epa), !is.na(posteam), !is.na(defteam))
if (use_garbage_filter) {
  gf <- function(d) d %>% filter(is.na(wp) | (wp > wp_cutoff & wp < 1 - wp_cutoff))
  plays <- gf(plays); st_plays <- gf(st_plays)
}
stopifnot(!any(plays$play_type %in% st_types))

off <- plays %>%
  group_by(game_id, team = posteam) %>%
  summarise(plays = n(), epa = mean(epa), .groups = "drop")

# ST: kicking team gets +epa, receiving team -epa; summed per team-game (points)
st_game <- bind_rows(
  st_plays %>% transmute(game_id, team = posteam, st_epa =  epa),
  st_plays %>% transmute(game_id, team = defteam, st_epa = -epa)
) %>% group_by(game_id, team) %>% summarise(st_net = sum(st_epa), .groups = "drop")

gdat <- sched %>%
  inner_join(off %>% rename(home_team = team, n_h = plays, epa_h = epa),
             by = c("game_id", "home_team")) %>%
  inner_join(off %>% rename(away_team = team, n_a = plays, epa_a = epa),
             by = c("game_id", "away_team")) %>%
  left_join(st_game %>% rename(home_team = team, st_h = st_net),
            by = c("game_id", "home_team")) %>%
  mutate(st_h    = coalesce(st_h, 0),
         poss    = (n_h + n_a) / 2,
         margin  = home_points - away_points,         # home margin
         hca     = as.numeric(!neutral_site))         # 0 for neutral-site games

teams    <- sort(unique(c(gdat$home_team, gdat$away_team)))
n_team   <- length(teams)
n_season <- length(seasons)

cat("\n=== DATA CHECKS ===\n")
cat("Games:", nrow(gdat), "| teams:", n_team, "| seasons:", paste(seasons, collapse = ", "), "\n")
print(gdat %>% group_by(season) %>%
        summarise(games = n(), neutral = sum(hca == 0), avg_plays = mean(poss),
                  avg_epa = mean(c(epa_h, epa_a)), avg_st_home = mean(st_h),
                  avg_margin = mean(margin), .groups = "drop"))
games_per_team <- table(c(gdat$home_team, gdat$away_team))
cat("Games per team (min, max):", min(games_per_team), max(games_per_team), "\n")
stopifnot(!anyNA(gdat[, c("epa_h", "epa_a", "n_h", "n_a", "st_h", "poss", "margin")]))

stan_data <- list(
  G = nrow(gdat), n_team = n_team, n_season = n_season,
  home   = match(gdat$home_team, teams),
  away   = match(gdat$away_team, teams),
  season = match(gdat$season, seasons),
  hca    = gdat$hca,
  epa_h  = gdat$epa_h, epa_a = gdat$epa_a,
  n_h    = gdat$n_h,   n_a   = gdat$n_a,
  st_h   = gdat$st_h,
  poss   = gdat$poss,
  margin = gdat$margin
)

# ==================================================================
# STAN MODEL
#   EPA/play (home):  mu + off[home] - def[away] + h*hca        sd = sigma_e/sqrt(plays)
#   EPA/play (away):  mu + off[away] - def[home] - h*hca
#   ST points (home): st[home] - st[away]                       sd = sigma_st
#   Plays per team:   mu_p + tempo[home] + tempo[away]
# off/def/st are team-season effects with a hierarchical Normal prior and an
# AR(1) link between seasons (rho_*), which replaces prior_season_weight and
# the fixed ridge: the data now estimate how much last year carries over and
# how spread out true team talent is (tau_*), i.e. the amount of shrinkage.
# "def" is defensive QUALITY (higher = better), so a team's net EPA edge is
# off + def.  Points per game = plays * EPA edge (+ ST points).
# ==================================================================
stan_code <- "
data {
  int<lower=1> G;
  int<lower=1> n_team;
  int<lower=1> n_season;
  array[G] int<lower=1, upper=n_team> home;
  array[G] int<lower=1, upper=n_team> away;
  array[G] int<lower=1, upper=n_season> season;
  vector<lower=0, upper=1>[G] hca;
  vector[G] epa_h;
  vector[G] epa_a;
  vector<lower=1>[G] n_h;
  vector<lower=1>[G] n_a;
  vector[G] st_h;
  vector<lower=1>[G] poss;
  vector[G] margin;
}
parameters {
  real mu;
  real h;
  real<lower=0> sigma_e;
  real<lower=0> tau_o;
  real<lower=0> tau_d;
  real<lower=0> tau_s;
  real<lower=0, upper=1> rho_o;
  real<lower=0, upper=1> rho_d;
  real<lower=0, upper=1> rho_s;
  matrix[n_team, n_season] z_o;
  matrix[n_team, n_season] z_d;
  matrix[n_team, n_season] z_s;
  real<lower=0> sigma_st;
  real mu_p;
  real<lower=0> tau_p;
  real<lower=0> sigma_p;
  vector[n_team] z_p;
}
transformed parameters {
  matrix[n_team, n_season] off;
  matrix[n_team, n_season] dfn;
  matrix[n_team, n_season] st;
  vector[n_team] tempo = tau_p * z_p;
  for (t in 1:n_team) {
    off[t, 1] = tau_o * z_o[t, 1];
    dfn[t, 1] = tau_d * z_d[t, 1];
    st[t, 1]  = tau_s * z_s[t, 1];
    for (s in 2:n_season) {
      off[t, s] = rho_o * off[t, s - 1] + tau_o * sqrt(1 - square(rho_o)) * z_o[t, s];
      dfn[t, s] = rho_d * dfn[t, s - 1] + tau_d * sqrt(1 - square(rho_d)) * z_d[t, s];
      st[t, s]  = rho_s * st[t, s - 1]  + tau_s * sqrt(1 - square(rho_s)) * z_s[t, s];
    }
  }
}
model {
  // priors
  mu       ~ normal(0, 0.05);
  h        ~ normal(0.015, 0.02);       // EPA/play home edge per side (~1.7 pts total)
  sigma_e  ~ normal(0, 3);
  tau_o    ~ normal(0, 0.1);
  tau_d    ~ normal(0, 0.1);
  tau_s    ~ normal(0.93, 0.3);         // informed by the many-season ST estimate
  rho_o    ~ beta(3, 3);
  rho_d    ~ beta(3, 3);
  rho_s    ~ beta(3, 3);
  sigma_st ~ normal(5, 2);              // many-season estimate was ~4.96
  mu_p     ~ normal(58, 10);
  tau_p    ~ normal(0, 5);
  sigma_p  ~ normal(0, 10);
  to_vector(z_o) ~ std_normal();
  to_vector(z_d) ~ std_normal();
  to_vector(z_s) ~ std_normal();
  z_p ~ std_normal();

  // likelihood
  for (g in 1:G) {
    int i = home[g];
    int j = away[g];
    int s = season[g];
    epa_h[g] ~ normal(mu + off[i, s] - dfn[j, s] + h * hca[g], sigma_e / sqrt(n_h[g]));
    epa_a[g] ~ normal(mu + off[j, s] - dfn[i, s] - h * hca[g], sigma_e / sqrt(n_a[g]));
    st_h[g]  ~ normal(st[i, s] - st[j, s], sigma_st);
    poss[g]  ~ normal(mu_p + tempo[i] + tempo[j], sigma_p);
  }
}
generated quantities {
  vector[n_team] off_pts;
  vector[n_team] def_pts;
  vector[n_team] st_pts;
  vector[n_team] strength_pts;   // points/game vs an average team, current (last) season
  vector[G] exp_margin;          // model-expected home margin for every game played
  real hfa_bias;
  real sigma_game;
  {
    vector[G] resid;
    real bias_sum = 0;
    real n_home = 0;
    for (t in 1:n_team) {
      off_pts[t] = mu_p * off[t, n_season];
      def_pts[t] = mu_p * dfn[t, n_season];
      st_pts[t]  = st[t, n_season];
      strength_pts[t] = off_pts[t] + def_pts[t] + st_pts[t];
    }
    for (g in 1:G) {
      int i = home[g];
      int j = away[g];
      int s = season[g];
      real poss_e = mu_p + tempo[i] + tempo[j];
      exp_margin[g] = poss_e * ((off[i, s] + dfn[i, s]) - (off[j, s] + dfn[j, s]) + 2 * h * hca[g])
                      + st[i, s] - st[j, s];
      resid[g] = margin[g] - exp_margin[g];
      if (hca[g] > 0.5) {
        bias_sum += resid[g];
        n_home += 1;
      }
    }
    hfa_bias   = bias_sum / n_home;                       // home margin the EPA model misses
    sigma_game = sqrt(mean(square(resid - hca * hfa_bias)));  // game-level noise (points)
  }
}
"

# ==================================================================
# FIT
# ==================================================================
stan_fit <- stan(model_code = stan_code, data = stan_data,
                 chains = n_chains, iter = n_iter, warmup = n_iter / 2,
                 seed = 2026, control = list(adapt_delta = 0.95))

cat("\n=== SAMPLER DIAGNOSTICS ===\n")
rstan::check_hmc_diagnostics(stan_fit)
hyper <- c("mu", "h", "sigma_e", "tau_o", "tau_d", "tau_s", "rho_o", "rho_d", "rho_s",
           "sigma_st", "mu_p", "tau_p", "sigma_p", "hfa_bias", "sigma_game")
hs <- summary(stan_fit, pars = hyper)$summary
print(round(hs[, c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], 3))
cat("Max Rhat over all parameters:",
    round(max(summary(stan_fit)$summary[, "Rhat"], na.rm = TRUE), 3), "(want < 1.01)\n")
cat("Min effective sample size:",
    round(min(summary(stan_fit)$summary[, "n_eff"], na.rm = TRUE)), "\n")

ext <- rstan::extract(stan_fit,
                      pars = c("off", "dfn", "st", "tempo", "mu_p", "h", "hfa_bias",
                               "sigma_game", "strength_pts", "off_pts", "def_pts",
                               "st_pts", "exp_margin"))
n_draws <- length(ext$mu_p)

# ------------------------------------------------------------------
# Model checks
# ------------------------------------------------------------------
em_mean <- colMeans(ext$exp_margin)
cal <- lm(gdat$margin ~ 0 + em_mean)
cat("\n=== CALIBRATION CHECK ===\n")
cat("Slope of actual margin on model-expected margin:", round(coef(cal), 3),
    "(1 = calibrated; in-sample so mildly optimistic)\n")
cat("Median game-level SD (sigma_game):", round(median(ext$sigma_game), 2), "pts\n")
cat("Median home margin missed by EPA model (hfa_bias):", round(median(ext$hfa_bias), 2), "pts\n")
cat("Mean strength across teams (should be ~0):", round(mean(ext$strength_pts), 3), "\n")

# ==================================================================
# RANKING WITH UNCERTAINTY
# ==================================================================
draw_summary <- function(M) tibble(
  team = teams,
  med  = apply(M, 2, median),
  lo50 = apply(M, 2, quantile, 0.25), hi50 = apply(M, 2, quantile, 0.75),
  lo90 = apply(M, 2, quantile, 0.05), hi90 = apply(M, 2, quantile, 0.95))

ranks <- t(apply(-ext$strength_pts, 1, rank, ties.method = "min"))   # draws x teams

cur <- max(seasons)
rec <- bind_rows(
  gdat %>% filter(season == cur) %>% transmute(team = home_team, w = margin > 0, l = margin < 0),
  gdat %>% filter(season == cur) %>% transmute(team = away_team, w = margin < 0, l = margin > 0)
) %>% group_by(team) %>%
  summarise(record = paste0(sum(w), "-", sum(l)), .groups = "drop")

ranking <- draw_summary(ext$strength_pts) %>%
  mutate(rank_med = apply(ranks, 2, median),
         rank_lo  = apply(ranks, 2, quantile, 0.05),
         rank_hi  = apply(ranks, 2, quantile, 0.95),
         p_top8   = colMeans(ranks <= 8),
         sd       = apply(ext$strength_pts, 2, sd),
         off      = apply(ext$off_pts, 2, median),
         def      = apply(ext$def_pts, 2, median),
         st       = apply(ext$st_pts, 2, median)) %>%
  left_join(rec, by = "team") %>%
  arrange(desc(med)) %>%
  mutate(rank = row_number())

options(pillar.sigfig = 3)
cat("\n=== TEAM STRENGTH (points/game vs average team, neutral field; 90% intervals) ===\n")
print(ranking %>% dplyr::select(rank, team, record, strength = med, sd, lo90, hi90,
                                rank_lo, rank_hi, p_top8, off, def, st), n = Inf, width = Inf)

# 1) Teams ranked by strength, with 50% (thick) and 90% (thin) intervals
ggplot(ranking, aes(x = reorder(team, med), y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lo90, ymax = hi90), linewidth = 0.6, alpha = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.8) +
  geom_point(size = 2) +
  coord_flip() +
  labs(x = NULL, y = "Strength (points per game vs average team)",
       title = paste0("Team strength, ", cur, " season"),
       subtitle = "Posterior median; thick = 50% interval, thin = 90% interval") +
  theme_minimal()

# 2) Rank intervals
ggplot(ranking, aes(x = reorder(team, -rank_med), y = rank_med)) +
  geom_linerange(aes(ymin = rank_lo, ymax = rank_hi), linewidth = 1) +
  geom_point(size = 2) +
  scale_y_reverse(breaks = seq(1, n_team, by = 4)) +
  coord_flip() +
  labs(x = NULL, y = "Rank (1 = best)",
       title = "How sure are we of each team's rank?",
       subtitle = "Posterior median rank and 90% interval") +
  theme_minimal()

# 3) Components: offense, defense, special teams (points per game)
comp <- bind_rows(draw_summary(ext$off_pts) %>% mutate(component = "Offense"),
                  draw_summary(ext$def_pts) %>% mutate(component = "Defense"),
                  draw_summary(ext$st_pts)  %>% mutate(component = "Special teams")) %>%
  mutate(team = factor(team, levels = rev(ranking$team)),
         component = factor(component, levels = c("Offense", "Defense", "Special teams")))
ggplot(comp, aes(x = team, y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lo90, ymax = hi90), linewidth = 0.5, alpha = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.4) +
  geom_point(size = 1.4) +
  coord_flip() + facet_wrap(~component) +
  labs(x = NULL, y = "Points per game vs average",
       title = "Strength by component (teams ordered by total)") +
  theme_minimal() + theme(axis.text.y = element_text(size = 7))

# ==================================================================
# MATCHUP: point estimate + uncertainty for the margin
# ==================================================================
matchup <- function(home, away, neutral = FALSE, market_spread = NULL, level = 0.90) {
  stopifnot(home %in% teams, away %in% teams, home != away)
  ih  <- match(home, teams); ia <- match(away, teams)
  hca <- as.numeric(!neutral)
  S   <- n_season
  
  net    <- function(i) ext$off[, i, S] + ext$dfn[, i, S]          # EPA/play edge
  poss_d <- ext$mu_p + ext$tempo[, ih] + ext$tempo[, ia]            # expected plays per team
  
  c_ed  <- poss_d * (net(ih) - net(ia))                             # offense + defense
  c_st  <- ext$st[, ih, S] - ext$st[, ia, S]                        # special teams
  c_hfa <- (poss_d * 2 * ext$h + ext$hfa_bias) * hca                # home field
  line  <- c_ed + c_st + c_hfa                                      # expected home margin
  
  reps <- 20                                                        # game-noise draws per posterior draw
  game <- rep(line, each = reps) + rnorm(length(line) * reps, 0, rep(ext$sigma_game, each = reps))
  
  p <- c((1 - level) / 2, 1 - (1 - level) / 2)
  qi <- function(x) quantile(x, p)
  
  cat(sprintf("\n=== %s (home) vs %s (away)%s ===\n", home, away,
              if (neutral) " [neutral site]" else ""))
  cat(sprintf("Expected home margin: %+.1f  |  model spread for %s: %+.1f\n",
              mean(line), home, -mean(line)))
  cat(sprintf("Uncertainty in the expected margin (parameters only): SD %.2f, %d%% interval [%+.1f, %+.1f]\n",
              sd(line), level * 100, qi(line)[1], qi(line)[2]))
  cat(sprintf("  SD by source: off/def %.2f | special teams %.2f | home field %.2f\n",
              sd(c_ed), sd(c_st), sd(c_hfa)))
  cat(sprintf("Predicted game result (incl. game noise): SD %.1f, %d%% interval [%+.1f, %+.1f]\n",
              sd(game), level * 100, qi(game)[1], qi(game)[2]))
  cat(sprintf("Win probability: %s %.1f%% | %s %.1f%%\n",
              home, 100 * mean(pnorm(line / ext$sigma_game)),
              away, 100 * (1 - mean(pnorm(line / ext$sigma_game)))))
  cat(sprintf("P(%s is the better team on a neutral field): %.1f%%\n", home,
              100 * mean(c_ed + c_st > 0)))
  if (!is.null(market_spread)) {
    cat(sprintf("Market spread for %s %+.1f: P(%s covers) = %.1f%% | P(expected margin beats market) = %.1f%%\n",
                home, market_spread, home, 100 * mean(game > -market_spread),
                100 * mean(line > -market_spread)))
  }
  
  dd <- bind_rows(tibble(value = line, what = "Expected margin (parameter uncertainty only)"),
                  tibble(value = game, what = "Game result (parameters + game noise)"))
  g <- ggplot(dd, aes(x = value, fill = what)) +
    geom_density(alpha = 0.45, color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = mean(line), color = "black") +
    labs(x = paste0("Home margin (", home, " minus ", away, ")"), y = NULL, fill = NULL,
         title = paste0(home, " vs ", away, if (neutral) " (neutral site)" else ""),
         subtitle = sprintf("Expected margin %+.1f (spread %+.1f); win prob %.0f%%",
                            mean(line), -mean(line),
                            100 * mean(pnorm(line / ext$sigma_game)))) +
    theme_minimal() + theme(legend.position = "bottom", axis.text.y = element_blank())
  if (!is.null(market_spread)) g <- g + geom_vline(xintercept = -market_spread, color = "red")
  print(g)
  invisible(list(line = line, game = game,
                 parts = tibble(offdef = c_ed, st = c_st, hfa = c_hfa)))
}

# ---- Choose a matchup ----
home_team     <- "GB"
away_team     <- "DAL"
neutral_game  <- FALSE
market_spread <- NULL        # e.g. -3.5 if the home team is a 3.5-point favourite
use_prompt    <- FALSE       # TRUE = ask in the console (interactive sessions only)

if (use_prompt && interactive()) {
  cat("Teams:", paste(teams, collapse = " "), "\n")
  home_team <- toupper(readline("Home team: "))
  away_team <- toupper(readline("Away team: "))
  neutral_game <- tolower(readline("Neutral site? (y/n): ")) == "y"
  ms <- readline("Market spread for the home team (blank = none): ")
  market_spread <- if (nzchar(ms)) as.numeric(ms) else NULL
}

res <- matchup(home_team, away_team, neutral = neutral_game, market_spread = market_spread)
