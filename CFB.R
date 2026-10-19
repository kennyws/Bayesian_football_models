# ==================================================================
# COLLEGE FOOTBALL HIERARCHICAL TEAM STRENGTH MODEL
# Bayesian offense / defense / special teams / tempo model using cfbfastR
# ==================================================================

library(cfbfastR)
library(dplyr)
library(ggplot2)
library(tibble)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ==================================================================
# SETTINGS
# ==================================================================
# Use completed seasons plus the current season if desired. The current
# season is partial until the season ends, so its ratings will be noisier.
seasons <- 2025:2026
include_playoffs <- FALSE
st_pattern <- "kickoff|punt|field goal|extra point"
use_garbage_filter <- TRUE
wp_cutoff <- 0.001

# Remove opponents that appear in only a small number of games across the
# selected seasons. In FBS play-by-play, FCS opponents can appear in the
# schedule but usually have far fewer appearances than FBS teams.
min_team_games <- 5

# Stan settings. Start with 2 chains and 1000 iterations to test the script.
n_chains <- 2
n_iter <- 1000

# Matchup simulation settings. A few game-noise simulations per posterior
# draw are sufficient; posterior uncertainty is already represented by draws.
game_reps <- 5

# ==================================================================
# LOAD AND STANDARDIZE PLAY-BY-PLAY
# ==================================================================
pbp <- cfbfastR::load_cfb_pbp(seasons = seasons) %>%
  mutate(game_id = as.character(game_id))

# Standardize names across cfbfastR data versions.
if ("EPA" %in% names(pbp) && !"epa" %in% names(pbp)) {
  pbp <- rename(pbp, epa = EPA)
}
if ("wp_before" %in% names(pbp) && !"wp" %in% names(pbp)) {
  pbp <- rename(pbp, wp = wp_before)
}
if ("posteam" %in% names(pbp) && !"pos_team" %in% names(pbp)) {
  pbp <- rename(pbp, pos_team = posteam)
}
if ("defteam" %in% names(pbp) && !"def_pos_team" %in% names(pbp)) {
  pbp <- rename(pbp, def_pos_team = defteam)
}
if (!"season" %in% names(pbp) && "year" %in% names(pbp)) {
  pbp <- rename(pbp, season = year)
}
if (!"home_team" %in% names(pbp) && "home" %in% names(pbp)) {
  pbp <- rename(pbp, home_team = home)
}
if (!"away_team" %in% names(pbp) && "away" %in% names(pbp)) {
  pbp <- rename(pbp, away_team = away)
}

# ------------------------------------------------------------------
# Restrict the model to FBS teams; pool non-FBS opponents as "FCS"
# ------------------------------------------------------------------

fbs_teams <- bind_rows(
  pbp %>%
    filter(tolower(home_team_division) == "fbs") %>%
    distinct(team = home_team),
  
  pbp %>%
    filter(tolower(away_team_division) == "fbs") %>%
    distinct(team = away_team)
) %>%
  distinct(team) %>%
  pull(team)

pbp <- pbp %>%
  # Keep games involving at least one FBS team
  filter(home_team %in% fbs_teams | away_team %in% fbs_teams) %>%
  # Pool every non-FBS team into a single opponent category
  mutate(
    across(
      c(home_team, away_team, pos_team, def_pos_team),
      ~ if_else(.x %in% fbs_teams, .x, "FCS")
    )
  )

# This cfbfastR version uses pos_team_score / def_pos_team_score rather
# than home_score / away_score. Convert each play's score to home/away
# perspective below.
required_cols <- c(
  "game_id", "season", "home_team", "away_team",
  "pos_team", "def_pos_team", "pos_team_score", "def_pos_team_score",
  "epa", "play_type"
)
missing_cols <- setdiff(required_cols, names(pbp))
if (length(missing_cols) > 0) {
  stop(
    "Required columns were not found in load_cfb_pbp(): ",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are:\n", paste(names(pbp), collapse = ", ")
  )
}

if (!include_playoffs && "season_type" %in% names(pbp)) {
  pbp <- pbp %>%
    filter(tolower(as.character(season_type)) %in%
             c("reg", "regular", "regular season"))
}

# Normalize play types because this data version uses labels such as
# "Field Goal" and "Extra Point", not snake_case labels.
pbp <- pbp %>%
  mutate(
    play_type_norm = tolower(gsub("_", " ", as.character(play_type))),
    is_st_play = grepl(st_pattern, play_type_norm),
    home_score_row = case_when(
      pos_team == home_team ~ as.numeric(pos_team_score),
      def_pos_team == home_team ~ as.numeric(def_pos_team_score),
      TRUE ~ NA_real_
    ),
    away_score_row = case_when(
      pos_team == away_team ~ as.numeric(pos_team_score),
      def_pos_team == away_team ~ as.numeric(def_pos_team_score),
      TRUE ~ NA_real_
    ),
    neutral_flag = if ("neutral_site" %in% names(.)) {
      !is.na(neutral_site) & neutral_site
    } else if ("location" %in% names(.)) {
      tolower(as.character(location)) %in% c("neutral", "neutral site")
    } else {
      FALSE
    }
  )

# ==================================================================
# ONE ROW PER GAME, HOME TEAM PERSPECTIVE
# ==================================================================
sched <- pbp %>%
  group_by(game_id) %>%
  summarise(
    season = first(season),
    home_team = first(home_team),
    away_team = first(away_team),
    neutral_site = any(neutral_flag %in% TRUE, na.rm = TRUE),
    home_points = if (all(is.na(home_score_row))) NA_real_ else max(home_score_row, na.rm = TRUE),
    away_points = if (all(is.na(away_score_row))) NA_real_ else max(away_score_row, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(home_points), is.finite(away_points),
    !is.na(home_team), !is.na(away_team), home_team != away_team
  )

# Identify teams with enough appearances to be plausibly FBS programs.
# This is a schedule-frequency filter, not a formal division lookup.
team_appearances <- bind_rows(
  sched %>% transmute(season, team = home_team),
  sched %>% transmute(season, team = away_team)
) %>%
  count(team, name = "n_games")

fbs_like_teams <- team_appearances %>%
  filter(n_games >= min_team_games) %>%
  pull(team)

sched <- sched %>%
  filter(home_team %in% fbs_like_teams, away_team %in% fbs_like_teams)

# ==================================================================
# OFFENSE / DEFENSE: SCRIMMAGE PLAYS ONLY
# ==================================================================
if (!all(c("rush", "pass") %in% names(pbp))) {
  stop("The play-by-play data must contain rush and pass indicators. Available columns: ",
       paste(names(pbp), collapse = ", "))
}

plays <- pbp %>%
  filter(
    (rush == 1 | pass == 1),
    !is_st_play,
    !is.na(epa), !is.na(pos_team), !is.na(def_pos_team),
    pos_team %in% fbs_like_teams,
    def_pos_team %in% fbs_like_teams
  )

# ==================================================================
# SPECIAL TEAMS
# ==================================================================
st_plays <- pbp %>%
  filter(
    is_st_play,
    !is.na(epa), !is.na(pos_team), !is.na(def_pos_team),
    pos_team %in% fbs_like_teams,
    def_pos_team %in% fbs_like_teams
  )

if (use_garbage_filter && "wp" %in% names(plays)) {
  gf <- function(d) {
    d %>% filter(is.na(wp) | (wp > wp_cutoff & wp < 1 - wp_cutoff))
  }
  plays <- gf(plays)
  st_plays <- gf(st_plays)
}

stopifnot(!any(plays$is_st_play))

off <- plays %>%
  group_by(game_id, team = pos_team) %>%
  summarise(n_plays = n(), epa = mean(epa), .groups = "drop")

# EPA is credited to the kicking team and debited from the opponent.
# Aggregate to a team-game total, then construct home-minus-away ST EPA.
st_game <- bind_rows(
  st_plays %>% transmute(game_id, team = pos_team, st_epa = epa),
  st_plays %>% transmute(game_id, team = def_pos_team, st_epa = -epa)
) %>%
  group_by(game_id, team) %>%
  summarise(st_net = sum(st_epa, na.rm = TRUE), .groups = "drop")

# ==================================================================
# JOIN GAME-LEVEL DATA
# ==================================================================
gdat <- sched %>%
  inner_join(
    off %>% rename(home_team = team, n_h = n_plays, epa_h = epa),
    by = c("game_id", "home_team")
  ) %>%
  inner_join(
    off %>% rename(away_team = team, n_a = n_plays, epa_a = epa),
    by = c("game_id", "away_team")
  ) %>%
  left_join(
    st_game %>% rename(home_team = team, st_home = st_net),
    by = c("game_id", "home_team")
  ) %>%
  left_join(
    st_game %>% rename(away_team = team, st_away = st_net),
    by = c("game_id", "away_team")
  ) %>%
  mutate(
    st_h = coalesce(st_home, 0) - coalesce(st_away, 0),
    poss = (n_h + n_a) / 2,
    margin = home_points - away_points,
    hca = as.numeric(!neutral_site)
  ) %>%
  filter(season %in% seasons) %>%
  filter(is.finite(epa_h), is.finite(epa_a), is.finite(poss), poss > 0)

if (nrow(gdat) == 0) {
  stop("No games remain after filtering. Check the play-by-play schema and min_team_games.")
}

teams <- sort(unique(c(gdat$home_team, gdat$away_team)))
n_team <- length(teams)
n_season <- length(seasons)

cat("\n=== DATA CHECKS ===\n")
cat("Games:", nrow(gdat), "| teams:", n_team,
    "| seasons:", paste(seasons, collapse = ", "), "\n")
print(
  gdat %>% group_by(season) %>%
    summarise(
      games = n(), neutral = sum(hca == 0), avg_plays_per_team = mean(poss),
      avg_epa = mean(c(epa_h, epa_a)), avg_st_home_minus_away = mean(st_h),
      avg_margin = mean(margin), .groups = "drop"
    )
)
games_per_team <- table(c(gdat$home_team, gdat$away_team))
cat("Games per team (min, max):", min(games_per_team), max(games_per_team), "\n")
stopifnot(!anyNA(gdat[, c("epa_h", "epa_a", "n_h", "n_a", "st_h", "poss", "margin")]))

stan_data <- list(
  G = nrow(gdat), n_team = n_team, n_season = n_season,
  home = match(gdat$home_team, teams),
  away = match(gdat$away_team, teams),
  season = match(gdat$season, seasons),
  hca = gdat$hca,
  epa_h = gdat$epa_h, epa_a = gdat$epa_a,
  n_h = gdat$n_h, n_a = gdat$n_a,
  st_h = gdat$st_h,
  poss = gdat$poss,
  margin = gdat$margin
)

# ==================================================================
# STAN MODEL
# Offense/defense effects are EPA per play. Higher defense effect means
# better defense. Special teams are modeled on the summed EPA scale.
# Team-season effects evolve through an AR(1) relationship across seasons.
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
  mu       ~ normal(0, 0.05);
  h        ~ normal(0.02, 0.025);
  sigma_e  ~ normal(0, 3);
  tau_o    ~ normal(0, 0.15);
  tau_d    ~ normal(0, 0.15);
  tau_s    ~ normal(0, 3);
  rho_o    ~ beta(3, 3);
  rho_d    ~ beta(3, 3);
  rho_s    ~ beta(3, 3);
  sigma_st ~ normal(0, 8);
  mu_p     ~ normal(65, 12);
  tau_p    ~ normal(0, 6);
  sigma_p  ~ normal(0, 12);
  to_vector(z_o) ~ std_normal();
  to_vector(z_d) ~ std_normal();
  to_vector(z_s) ~ std_normal();
  z_p ~ std_normal();

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
  vector[n_team] strength_pts;
  vector[G] exp_margin;
  real hfa_bias;
  real sigma_game;
  {
    vector[G] resid;
    real bias_sum = 0;
    real n_home = 0;
    for (t in 1:n_team) {
      off_pts[t] = mu_p * off[t, n_season];
      def_pts[t] = mu_p * dfn[t, n_season];
      st_pts[t] = st[t, n_season];
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
    hfa_bias = n_home > 0 ? bias_sum / n_home : 0;
    sigma_game = sqrt(mean(square(resid - hca * hfa_bias)));
  }
}
"

# ==================================================================
# FIT
# ==================================================================
stan_fit <- stan(
  model_code = stan_code,
  data = stan_data,
  chains = n_chains,
  iter = n_iter,
  warmup = floor(n_iter / 2),
  seed = 2026,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

cat("\n=== SAMPLER DIAGNOSTICS ===\n")
rstan::check_hmc_diagnostics(stan_fit)
hyper <- c("mu", "h", "sigma_e", "tau_o", "tau_d", "tau_s", "rho_o", "rho_d", "rho_s",
           "sigma_st", "mu_p", "tau_p", "sigma_p", "hfa_bias", "sigma_game")
hs <- summary(stan_fit, pars = hyper)$summary
print(round(hs[, c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")], 3))
all_summary <- summary(stan_fit)$summary
cat("Max Rhat over all parameters:", round(max(all_summary[, "Rhat"], na.rm = TRUE), 3), "(want < 1.01)\n")
cat("Min effective sample size:", round(min(all_summary[, "n_eff"], na.rm = TRUE)), "\n")

ext <- rstan::extract(
  stan_fit,
  pars = c("off", "dfn", "st", "tempo", "mu_p", "h", "hfa_bias",
           "sigma_game", "strength_pts", "off_pts", "def_pts", "st_pts", "exp_margin")
)
n_draws <- length(ext$mu_p)

# ==================================================================
# MODEL CHECKS
# ==================================================================
em_mean <- colMeans(ext$exp_margin)
cal <- lm(gdat$margin ~ 0 + em_mean)
cat("\n=== CALIBRATION CHECK ===\n")
cat("Slope of actual margin on model-expected margin:", round(coef(cal), 3),
    "(1 = calibrated; in-sample, so mildly optimistic)\n")
cat("Median game-level SD (sigma_game):", round(median(ext$sigma_game), 2), "points\n")
cat("Median home margin missed by EPA model (hfa_bias):", round(median(ext$hfa_bias), 2), "points\n")
cat("Mean strength across teams:", round(mean(ext$strength_pts), 3), "points\n")

# ==================================================================
# RANKING WITH UNCERTAINTY
# ==================================================================
draw_summary <- function(M) {
  tibble(
    team = teams,
    med = apply(M, 2, median),
    lo50 = apply(M, 2, quantile, probs = 0.25),
    hi50 = apply(M, 2, quantile, probs = 0.75),
    lo90 = apply(M, 2, quantile, probs = 0.05),
    hi90 = apply(M, 2, quantile, probs = 0.95)
  )
}

ranks <- t(apply(-ext$strength_pts, 1, rank, ties.method = "min"))
cur <- max(seasons)
rec <- bind_rows(
  gdat %>% filter(season == cur) %>% transmute(team = home_team, w = margin > 0, l = margin < 0),
  gdat %>% filter(season == cur) %>% transmute(team = away_team, w = margin < 0, l = margin > 0)
) %>%
  group_by(team) %>%
  summarise(record = paste0(sum(w, na.rm = TRUE), "-", sum(l, na.rm = TRUE)), .groups = "drop")

ranking <- draw_summary(ext$strength_pts) %>%
  mutate(
    rank_med = apply(ranks, 2, median),
    rank_lo = apply(ranks, 2, quantile, probs = 0.05),
    rank_hi = apply(ranks, 2, quantile, probs = 0.95),
    p_top8 = colMeans(ranks <= 8),
    sd = apply(ext$strength_pts, 2, sd),
    off = apply(ext$off_pts, 2, median),
    def = apply(ext$def_pts, 2, median),
    st = apply(ext$st_pts, 2, median)
  ) %>%
  left_join(rec, by = "team") %>%
  arrange(desc(med)) %>%
  mutate(rank = row_number())

options(pillar.sigfig = 3)
cat("\n=== TEAM STRENGTH (points/game vs average team, neutral field; 90% intervals) ===\n")
print(
  ranking %>%
    dplyr::select(rank, team, record, strength = med, sd, lo90, hi90, rank_lo, rank_hi, p_top8, off, def, st),
  n = Inf, width = Inf
)

# 1) Strength rankings with 50% and 90% intervals
p_strength <- ggplot(slice(ranking, 1:25), aes(x = reorder(team, med), y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lo90, ymax = hi90), linewidth = 0.6, alpha = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.8) +
  geom_point(size = 2) +
  coord_flip() +
  labs(
    x = NULL, y = "Strength (points per game vs average team)",
    title = paste0("College football team strength, ", cur, " season"),
    subtitle = "Posterior median; thick = 50% interval, thin = 90% interval"
  ) +
  theme_minimal()
print(p_strength)

# 2) Rank intervals
p_rank <- ggplot(ranking, aes(x = reorder(team, -rank_med), y = rank_med)) +
  geom_linerange(aes(ymin = rank_lo, ymax = rank_hi), linewidth = 1) +
  geom_point(size = 2) +
  scale_y_reverse(breaks = seq(1, n_team, by = 5)) +
  coord_flip() +
  labs(
    x = NULL, y = "Rank (1 = best)",
    title = "Uncertainty in college football team rankings",
    subtitle = "Posterior median rank and 90% interval"
  ) +
  theme_minimal()
print(p_rank)

# 3) Components: offense, defense, and special teams
comp <- bind_rows(
  draw_summary(ext$off_pts) %>% mutate(component = "Offense"),
  draw_summary(ext$def_pts) %>% mutate(component = "Defense"),
  draw_summary(ext$st_pts) %>% mutate(component = "Special teams")
) %>%
  mutate(
    team = factor(team, levels = rev(ranking$team)),
    component = factor(component, levels = c("Offense", "Defense", "Special teams"))
  )

p_components <- ggplot(comp, aes(x = team, y = med)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = lo90, ymax = hi90), linewidth = 0.5, alpha = 0.6) +
  geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 1.4) +
  geom_point(size = 1.4) +
  coord_flip() +
  facet_wrap(~component) +
  labs(
    x = NULL, y = "Points per game vs average",
    title = "Strength components (teams ordered by total strength)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))
print(p_components)

# ==================================================================
# MATCHUP SIMULATOR
# ==================================================================
matchup <- function(home, away, neutral = FALSE, market_spread = NULL, level = 0.90) {
  stopifnot(home %in% teams, away %in% teams, home != away)
  ih <- match(home, teams)
  ia <- match(away, teams)
  hca <- as.numeric(!neutral)
  S <- n_season

  net <- function(i) ext$off[, i, S] + ext$dfn[, i, S]
  poss_d <- ext$mu_p + ext$tempo[, ih] + ext$tempo[, ia]

  c_ed <- poss_d * (net(ih) - net(ia))
  c_st <- ext$st[, ih, S] - ext$st[, ia, S]
  c_hfa <- (poss_d * 2 * ext$h + ext$hfa_bias) * hca
  line <- c_ed + c_st + c_hfa

  # Limit repeated game-noise simulations to keep college matchups quick.
  game <- rep(line, each = game_reps) +
    rnorm(length(line) * game_reps, 0, rep(ext$sigma_game, each = game_reps))

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  qi <- function(x) quantile(x, probs = probs, na.rm = TRUE)
  win_prob <- mean(pnorm(line / pmax(ext$sigma_game, 1e-6)))

  cat(sprintf("\n=== %s (home) vs %s (away)%s ===\n", home, away,
              if (neutral) " [neutral site]" else ""))
  cat(sprintf("Expected home margin: %+.1f | model spread for %s: %+.1f\n",
              mean(line), home, -mean(line)))
  cat(sprintf("Uncertainty in expected margin (parameters only): SD %.2f, %d%% interval [%+.1f, %+.1f]\n",
              sd(line), level * 100, qi(line)[1], qi(line)[2]))
  cat(sprintf("SD by source: offense/defense %.2f | special teams %.2f | home field %.2f\n",
              sd(c_ed), sd(c_st), sd(c_hfa)))
  cat(sprintf("Predicted game result (incl. game noise): SD %.1f, %d%% interval [%+.1f, %+.1f]\n",
              sd(game), level * 100, qi(game)[1], qi(game)[2]))
  cat(sprintf("Win probability: %s %.1f%% | %s %.1f%%\n",
              home, 100 * win_prob, away, 100 * (1 - win_prob)))
  cat(sprintf("P(%s is better on a neutral field): %.1f%%\n",
              home, 100 * mean(c_ed + c_st > 0)))

  if (!is.null(market_spread)) {
    cat(sprintf(
      "Market spread for %s %+.1f: P(%s covers) = %.1f%% | P(expected margin beats market) = %.1f%%\n",
      home, market_spread, home, 100 * mean(game > -market_spread),
      100 * mean(line > -market_spread)
    ))
  }

  dd <- bind_rows(
    tibble(value = line, what = "Expected margin (parameter uncertainty only)"),
    tibble(value = game, what = "Game result (parameters + game noise)")
  )
  g <- ggplot(dd, aes(x = value, fill = what)) +
    geom_density(alpha = 0.45, color = NA) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = mean(line), color = "black") +
    labs(
      x = paste0("Home margin (", home, " minus ", away, ")"), y = NULL, fill = NULL,
      title = paste0(home, " vs ", away, if (neutral) " (neutral site)" else ""),
      subtitle = sprintf("Expected margin %+.1f (spread %+.1f); win probability %.0f%%",
                         mean(line), -mean(line), 100 * win_prob)
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", axis.text.y = element_blank())
  if (!is.null(market_spread)) g <- g + geom_vline(xintercept = -market_spread, color = "red")
  print(g)

  invisible(list(line = line, game = game,
                 parts = tibble(offdef = c_ed, st = c_st, hfa = c_hfa)))
}

# ==================================================================
# CHOOSE A MATCHUP
# Team names must exactly match the names printed in `teams` / `ranking`.
# ==================================================================
home_team <- "Auburn"
away_team <- "LSU"
neutral_game <- FALSE
market_spread <- NULL       # e.g. -3.5 if home team is a 3.5-point favourite
use_prompt <- FALSE         # TRUE asks for teams in an interactive console

if (use_prompt && interactive()) {
  cat("Teams:\n", paste(teams, collapse = ", "), "\n")
  home_team <- readline("Home team: ")
  away_team <- readline("Away team: ")
  neutral_game <- tolower(readline("Neutral site? (y/n): ")) == "y"
  ms <- readline("Market spread for the home team (blank = none): ")
  market_spread <- if (nzchar(ms)) as.numeric(ms) else NULL
}

res <- matchup(home_team, away_team, neutral = neutral_game,
               market_spread = market_spread)
