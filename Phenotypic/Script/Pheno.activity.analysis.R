setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls())
if (dev.cur() > 1) dev.off()
gc()

library(data.table)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

## 1. Settings ----

cold_col <- "steelblue"
hot_col  <- "maroon"

theme_activity <- theme_classic(base_size = 22) +
  theme(
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    plot.margin = margin(10, 15, 15, 10)
  )


## 2. Helper functions ----

parse_meta_datetime <- function(x, end_of_day = FALSE) {
  x <- trimws(as.character(x))
  out <- as.POSIXct(rep(NA_character_, length(x)))
  has_time <- !is.na(x) & x != "" & grepl(":", x)
  date_only <- !is.na(x) & x != "" & !has_time
  
  out[has_time] <- as.POSIXct(x[has_time], format = "%m/%d/%y %H:%M")
  out[date_only] <- as.POSIXct(x[date_only], format = "%m/%d/%y")
  if (end_of_day) out[date_only] <- out[date_only] + 24 * 60 * 60 - 1
  
  out
}

read_monitor_file <- function(path) {
  x <- fread(path, header = FALSE)
  chamber_cols <- paste0("chamber_", seq_len(ncol(x) - 10))
  
  setnames(x, c("record_id", "date_raw", "time_raw", "col4", "col5",
                "monitor_id", "col7", "status", "col9", "light", chamber_cols))
  
  x[, file := basename(path)]
  x[, datetime := as.POSIXct(paste(date_raw, time_raw), format = "%d %b %y %H:%M:%S")]
  x[]
}


## 3. Read metadata ----

meta <- fread("../Data/activity.metadata.csv", na.strings = c("", "NA", "NaN"))
setnames(meta, sub("^\\ufeff", "", names(meta)))

meta[, file := basename(file)]
meta[, region_id := as.integer(region_id)]
meta[, replicate := as.integer(replicate)]
meta[, survivors := as.integer(survivors)]
meta[, start_datetime := parse_meta_datetime(start_datetime)]
meta[, stop_datetime  := parse_meta_datetime(stop_datetime, end_of_day = TRUE)]

meta[, evo := fifelse(grepl("^HOT", genotype, ignore.case = TRUE), "Hot",
                      fifelse(grepl("^COLD", genotype, ignore.case = TRUE), "Cold", NA_character_))]
meta[, genotype_id := as.integer(gsub("\\D+", "", genotype))]
meta[, sample_id := fifelse(is.na(genotype), NA_character_, paste(genotype, replicate, sep = "_rep"))]

meta_used <- meta[!is.na(genotype)]


## 4. Read monitor files ----

monitor_paths <- file.path("../Data", unique(meta$file))
monitor_wide <- rbindlist(lapply(monitor_paths, read_monitor_file), use.names = TRUE, fill = TRUE)
monitor_wide <- monitor_wide[, -c(4, 5, 7, 8, 9)]


## 5. Make long activity table ----

chamber_cols <- grep("^chamber_", names(monitor_wide), value = TRUE)

activity_long <- melt(
  monitor_wide,
  id.vars = setdiff(names(monitor_wide), chamber_cols),
  measure.vars = chamber_cols,
  variable.name = "region_id",
  value.name = "activity_count"
)

activity_long[, region_id := as.integer(sub("^chamber_", "", region_id))]
activity_long[, activity_count := as.integer(activity_count)]

activity_long <- merge(activity_long, meta_used, by = c("file", "region_id"), all = FALSE)
activity_long <- activity_long[datetime >= start_datetime & datetime <= stop_datetime]

activity_long[, time_from_start_h := as.numeric(difftime(datetime, start_datetime, units = "hours"))]
activity_long[, time_from_start_d := time_from_start_h / 24]
activity_long[, activity_per_survivor := activity_count / survivors]
activity_long[, evo := factor(evo, levels = c("Cold", "Hot"))]

setorder(activity_long, file, region_id, datetime)

fwrite(activity_long, "./activity_long_with_metadata.csv")


## 6. Sanity checks ----

cat("Metadata rows:", nrow(meta), "\n")
cat("Assigned chambers:", nrow(meta_used), "\n")
cat("Wide monitor rows:", nrow(monitor_wide), "\n")
cat("Long curated rows:", nrow(activity_long), "\n\n")

print(activity_long[, .(
  n_genotypes = uniqueN(genotype),
  n_chambers = uniqueN(sample_id),
  n_observations = .N
), by = evo])

print(activity_long[, .(
  n_chambers = uniqueN(region_id),
  n_observations = .N,
  first_time = min(datetime),
  last_time = max(datetime)
), by = file])


## 7. Make 1-hour bins ----

activity_long <- activity_long[survivors >= 6]

activity_long[, sample_time := datetime - 10 * 60]
activity_long[, bin_start := as.POSIXct(format(sample_time, "%Y-%m-%d %H:00:00"))]
activity_long[, clock_hour := as.integer(format(bin_start, "%H"))]

first_keep <- activity_long[light == 1 & clock_hour == 8, .(first_keep = min(bin_start)), by = file]
activity_long <- merge(activity_long, first_keep, by = "file")
activity_long <- activity_long[bin_start >= first_keep]

activity_long[, first_bin := min(bin_start), by = file]
activity_long[, time_1h := as.numeric(difftime(bin_start, first_bin, units = "hours"))]
activity_long[, assay_day := floor(time_1h / 24) + 1]
activity_long <- activity_long[time_from_start_h <= 247, ]
activity_long[, phase := fcase(
  clock_hour %in% c(8, 20), "shift",
  clock_hour >= 9 & clock_hour < 20, "light",
  default = "dark"
)]

activity_long[, phase := factor(phase, levels = c("dark", "shift", "light"))]

activity_1h <- activity_long[, .(
  activity_per_h = sum(activity_per_survivor),
  n_records = .N,
  light_prop = mean(light)
), by = .(evo, genotype, genotype_id, sample_id, replicate, time_1h, assay_day, phase)]

print(activity_1h[, .(
  n_rows = .N,
  n_samples = uniqueN(sample_id),
  min_records = min(n_records),
  max_records = max(n_records)
), by = .(assay_day, phase)])


## 8. Time-series summaries ----

activity_1h_genotype <- activity_1h[, .(
  activity_per_h = mean(activity_per_h)
), by = .(evo, genotype, genotype_id, time_1h, assay_day, phase)]

activity_1h_evo <- activity_1h_genotype[, .(
  mean_activity = mean(activity_per_h),
  se_activity = sd(activity_per_h) / sqrt(.N),
  assay_day = first(assay_day),
  phase = first(phase)
), by = .(evo, time_1h)]

p_genotype <- ggplot(activity_1h_genotype, aes(time_1h, activity_per_h, color = evo, group = genotype)) +
  geom_line(alpha = 0.45, linewidth = 0.8) +
  scale_color_manual(values = c("Cold" = cold_col, "Hot" = hot_col)) +
  labs(x = "Hours since first analyzed hour",
       y = "Activity per survivor per hour",
       color = "Regime") +
  theme_activity

ggsave("../Plot/activity_genotype_time_series.png", p_genotype, width = 14, height = 7)


## 9. Day-phase phenotype ----

day_phase_sample <- activity_1h[phase != "shift", .(
  activity_per_h = mean(activity_per_h),
  n_hours = .N
), by = .(evo, genotype, sample_id, assay_day, phase)]

print(day_phase_sample[, .(
  n_rows = .N,
  min_hours = min(n_hours),
  max_hours = max(n_hours)
), by = .(assay_day, phase)])

day_phase_sample <- day_phase_sample[n_hours >= 6]
day_phase_sample <- day_phase_sample[, if (uniqueN(evo) == 2) .SD, by = .(assay_day, phase)]

day_phase_sample[, assay_day := factor(assay_day)]
day_phase_sample[, phase := factor(as.character(phase), levels = c("dark", "light"))]

day_phase_genotype <- day_phase_sample[, .(
  activity_per_h = mean(activity_per_h)
), by = .(evo, genotype, assay_day, phase)]

day_phase_summary <- day_phase_genotype[, .(
  mean_activity = mean(activity_per_h),
  se_activity = sd(activity_per_h) / sqrt(.N)
), by = .(evo, assay_day, phase)]

p_day_phase <- ggplot(day_phase_summary, aes(assay_day, mean_activity, fill = evo)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = mean_activity - se_activity,
                    ymax = mean_activity + se_activity),
                position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.7) +
  facet_wrap(~phase) +
  scale_fill_manual(values = c("Cold" = cold_col, "Hot" = hot_col)) +
  labs(x = "Assay day",
       y = "Mean activity per survivor per hour",
       fill = "Regime") +
  theme_activity +
  theme(panel.spacing = unit(1.2, "lines"))

ggsave("../Plot/activity_day_phase_barplot.png", p_day_phase, width = 12, height = 7)


## 10. Day-phase model ----

m_day_phase <- lmer(log1p(activity_per_h) ~ evo * assay_day * phase +
                      (1 | genotype) + (1 | sample_id),
                    data = day_phase_sample)

anova(m_day_phase)

day_phase_test <- as.data.table(summary(
  contrast(emmeans(m_day_phase, ~ evo | assay_day + phase), method = "revpairwise")
))

day_phase_test[, q := p.adjust(p.value, method = "BH")]
day_phase_test[, sig := q < 0.05]
day_phase_test[, estimate_hot_minus_cold := estimate]

print(day_phase_test)




## 11. Significant pure dark/light windows ----

period_windows <- unique(activity_1h[, .(time_1h, assay_day, phase)])
setorder(period_windows, time_1h)

period_windows[, new_block := phase != shift(phase, fill = first(phase))]
period_windows[, block := cumsum(new_block)]

period_windows <- period_windows[, .(
  xmin = min(time_1h),
  xmax = max(time_1h) + 1,
  assay_day = first(assay_day),
  phase = first(as.character(phase))
), by = block]

sig_windows <- day_phase_test[sig == TRUE, .(
  assay_day = as.integer(as.character(assay_day)),
  phase = as.character(phase),
  estimate_hot_minus_cold,
  p.value,
  q
)]

sig_windows <- merge(sig_windows, period_windows, by = c("assay_day", "phase"), all.x = TRUE)
setorder(sig_windows, assay_day, phase)

print(sig_windows)




## 12. Mean curve with phase bar ----
y_top <- max(activity_1h_evo$mean_activity + activity_1h_evo$se_activity, na.rm = TRUE)

if (nrow(sig_windows) > 0) {
  sig_windows[, sig_label := fcase(
    q < 0.0001, "****",
    q < 0.001,  "***",
    q < 0.01,   "**",
    q < 0.05,   "*",
    default = ""
  )]
} else {
  sig_windows[, sig_label := character()]
}

p_mean_curve <- ggplot(activity_1h_evo, aes(time_1h, mean_activity)) +
  geom_rect(data = sig_windows,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.14) +
  geom_rect(data = period_windows,
            aes(xmin = xmin, xmax = xmax,
                ymin = 1.015 * y_top, ymax = 1.105 * y_top,
                fill = phase),
            inherit.aes = FALSE, color = NA) +
  geom_text(data = sig_windows,
            aes(x = (xmin + xmax) / 2, y = 1.055 * y_top, label = sig_label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +
  geom_ribbon(data = activity_1h_evo[evo == "Cold"],
              aes(x = time_1h,
                  ymin = pmax(0, mean_activity - se_activity),
                  ymax = mean_activity + se_activity),
              inherit.aes = FALSE, fill = cold_col, alpha = 0.15) +
  geom_ribbon(data = activity_1h_evo[evo == "Hot"],
              aes(x = time_1h,
                  ymin = pmax(0, mean_activity - se_activity),
                  ymax = mean_activity + se_activity),
              inherit.aes = FALSE, fill = hot_col, alpha = 0.15) +
  geom_line(aes(color = evo), linewidth = 1.1) +
  scale_color_manual(values = c("Cold" = cold_col, "Hot" = hot_col)) +
  scale_fill_manual(values = c("light" = "gold", "dark" = "grey65", "shift" = "orange")) +
  coord_cartesian(ylim = c(0, 1.10 * y_top), clip = "off") +
  labs(
    x = "Time (hr)",
    y = "Activity / fly / hour",
    caption = "**** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05"
  ) +
  theme_activity +
  theme(
    legend.position = "none",
    plot.caption = element_text(size = 16, hjust = 0.5),
    plot.margin = margin(12, 15, 25, 10)
  )

png(filename = "../Plot/activity_mean_time_series_with_phase_bar.png", width = 18, height =8, unit = "in", res = 450)
print(p_mean_curve)
dev.off()

saveRDS(p_mean_curve, "./activity.curve.plot.RDS")
saveRDS(list(
  activity_1h_evo = activity_1h_evo,
  period_windows = period_windows,
  sig_windows = sig_windows
), "./activity.curve.data.RDS")
## 13. Day-phase p-values ----

p_pvalue <- ggplot(day_phase_test, aes(as.numeric(as.character(assay_day)), -log10(q), color = phase)) +
  geom_point(size = 3.2, alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.8) +
  scale_color_manual(values = c("dark" = "grey40", "light" = "goldenrod3")) +
  scale_x_continuous(breaks = sort(unique(as.numeric(as.character(day_phase_test$assay_day))))) +
  labs(x = "Assay day",
       y = "-log10(FDR-adjusted p-value)",
       color = "Phase") +
  theme_activity

ggsave("../Plot/activity_day_phase_pvalues.png", p_pvalue, width = 10, height = 7)