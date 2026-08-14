library(tern, warn.conflicts = FALSE)
library(chevron)
library(citril)
library(haven)
library(dplyr, warn.conflicts = FALSE)

# setup ----
molecule <- Sys.getenv("AEGIS_MOLECULE")
study <- Sys.getenv("AEGIS_STUDY")
rep_event <- unlist(strsplit(getwd(), .Platform$file.sep))[4]
run_area <- ifelse(interactive(), "dev", Sys.getenv("AEGIS_ENV"))
rep_path <- file.path("", molecule, study, rep_event)

args <- commandArgs(trailingOnly = TRUE)

## prepare env ----
paths <- list(
  wd = getwd(),
  meddra = "/aegis/references/MedDRA",
  adam = file.path(rep_path, run_area, "data/adam"),
  metadataout = file.path(rep_path, run_area, "data/metadataout"),
  entimice_adam = file.path(rep_path, "ref/entimice/adam"),
  repout = file.path(rep_path, run_area, "repout"),
  output = file.path(rep_path, run_area, "output")
)

## filters ----
filters <- args[1]
if (interactive()) {
  filters <- "SE"
}
filter_labels <- unlist(strsplit(filters, "_"))

## data input ----
options(
  citril.path.adam = paths$entimice_adam,
  citril.path.reporting = file.path(paths$wd, "activity.json"),
  citril.path.filters = file.path(paths$wd, "filters.json"),
  citril.path.formats = file.path(paths$wd, "meta_formats.yaml"),
  citril.path.rules = file.path(paths$wd, "meta_rules.yaml"),
  citril.path.whiskers = file.path(paths$wd, "meta_whiskers.yaml"),
  citril.path.specs = file.path(paths$wd, "arg_specs.yaml"),
  citril.path.rds = paths$repout,
  citril.path.att = paths$metadataout
)

adam_db <- read_adam(sas7bdat = read_sas, select = "adsl")

reporting_info <- read_reporting_info()
filters_all <- read_filters(default_filters = list())
formats_all <- read_formats()
rules_all <- read_rules()
read_whiskers()

prog_name <- "t_dm"
output <- file.path(paths$output, paste0(prog_name, "_", filters, ".out"))

# 削除
unlink(output, force = TRUE)

tlg_name <- parse_tlg_name(prog_name, filters)

arg_spec <- read_specs(entry = tlg_name)
tlg_spec <- list(
  lopo_titles = "Demographics"
)

# filter data
adam_db <- filter_with_suffix(
  adam_db,
  filters,
  suffix_list = filters_all
)

# reformat data
adam_db <- extract_rules(
  format_list = formats_all[[tlg_name]],
  rules_list = rules_all,
  default = formats_all[["all"]]
) |>
  reformat(obj = adam_db, format = _)


# craete table ----

adsl <- adam_db$adsl |>
  mutate(
    RACE = as.factor(RACE)
  )

lyt <- basic_table(show_colcounts = TRUE) |>
  split_cols_by(var = "ACTARM") |>
  add_overall_col("All Patients") |>
  analyze_vars(
    vars = c("AGE", "AGEGR1", "SEX", "RACE"),
    var_labels = c("Age (yr)", "Age Group", "Sex", "Race")
  )

tlg_output <- build_table(lyt, adsl)

final_tlg <- decorate_tlg(
  tlg_output,
  main = append_filter_title(
    main = tlg_spec$lopo_titles,
    suffix = filter_labels,
    filter_list = filters_all
  ),
  sub = reporting_info$title,
  footnote = c(tlg_spec$lopo_footnotes, reporting_info$footnote),
  sup_footnote = run_information(output),
  pagesize = "L7"
)

export_tlg(final_tlg, file = output)

output_file <- parse_output_path(output)

unlink(output_file$att, force = TRUE)

make_att(final_tlg, output_file$att)
