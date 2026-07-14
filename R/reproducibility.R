# ============================================================
# RNA-Seq Project Setup Script (Non-interactive, GitHub-enabled)
# ============================================================

# ---- 0. Load or install required packages ----
required_pkgs <- c("usethis", "renv", "gitcreds")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(usethis)
library(renv)
library(gitcreds)

# ---- 1. Create new R project ----
#project_name <- "rnaseq_ProjectName4"
project_path <- "/home/rstudio/project"

setwd(project_path)

usethis::create_project(path=".", open = FALSE)


# ---- 2. Initialize renv ----
renv::init(bare = TRUE, bioconductor = TRUE, restart = FALSE)
#renv::snapshot(prompt = FALSE)
# ---- 3. Install essential packages ----

renv::install(packages = read.table("R/lib.list")$V1, repos = c(BiocManager::repositories(), getOption("repos")), prompt = FALSE)

# ---- 4. Snapshot renv ----
renv::snapshot(prompt = FALSE)

# ---- 5. Create directory structure ----
dirs <- c(
  "report",
  "out",
  "scripts",
  "logs"
)
lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# ---- 6. RNA-Seq specific .gitignore ----
gitignore_lines <- c(
  "# RNA-Seq specific ignores",
  "data/raw_data/*",
  "data/reference_data/*",
  "data/meta_data/*",
  "data/processed_data/*",
  "results/*",
  "logs/",
  ".Rhistory",
  ".RData",
  ".Rproj.user",
  "renv/library/",
  "renv/staging/",
  ".Renviron"
)
writeLines(gitignore_lines, ".gitignore")

# ---- 7. Create README ----
readme_content <- c(
  "# RNA-Seq Analysis Project",
  "",
  "This repository contains a reproducible RNA-Seq analysis workflow.",
  "",
  "## Directory structure",
  "- `data/raw_data/` — Raw sequencing files (FASTQ, etc.)",
  "- `data/reference_data/` — Reference genome, annotation files",
  "- `data/meta_data/` — Sample metadata",
  "- `data/processed_data/` — Processed data: trimming, alignment, counts",
  "- `results/` — QC, differential expression, and functional profiling results",
  "- `report/` — Analysis reports",
  "- `scripts/` — R scripts for analysis steps",
  "- `R/` — Internal functions",
  "- `logs/` — Log files from pipelines",
  "",
  "## Environment",
  "This project uses **renv** for reproducibility.",
  "Run `renv::restore()` to install the correct package versions.",
  "",
  "## Version control",
  "Git and GitHub are configured automatically by this setup script."
)
writeLines(readme_content, "README.md")

# ---- 8. Setup GitHub ----
# You must have a GitHub Personal Access Token (PAT) stored securely.
# If not already configured, run this once manually before using the script:
gitcreds::gitcreds_set()

# ---- 9. Initialize Git ----
usethis::use_git(message = "Initial commit: setup RNA-Seq project")

# # Retrieve token from the gitcreds store
# token <- gitcreds::gitcreds_get()$password
#
# # Set environment variable for GitHub authentication
# Sys.setenv(GITHUB_PAT = token)

# # Create and push to GitHub repo automatically (non-interactive)
# usethis::use_github(private = FALSE)

# ---- 10. Done ----
cat("✅ RNA-Seq project fully initialized and pushed to GitHub at:", project_path, "\n")


rstudioapi::openProject(project_path)
