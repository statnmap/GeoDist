usethis::use_build_ignore("dev_history.R")

# Doc
usethis::use_readme_rmd()
usethis::use_vignette("name")
usethis::use_code_of_conduct()

# Dependencies
attachment::att_amend_desc()

# CI
# remotes::install_github("ropensci/tic")
# tic::use_tic()
# tic::use_tic(wizard = FALSE, linux = "travis", mac = "none", windows = "appveyor", deploy = "travis",
#              matrix = "none",
#              travis_endpoint = ".org")
# travis::travis_enable(endpoint = ".org")
# travis::browse_travis_token(endpoint = '.org')


