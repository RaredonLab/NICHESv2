require(roxygen2)
require(devtools)
require(pkgdown)

# Document
setwd("~/GitHub/NICHESv2-msbr")
document()
check()

# Check installation from local
setwd("~/GitHub")
install('NICHESv2-msbr')

# Update website
setwd("~/GitHub/NICHESv2-msbr")
pkgdown::build_site(lazy = T)

# pkgdown::build_article(name = "07 Spatiotemporal NICHES",lazy = T)

# Install from release
# install_github('RaredonLab/NICHESv2',
#                ref = 'v2.0.0',
#                auth_token = 'f7ea5d8790fe721ac0c9d5ef115d04068b19ed6d',
#                force = T)
