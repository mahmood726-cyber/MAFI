# Create example dataset for MAFI package
# Uses classic BCG vaccine meta-analysis from metafor

library(metafor)

# BCG vaccine effectiveness
data(dat.bcg)
bcg <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)

mafi_example_bcg <- data.frame(
  study = dat.bcg$author,
  year = dat.bcg$year,
  yi = as.numeric(bcg$yi),
  vi = as.numeric(bcg$vi),
  sei = sqrt(as.numeric(bcg$vi))
)

# Save
save(mafi_example_bcg, file = "data/mafi_example_bcg.rda", compress = "xz")

cat("Created mafi_example_bcg with", nrow(mafi_example_bcg), "studies\n")
