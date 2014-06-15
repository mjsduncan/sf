### regress expression value on age

fits <- function (dat) {
  lm.fit(cbind(1, seq_len(nrow(dat))), t(dat))
  dat <- cbind(dat, t(coef(fits)))
  names(dat)[-(1:3)] <- c("Intercept","Slope")
  t(coef(fits))
}