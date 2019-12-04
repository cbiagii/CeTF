bivar.awk <- function(df) {
  cat(green(
    "########## 1th Column ##########" %+% '\n' %+%
      sprintf("N     =  %s", length(df[,1])) %+% '\n' %+%
      sprintf("Mean  =  %s", format(round(mean(df[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Std   =  %s", format(round(sd(df[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Min   =  %s", format(round(min(df[,1]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Max   =  %s", format(round(max(df[,1]), 4), nsmall = 4)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  cat(red(
    "########## 2th Column ##########" %+% '\n' %+%
      sprintf("N     =  %s", length(df[,2])) %+% '\n' %+%
      sprintf("Mean  =  %s", format(round(mean(df[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Std   =  %s", format(round(sd(df[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Min   =  %s", format(round(min(df[,2]), 4), nsmall = 4)) %+% '\n' %+%
      sprintf("Max   =  %s", format(round(max(df[,2]), 4), nsmall = 4)) %+% '\n' %+%
      "################################" %+% '\n'
  ))

  cat(blue(
    sprintf("################################" %+% '\n' %+%
              "Correlation  =  %s", format(round(cor(df[,1], df[,2]), 7), nsmall = 7)) %+% '\n' %+%
      sprintf("Regression   =  %s", format(round(lm(df[,2] ~ df[,1], data = df)[[1]][[2]], 7), nsmall = 7)) %+% '\n' %+%
      "################################" %+% '\n'
  ))
}
