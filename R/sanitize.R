clip_outliers <- function(Y, thresh=4) {
  
  tY <- t(Y)
  med.y <- apply(Y, 2, median)
  mad.y <- apply(Y, 2, mad)
  
  outlier_rows <- function(Y, med.y, mad.y) {
    rowSums(abs((Y - med.y)/mad.y) > thresh + 1e-4) >= 1
  }
  
  clip_row <- function(Y, med.y, mad.y) {
    Y <- pmin(Y, med.y + thresh*mad.y)
    Y <- pmax(Y, med.y - thresh*mad.y)
    Y
  }
  
  outliers <- outlier_rows(tY, med.y, mad.y)
  
  while(any(outliers)) {
    tY[outliers, ] <- clip_row(tY[outliers, ,drop=FALSE], 
                               med.y[outliers], 
                               mad.y[outliers])
    mad.y[outliers] <- apply(tY[outliers, ,drop=FALSE], 1, mad)
    med.y[outliers] <- apply(tY[outliers, ,drop=FALSE], 1, median)
    outliers[outliers] <- outlier_rows(
      tY[outliers, ,drop=FALSE],
      med.y[outliers],
      mad.y[outliers]
    ) 
  }
  
  t(tY)
}