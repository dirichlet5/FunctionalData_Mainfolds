evaluation = function(n, results, true){
  tp = length(intersect(results, true))/length(true)
  fp = length(setdiff(results, true))/(n - length(true))
  list(tp = tp, fp = fp)
}