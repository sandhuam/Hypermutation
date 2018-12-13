# 6 sided die, 1 comes up 11 times, odds that the die is fair (percent)
fairdie = integer()
fairodds = numeric()
avgchance = numeric()
for (k in 1:100){
for (j in 1:100){
  fairdie = replicate(100, sample(1:6, 1, T))
if(length(which(fairdie == 1)) == 11){
  fairodds[j] = 1
} else {
    fairodds[j] = 0
  }
}
fairchance = length(which(fairodds == 1))
avgchance[k] = fairchance
}
mean(avgchance)
# How many 1's to say die is unfair with 95% certainty
list1 = numeric()
for (i in 1:10000){
fairdie2 = replicate(100, sample(1:6, 1, T))
list1[i] = length(which(fairdie2 == 1))
}
d = head(sort(list1), length(list1)*0.05)
if(max(d) %in% list1[-d]) {
  max(d) - 1
} else {
  max(d)
}
# 1 came up 11 times but die is weighted {1:0.3,2:0.3,3:0.1,4:0.1,5:0.1,6:0.1}
unfairdie = replicate(100, sample(1:6, 1, T, prob = as.numeric(c("0.3", "0.3", "0.1", "0.1", "0.1", "0.1"))))
if (length(unfairdie[unfairdie == 1]) == 11) {
  print("consistent")
} else if (length(unfairdie[unfairdie == 1]) %in% 8:14) {
  print("relatively consistent")
} else {
  print("inconsistent")
}
# True odds of rolling a "1"; should be 30% based on the code above
unfairdie2 = replicate(100, replicate(100, sample(1:6, 1, T, prob = as.numeric(c("0.3", "0.3", "0.1", "0.1", "0.1", "0.1")))))
length(unfairdie2[unfairdie2 == 1])/10000
                       