
sp <- qread("./simulation/tmp.qs")
sp2 <- qread("./simulation/tmp.qs")
setDT(sp$pop)
setDT(sp2$pop)
all.equal(sp$pop, sp2$pop)
l <- mk_scenario_init2("", diseases, sp, design)
l2 <- mk_scenario_init2("", diseases, sp2, design)
all.equal(l, l2)
simcpp(sp$pop, l, sp$mc_aggr)
simcpp(sp2$pop, l2, sp2$mc_aggr)
all.equal(sp$pop, sp2$pop)

i = 1L
while (identical(sp$pop[seq(i)], sp2$pop[seq(i)])) {
  i <- i + 1e4
  print(i)
}
i <- i - 1e4
while (identical(sp$pop[seq(i)], sp2$pop[seq(i)])) {
i <- i + 1e3
print(i)
}
i <- i - 1e3
while (identical(sp$pop[seq(i)], sp2$pop[seq(i)])) {
  i <- i + 1e2
  print(i)
}
i <- i - 1e2
while (identical(sp$pop[i], sp2$pop[i])) {
  i <- i + 1
  print(i)
}
all.equal(sp$pop[i], sp2$pop[i])
for (j in seq_along(sp$pop)) {
  if (!identical(sp$pop[i, j, with = F], sp2$pop[i, j, with = F]))
    print(rbind(sp$pop[i, j, with = F], sp2$pop[i, j, with = F]))
}
tt <- sp$pop[i, pid]
sp$pop[pid == tt, .(breast_ca_prvl, all_cause_mrtl)]
sp2$pop[pid == tt, .(breast_ca_prvl, all_cause_mrtl)]
identical(sp$pop[pid == tt, t2dm_prvl], sp2$pop[pid == tt, t2dm_prvl])
identical(sp$pop[pid == tt, prb_breast_ca_mrtl1], sp2$pop[pid == tt, prb_breast_ca_mrtl1])
identical(sp$pop[pid == tt, prb_breast_ca_mrtl2], sp2$pop[pid == tt, prb_breast_ca_mrtl2])

all.equal(sp$pop[i], sp2$pop[i])
View(rbind(sp$pop[i], sp2$pop[i]))
