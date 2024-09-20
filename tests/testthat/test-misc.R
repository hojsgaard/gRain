context("gRain misc")
##############################

yn   <- c("yes", "no")
a    <- cpt(~asia,                  values=c(1,99), levels=yn)
t.a  <- cpt(~tub + asia,            values=c(5,95,1,99), levels=yn)
s    <- cpt(~smoke,                 values=c(5,5), levels=yn)
l.s  <- cpt(~lung + smoke,          values=c(1,9,1,99), levels=yn)
b.s  <- cpt(~bronc + smoke,         values=c(6,4,3,7), levels=yn)
e.lt <- cpt(~either + lung + tub,   values=c(1,0,1,0,1,0,0,1), levels=yn)
x.e  <- cpt(~xray + either,         values=c(98,2,5,95), levels=yn)
d.be <- cpt(~dysp + bronc + either, values=c(9,1,7,3,8,2,1,9), levels=yn)
cpt_list  <- list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be)
chest_cpt <- compile_cpt(cpt_list)
bn <- grain(chest_cpt)

disease <- c("tub", "lung", "bronc")
evidence <- list(asia="yes", dysp="yes")

test_that("nodeNames()", {
    vn <- c("asia", "tub", "smoke", "lung", "bronc", "either", "xray", "dysp")
    expect_setequal(nodeNames(bn), vn)

    pr <- querygrain(bn, nodes=disease)
    bn2 <- evidence_add(bn, evidence)
    
    expect_equal(as.numeric(pr$tub)[1], 0.0104)
    expect_equal(as.numeric(pr$lung)[1], 0.055)
    expect_equal(as.numeric(pr$bronc)[1], 0.45)

    expect_null(evidence_prob(bn))
    expect_equal(evidence_prob(bn2), 0.004501375)
    
})

