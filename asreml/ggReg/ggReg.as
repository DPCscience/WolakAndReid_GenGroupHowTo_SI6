Fixed regression for Explicit genetic groups
 id !P
 p
 is !I 0 1
 gen 15
 f
 foc0
 g1
reg.ped !ALPHA !SKIP 1 !MAKE !GIV !DIAG
ggTutorial.dat !SKIP 1
p ~ mu f g1 foc0 !r id
