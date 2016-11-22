ASReml's Astar for fixed genetic groups implicitly within the random effects
 id !P
 p
 is !I 0 1
 gen 15
 f
 foc0
 g1
asremlGG.ped !ALPHA !SKIP 1 !MAKE !GIV !DIAG !GROUPS 2
ggTutorial.dat !SKIP 1
!LAST id 2
p ~ mu f !r id
