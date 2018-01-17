
PH_MIN = 0  # minimum pH value
PH_MAX = 14  # maximum pH value
MAXLOOP = 2000  # maximum number of iterations
EPSI = 0.0001  # desired precision

pk = list(A = c(3.55, 7.59, 0.0  , 0.0  , 0.0),
B = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
C = c(3.55, 7.50, 9.00 , 9.00 , 9.00),
D = c(4.55, 7.50, 4.05 , 4.05 , 4.05),
E = c(4.75, 7.70, 4.45 , 4.45 , 4.45),
F = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
G = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
H = c(3.55, 7.50, 5.98 , 5.98 , 5.98),
I = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
J = c(0.00, 0.00, 0.0  , 0.0  , 0.0),
K = c(3.55, 7.50, 10.00, 10.00, 10.00),
L = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
M = c(3.55, 7.00, 0.0  , 0.0  , 0.0),
N = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
O = c(0.00, 0.00, 0.0  , 0.0  , 0.0),
P = c(3.55, 8.36, 0.0  , 0.0  , 0.0),
Q = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
R = c(3.55, 7.50, 12.0 , 12.0 , 12.0),
S = c(3.55, 6.93, 0.0  , 0.0  , 0.0),
T = c(3.55, 6.82, 0.0  , 0.0  , 0.0),
U = c(0.00, 0.00, 0.0  , 0.0  , 0.0),
V = c(3.55, 7.44, 0.0  , 0.0  , 0.0),
W = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
X = c(3.55, 7.50, 0.0  , 0.0  , 0.0),
Y = c(3.55, 7.50, 10.00, 10.00, 10.00),
Z = c(3.55, 7.50, 0.0  , 0.0  , 0.0) )

hydrophobicityIndex = c(A = 1.8,
  R = -4.5,
  N = -3.5,
  D=-3.5,
  C=2.5,
  E=-3.5,
  Q=-3.5,
  G=-0.4,
  H=-3.2,
  I=4.5,
  L=3.8,
  K=-3.9,
  M=1.9,
  F=2.8,
  P=-1.6,
  S=-0.8,
  T=-0.7,
  W=-0.9,
  Y=-1.3,
  V=4.2)

pI <- function(peptide = AAStringSet("ONLYCAPITALLETTERSXLE"), charge_increment = 0.0) {
  res <- rep(0.0, length(peptide))
  composition = t(sapply(strsplit(as.character(peptide), split=""), function(x) { table(factor(x,levels=names(pk))) } ))
  for (k in 1:nrow(composition)) {
#     cat("Calculate pI for peptide ",k," / ",nrow(composition),"\r")
    comp = composition[k,]
    nterm_res = substr(peptide[k], 1, 1)
    cterm_res = substr(peptide[k], width(peptide[k]), width(peptide[k]))
    ph_min = PH_MIN
    ph_max = PH_MAX

    i = 1
    while ((i < MAXLOOP) & ((ph_max - ph_min)>EPSI)) {
      ph_mid = ph_min + (ph_max - ph_min) / 2.0

      cter = 10^(-pk[[cterm_res]][1]) / (10^(-pk[[cterm_res]][1]) + 10^(-ph_mid));
      nter = 10^(-ph_mid) / (10^(-pk[[nterm_res]][2]) + 10^(-ph_mid));

      carg = comp[["R"]] * 10^(-ph_mid) / (10^(-pk[["R"]][3]) + 10^(-ph_mid));
      chis = comp[["H"]] * 10^(-ph_mid) / (10^(-pk[["H"]][3]) + 10^(-ph_mid));
      clys = comp[["K"]] * 10^(-ph_mid) / (10^(-pk[["K"]][3]) + 10^(-ph_mid));

      casp = comp[["D"]] * 10^(-pk[["D"]][3]) / (10^(-pk[["D"]][3]) + 10^(-ph_mid));
      cglu = comp[["E"]] * 10^(-pk[["E"]][3]) / (10^(-pk[["E"]][3]) + 10^(-ph_mid));

      ccys = comp[["C"]] * 10^(-pk[["C"]][3]) / (10^(-pk[["C"]][3]) + 10^(-ph_mid));
      ctyr = comp[["Y"]] * 10^(-pk[["Y"]][3]) / (10^(-pk[["Y"]][3]) + 10^(-ph_mid));

      charge = carg + clys + chis + nter + charge_increment - (casp + cglu + ctyr + ccys + cter);
      if (charge > 0.0) {
        ph_min = ph_mid;
      } else {
        ph_max = ph_mid;
      }   
      i = i+1
    }
    res[k] = ph_mid
  }

  return(res);
}

hydrophobicity = function(peptide) {
  sapply(strsplit(as.character(peptide), split=""), function(x) {
    t = table(factor(x,levels=names(hydrophobicityIndex)))[]
    sum(hydrophobicityIndex * t) / sum(t)
  } )
}
