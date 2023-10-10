(TeX-add-style-hook "manual-BOX-ggH"
 (lambda ()
    (LaTeX-add-environments
     "enumeratenumeric")
    (LaTeX-add-bibitems
     "Nason:2004rx"
     "Frixione:2007vw"
     "Alioli:2010xd"
     "Frixione:2002ik"
     "Frixione:2006gn"
     "Dawson:1990zj"
     "Djouadi:1991tka"
     "Spira:1995rr"
     "Alioli:2008tz"
     "mcfm"
     "Cacciari:2005hq"
     "Boos:2001cv"
     "Alwall:2006yp"
     "Altarelli:1989wu"
     "Whalley:2005nh"
     "Alioli:2008gx"
     "Bagnaschi:2011tu")
    (TeX-add-symbols
     '("tmtexttt" 1)
     '("tmtextit" 1)
     "sss"
     "as"
     "POWHEG"
     "POWHEGBOX"
     "HERWIG"
     "PYTHIA"
     "MCatNLO"
     "kt"
     "pt"
     "LambdaQCD")
    (TeX-run-style-hooks
     "url"
     "enumerate"
     "amssymb"
     "amsmath"
     ""
     "latex2e"
     "JHEP310"
     "JHEP3"
     "paper")))

