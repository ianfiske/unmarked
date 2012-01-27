(TeX-add-style-hook "cap-recap"
 (lambda ()
    (LaTeX-add-bibliographies
     "unmarked")
    (LaTeX-add-labels
     "mod")
    (TeX-run-style-hooks
     "amsmath"
     "fullpage"
     "natbib"
     "Sweave"
     "fontenc"
     "OT1"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

