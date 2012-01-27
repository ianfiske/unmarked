(TeX-add-style-hook "cap-recap"
 (lambda ()
    (LaTeX-add-bibliographies
     "unmarked")
    (TeX-run-style-hooks
     "fullpage"
     "natbib"
     "Sweave"
     "fontenc"
     "OT1"
     "latex2e"
     "art10"
     "article"
     "a4paper")))

