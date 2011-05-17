(TeX-add-style-hook "unmarked"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "tab:models")
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

