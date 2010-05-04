(TeX-add-style-hook "unmarked"
 (lambda ()
    (LaTeX-add-bibliographies
     "/home/ian/Documents/bibtex/dissertation")
    (LaTeX-add-labels
     "sec:models-impl-unmark"
     "tab:models"
     "sec:occ"
     "sec:repeated-count-data"
     "eq:pc2"
     "eq:mp2"
     "sec:double-observ-sampl"
     "sec:unmarked-usage"
     "sec:data-requirements"
     "sec:fitting-models"
     "sec:examining-model-fits"
     "fig:preddet"
     "fig:pb"
     "sec:future-direct-unmark")
    (TeX-add-symbols
     "um"
     "rlang"
     "scovs"
     "ocovs")
    (TeX-run-style-hooks
     "rotating"
     "inputenc"
     "utf8"
     "amssymb"
     "amsmath"
     "latex2e"
     "jss10"
     "jss"
     "article"
     "shortnames")))

