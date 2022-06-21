(setq org-latex-pdf-process 
      '("latexmk -pdflatex='pdflatex -interaction nonstopmode' -bibtex -output-directory=export/ -pdf -bibtex -f %f"))
